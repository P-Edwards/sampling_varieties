import numpy as np
from search_space import Search_Space,_splitBox,search_box,_middleOfBox
from walking import bertiniWalk,bertiniMinimizer,bertiniEval
import multiprocessing
import os,sys
import resource
from shutil import rmtree,copyfile
import tempfile
from os.path import join
from math import sqrt

Process = multiprocessing.Process
ss = Search_Space

ALLOWED_MAX = 550
number_of_bad_boxes_limit = 1e6
number_of_eval_functions = 3


def check_for_nonzero_entry(input_list): 
	for entry in input_list:
		if entry > 0: 
			return True
	return False



def core_algorithm(local_template_location,global_template_location,mpi_executable_location,bertini_string,number_of_min_processors,number_of_cheap_processors,number_of_functions,bounds,global_bounds,density_parameter,rounding_precision,points_queue,rectangles_queue,sample_rectangles_queue,eval_template_location=False,skip_cheap=False,smaller_box_info=False):
	if smaller_box_info: 	
		space = ss(density_parameter,rounding_precision,bounds,global_bounds,smaller_box_info[0],smaller_box_info[1])
	else: 
		space = ss(density_parameter,rounding_precision,bounds,global_bounds)
	dimensionality = space.dimension

	# Step 1: Gather the local information and add it to the space
	process_list = list([])
	manager = multiprocessing.Manager() 
	return_queue = manager.Queue()
	
	heuristic_switch = True


	side_length = 2*density_parameter/sqrt(dimensionality)

	grid_points = list()
	print "side length", side_length
	for i in xrange(0,len(bounds),2): 
		dim = i/2
		start_point = bounds[i]
		end_point = bounds[i+1]
		# Add list for this dimension
		grid_points += [list([])]
		number_of_grid_points = int( (end_point-start_point)/side_length )
		if number_of_grid_points<=1: 
			grid_points[dim] = [start_point,(start_point+end_point)/2.0,end_point]
		else: 
			for j in xrange(number_of_grid_points+1): 
				grid_points[dim] += [start_point+float(j)*side_length]
			grid_points[dim]+= [end_point]
	new_boxes = _splitBox(bounds,grid_points)
	new_boxes = [search_box(box) for box in new_boxes]

	max_number_of_points = len(new_boxes)

	simple_run_switch = False
	if max_number_of_points <= ALLOWED_MAX: 
		space.bad_boxes = new_boxes
		heuristic_switch = False
		simple_run_switch = True

	if skip_cheap!=True and heuristic_switch==True:	
		number_of_walk_points = 0.08*max_number_of_points
		number_of_walk_points = max(number_of_walk_points/float(number_of_cheap_processors),1)
		print "Number of walk points: ", number_of_walk_points
		for dummy_variable in xrange(0,number_of_cheap_processors): 
			#(path_to_template,number_of_bertini_processes,number_of_functions,number_of_points,dimensionality,walk_length,bertini_executable,mpi_executable,space_bounds,return_queue)
			arguments = (local_template_location,1,1,number_of_walk_points,dimensionality,1.5*density_parameter,bertini_string,mpi_executable_location,bounds,return_queue)
			process = Process(target=bertiniWalk,args=arguments)
			process.start()
			process_list.append(process)

		for process in process_list: 
			process.join()

		results = return_queue.get()

		for point in results: 
			space.addPoint(point,False,True)


	run_path = tempfile.mkdtemp()
	if not os.path.exists(run_path):
		os.makedirs(run_path)

	if eval_template_location != False: 
		eval_path = join(run_path,"eval")
		os.makedirs(eval_path)
		eval_configuration = join(eval_template_location,"input")
		copyfile(eval_configuration,join(eval_path,"input"))


	configuration_file = join(global_template_location,"input_param")
	input_solutions = join(global_template_location,"start_points")
	start_parameters = join(global_template_location,"start_parameters")



	copyfile(configuration_file,join(run_path,"input_param"))
	copyfile(input_solutions,join(run_path,"start_points"))
	copyfile(start_parameters,join(run_path,"start_parameters"))


	test_point = space.checkCover()

	# Only start adding sample points once the exclusion boxes have gotten uselessly small
	# in order to reduce number of samples

	print "Test point: ", test_point
	if heuristic_switch==True and simple_run_switch==False: 
		while test_point!=True and space.current_max_length > space.max_length_limit and len(space.bad_boxes) < number_of_bad_boxes_limit:
			min_distance_results =  bertiniMinimizer(run_path,test_point,number_of_min_processors,number_of_functions,dimensionality,bertini_string,mpi_executable_location,bounds)
			
			# It's highly unlikely but technically possible to hit a singular test point 
			# value. This handles that case.
			if min_distance_results['distance'] == 0: 
				# Pick random unit vector     
				rand_vector = np.random.normal(size=dimensionality)
				rand_vector = rand_vector/np.linalg.norm(rand_vector)
				test_point = np.array(test_point) + rand_vector
				test_point = list(test_point)
				print "Test point (skipped): ", test_point
				continue

			print "Test point: ", test_point, "Distance: ", min_distance_results['distance']
			print "Max length: ", space.current_max_length, "Limit: ", space.max_length_limit
			# It's possible that the test box is simply a very small box very close to the variety.
			# Adding the solution box for the minimizing point is more efficient in this case than 
			# the smaller exclusion box
			if min_distance_results['distance'] < density_parameter: 
				space.addPoint(min_distance_results['min_point'],False,True)
			space.addPoint(test_point,min_distance_results['distance'])

			test_point = space.checkCover()
	

	if simple_run_switch==False: 
		while test_point!=True:
			min_distance_results =  bertiniMinimizer(run_path,test_point,number_of_min_processors,number_of_functions,dimensionality,bertini_string,mpi_executable_location,bounds)
			
			# It's highly unlikely but technically possible to hit a singular test point 
			# value. This handles that case.
			if min_distance_results['distance'] == 0: 
				# Pick random unit vector     
				rand_vector = np.random.normal(size=dimensionality)
				rand_vector = rand_vector/np.linalg.norm(rand_vector)
				test_point = np.array(test_point) + rand_vector
				test_point = list(test_point)
				print "Test point (skipped): ", test_point
				continue

			space.addPoint(test_point,min_distance_results['distance'])
			# Also add all the critical points that we found

			points = min_distance_results['points']
			if eval_template_location != False: 
				evaluated_points = bertiniEval(eval_path,points,number_of_min_processes,number_of_eval_functions,dimensionality,bertini_string,mpi_executable_location)
				filtered_points = [points[index] for index in xrange(len(points)) if check_for_nonzero_entry(evaluated_points[index])]
				points = filtered_points

			for point in points: 
				if heuristic_switch==True: 
					space.addPoint(point,False,True)
				else: 
					space.addPoint(point)
			test_point = space.checkCover()
			print "Test point(b): ", test_point, "Distance: ", min_distance_results['distance'],"Max length: ", space.current_max_length, "Limit: ", space.max_length_limit

	# In this case we can get a small enough sampling with a much simpler procedure
	else: 
		for box in space.bad_boxes: 
			min_distance_results = bertiniMinimizer(run_path,_middleOfBox(box.box),number_of_min_processors,number_of_functions,dimensionality,bertini_string,mpi_executable_location,bounds)
			if min_distance_results['distance'] < density_parameter: 
				space.addPoint(_middleOfBox(box.box),False,False)
			space.addPoint(_middleOfBox(box.box),min_distance_results['distance'])
		test_point = space.checkCover()
		print "Status",test_point

	rmtree(run_path)

	sample_points = space.outputSamplePoints()
	if len(sample_points) > 0: 
		points_queue.put([list(point) for point in sample_points])

	sample_rectangles = space.outputSampleRectangles()
	sample_rectangles_queue.put(sample_rectangles)
	min_rectangles = space.outputMinPointsMatlab()
	rectangles_queue.put(min_rectangles)


