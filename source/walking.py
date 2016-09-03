import numpy as np 
import tempfile
import os
import sys
from os.path import join
from shutil import copyfile,rmtree
from subprocess import call
from multiprocessing import Process, Queue



norm = np.linalg.norm


def readSolutionsFile(solutions_path,dimensionality):
	solution_points = np.empty([0,dimensionality])
	
	if not os.path.isfile(solutions_path): 
		return []
	# Read in solution points
	solutions_file = open(solutions_path,"r") 
	lines = solutions_file.readlines() 
	solutions_file.close()

	number_of_solutions = int(lines[0])


	# This shouldn't happen, but sometimes thresholds don't work
	if number_of_solutions == 0:
		return solution_points;

	# Removes unnecessary lines at the beginning of the file
	lines = lines[2:]
	coordinate = 0 
	current_point = np.empty(dimensionality)
	line_number = 0

	while line_number < len(lines): 
		# This triggers if we've read in all the non-auxiliary variables
		if len(solution_points) == number_of_solutions: 
			break 
		if coordinate == dimensionality:
			coordinate = 0 
			solution_points = np.vstack([solution_points,current_point])
			# Lines with auxiliary variables
			counter = 0 
			while (lines[line_number+counter]!='\n'):
				counter += 1
			line_number += counter+1
			continue

		solution_fragment = lines[line_number]
		solution_fragment = solution_fragment.split(" ")[0]
		solution_fragment = float(solution_fragment)
		current_point[coordinate] = solution_fragment
		coordinate += 1
		line_number += 1
	return solution_points

def bertiniMinimizer(run_path,input_point,number_of_bertini_processes,number_of_functions,dimensionality,bertini_executable,mpi_executable,space_bounds):
	# final_parameters contains only information about the test point
	# only need one
	# One parameter per variable, plus one if we're doing the smoothing thing
	def constructFinalParametersString(final_parameters,parameters_are_complex,number_of_functions):
		output_string = str((dimensionality+number_of_functions)) + "\n\n"
		for parameter in final_parameters: 
			if parameters_are_complex:
				output_string += str(parameter[0]) + " " + str(parameter[1]) + "\n"
			else: 
				output_string += str(parameter) + " " + "0\n"
		for dummy_variable in xrange(number_of_functions):
			output_string += "0 0\n"
		return output_string

	final_parameters_path = join(run_path,"final_parameters")

	solutions_file_path = join(run_path,"real_finite_solutions")

	devnull = open("/dev/null","w")

	final_parameters_file = open(final_parameters_path,"w")
	final_parameters_file.write(constructFinalParametersString(input_point,False,number_of_functions))
	final_parameters_file.close()
	call(args=[mpi_executable,"-np", str(number_of_bertini_processes), bertini_executable,"input_param","start_points"],stdout=devnull.fileno(),cwd=run_path)
	solution_points = readSolutionsFile(solutions_file_path,dimensionality)
	if len(solution_points) == 0: 
		return {"points": False, "distance": 0.0}
	solution_distances = np.array([norm(input_point-critical_point) for critical_point in solution_points])
	min_index = np.argmin(solution_distances)
	devnull.close()
	return {"points": solution_points, "distance": solution_distances[min_index],"min_point": solution_points[min_index]}
 

def bertiniEval(run_path,points,number_of_bertini_processes,number_of_functions,dimensionality,bertini_executable,mpi_executable):
	start_path = join(run_path,"start")
	solutions_file_path = join(run_path,"function")

	output_string = str((len(points))) + "\n\n"
	for point in points:
		for coordinate in point:  
			output_string += str(coordinate) + " " + "0\n"

	devnull = open("/dev/null","w")
	start_file = open(start_path,"w")
	start_file.write(output_string)
	start_file.close()	
	call(args=[mpi_executable,"-np", str(number_of_bertini_processes), bertini_executable],stdout=devnull.fileno(),cwd=run_path)
	devnull.close()
	# The second argument here is number of functions instead of the dimensionality
	# since evaluating number_of_functions polynomials on a single point produces that 
	# many answers
	evaluated_points = readSolutionsFile(solutions_file_path,number_of_functions)
	return evaluated_points

def bertiniWalk(path_to_template,number_of_bertini_processes,number_of_functions,number_of_points,dimensionality,walk_length,bertini_executable,mpi_executable,space_bounds,return_queue):	
	run_path = tempfile.mkdtemp()

	if not os.path.exists(run_path):
		os.makedirs(run_path)


	input_template_path = join(path_to_template,"input_template")
	input_config_path = join(run_path,"CONFIG")
	start_path = join(run_path,"start_points")
	solutions_path = join(run_path,"real_solutions")

	copyfile(input_template_path,input_config_path)


	def constructStartSolutionsString(test_point,number_of_functions,patch_variables): 
		output_string = "1\n"
		for coordinate in test_point: 
			output_string += str(coordinate) + " 0\n"
		output_string += str(patch_variables[0]) + " 0\n"
		for variable in patch_variables[1:]: 
			output_string += "0 0\n"
		return output_string

	# Set up the seed point for the walk
	seed_point = np.empty(dimensionality)
	for j in xrange(0,dimensionality):
		index = 2*j 
		seed_point[j] = np.random.uniform(low=space_bounds[index],high=space_bounds[index+1])

	test_point = seed_point

	devnull = open("/dev/null","w")

	output_points = np.empty([0,dimensionality])
	while len(output_points) < number_of_points:
		print "Number of points: ", len(output_points)
		# Set up bertini input files and run bertini
		start = open(start_path,"w")
		patch_variables = np.random.uniform(low=-1e1,high=1e1,size=number_of_functions+1)
		start.write(constructStartSolutionsString(test_point,number_of_functions,patch_variables))
		start.close()
		config = open(input_config_path,"w")
		old_data_file = open(input_template_path,"r")
		config_data = old_data_file.read()
		old_data_file.close()
		for coordinate in xrange(0,len(test_point)): 
			config_data = config_data.replace("CHANGE_ME_"+str(coordinate+1),str(test_point[coordinate]))
		for coordinate in xrange(len(patch_variables)): 
			config_data = config_data.replace("CHANGE_ME_A"+str(coordinate+1),str(patch_variables[coordinate]))

		config.write(config_data)
		config.close()

		call(args=[mpi_executable,"-np", str(number_of_bertini_processes), bertini_executable,"CONFIG","start_points"],cwd=run_path,stdout=devnull.fileno())
		
		solution_points = readSolutionsFile(solutions_path,dimensionality)
		
		# If for some reason the path wasn't trackable;
		if len(solution_points)==0:
			continue
		# Add solution points to the output points
		output_points = np.vstack([output_points,solution_points])
		current_point = solution_points[0]
		
		# Calculate next test test_point
		
		# Find unit normal vector
		normal_vector = test_point-current_point
		normal_vector = normal_vector/norm(normal_vector)


		# Pick (non-parallel to normal_vector) random unit vector     
		rand_vector = np.random.normal(size=dimensionality)
		rand_vector = rand_vector/np.linalg.norm(rand_vector)

		while (rand_vector == normal_vector).all(): 
			rand_vector = np.random.normal(size=dimensionality)
			rand_vector = rand_vector/np.linalg.norm(rand_vector)
		# Project it onto the hyperplane 
		rand_vector = rand_vector - np.dot(rand_vector,normal_vector)*normal_vector
		test_point = walk_length*(rand_vector/norm(rand_vector)) + current_point
	rmtree(run_path)
	devnull.close()
	return_queue.put(output_points)
	return 



