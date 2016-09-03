import numpy as np 
import tempfile
import os
import sys
from os.path import join
from shutil import copyfile,rmtree
from subprocess import call
from multiprocessing import Process

norm = np.linalg.norm


def readSolutionsFile(solutions_path):
	solution_points = np.empty([0,dimensionality])
	# Read in solution points
	solutions_file = open(solutions_path,"r") 
	lines = solutions_file.readlines() 
	solutions_file.close()

	number_of_solutions = int(lines[0])


	# This shouldn't happen, but sometimes thresholds don't work
	if number_of_solutions == 0: 
		continue

	test_points = np.vstack([test_points,test_point])

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
			line_number += number_of_functions+1+1
			continue

		solution_fragment = lines[line_number]
		solution_fragment = solution_fragment.split(" ")[0]
		solution_fragment = float(solution_fragment)
		current_point[coordinate] = solution_fragment
		coordinate += 1
		line_number += 1
	return solution_points

def bertiniMinimizer(path_to_template,input_point,number_of_bertini_processes,number_of_functions,dimensionality,bertini_executable,mpi_executable,space_bounds):
	run_path = tempfile.mkdtemp()
	if not os.path.exists(run_path):
		os.makedirs(run_path)

	configuration_file = join(path_to_template,"input_param")
	input_solutions = join(path_to_template,"finite_solutions")
	start_parameters = join(path_to_template,"start_parameters")

	copyfile(configuration_file,join(run_path,"input_param"))
	copyfile(input_solutions,join(run_path,"start_points"))
	copyfile(start_parameters,join(run_path,"start_parameters"))


	# final_parameters contains only information about the test point
	# only need one
	# One parameter per variable, plus one if we're doing the smoothing thing
	def constructFinalParametersString(final_parameters,parameters_are_complex,number_of_functions):
		output_string = str((dimensionality+1)) + "\n\n"
		for parameter in final_parameters: 
			if parameters_are_complex:
				output_string += str(parameter[0]) + " " + str(parameter[1]) + "\n"
			else: 
				output_string += str(parameter) + " " + "0\n"
		output_string += "0 0\n"
		return output_string

	final_parameters_path = join(run_path,"final_parameters")

	solutions_file_path = join(run_path,"real_finite_solutions")

	final_parameters_file = open(final_parameters_path,"w")
	final_parameters_file.write(constructFinalParametersString(input_point,False,number_of_functions))
	final_parameters_file.close()
	call(args=["mpiexec","-np", str(number_of_bertini_processes), bertini_executable,"input_param","start_points", "> /dev/null"],cwd=run_path)
	solution_points = readSolutionsFile(solutions_file_path)
	solution_distances = [norm(input_point-critical_point) for critical_point in solution_points]
	min_index = solution_distances.argmin()
	return solution_points[min_index],solution_distances[min_index]
 


def bertiniWalk(path_to_template,number_of_bertini_processes,number_of_functions,number_of_points,dimensionality,walk_length,bertini_executable,mpi_executable,space_bounds):	
	run_path = tempfile.mkdtemp()

	if not os.path.exists(run_path):
		os.makedirs(run_path)


	input_template_path = join(path_to_template,"input_template")
	input_config_path = join(run_path,"CONFIG")
	start_path = join(run_path,"start_points")
	solutions_path = join(run_path,"real_finite_solutions")

	copyfile(input_template_path,input_config_path)


	# final_parameters contains only information about the test point
	# only need one
	# One parameter per variable, plus one if we're doing the smoothing thing
	def constructFinalParametersString(final_parameters,parameters_are_complex,number_of_functions):
		output_string = str((dimensionality+1)) + "\n\n"
		for parameter in final_parameters: 
			if parameters_are_complex:
				output_string += str(parameter[0]) + " " + str(parameter[1]) + "\n"
			else: 
				output_string += str(parameter) + " " + "0\n"
		output_string += "0.0001 0\n"
		return output_string 

	def constructStartSolutionsString(test_point,number_of_functions): 
		output_string = "1\n"
		for coordinate in test_point: 
			output_string += str(coordinate) + " 0\n"
		for dummy_variable in xrange(0,number_of_functions+1): 
			output_string += "0.0001 0\n"
		return output_string

	# Set up the seed point for the walk
	seed_point = np.empty(dimensionality)
	for j in xrange(0,dimensionality):
		seed_point[j] = np.random.uniform(low=space_bounds[j][0],high=space_bounds[j][1])

	test_point = seed_point

	devnull = open("/dev/null","w")

	test_points = np.empty([0,dimensionality])
	output_points = np.empty([0,dimensionality])
	while len(output_points) < number_of_points:
		# Set up bertini input files and run bertini
		start = open(start_path,"w")
		start.write(constructStartSolutionsString(test_point,number_of_functions))
		start.close()
		config = open(input_config_path,"w")
		old_data_file = open(input_template_path,"r")
		config_data = old_data_file.read()
		old_data_file.close()
		for coordinate in xrange(0,len(test_point)): 
			config_data = config_data.replace("CHANGE_ME_"+str(coordinate+1),str(test_point[coordinate]))
		config.write(config_data)
		config.close()
		call(args=[mpi_executable,"-np", str(number_of_bertini_processes), bertini_executable,"CONFIG","start_points"],cwd=run_path,stdout=devnull.fileno())


		solution_points = readSolutionsFile(solutions_path)
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

	unique_id_modifier = str(int(np.random.uniform()*5000))
		
	# with open("output_x","a") as output_x:
	# 	np.savetxt(output_x,output_points[:,0],newline=",")
	# with open("output_y","a") as output_y:
	# 	np.savetxt(output_y,output_points[:,1],newline=",")
	# with open("test_x","a") as test_x:
	# 	np.savetxt(test_x,test_points[:,0],newline=",")
	# with open("test_y","a") as test_y:
	# 	np.savetxt(test_y,test_points[:,1],newline=",")
	rmtree(run_path)
	return output_points

# process_list = list([]) 
# for dummy_variable in xrange(0,8): 
# 	template_string = "/users/mf15pe/Documents/sampling_algorithm/GDH_example"
# 	executable_string = "/users/mf15pe/Documents/mpich-install/bin/mpirun"
# 	#(path_to_template,path_to_workplace,number_of_bertini_processes,number_of_functions,number_of_points,dimensionality,walk_length,bertini_executable,mpi_executable,space_bounds)
# 	arguments = (template_string,tempfile.mkdtemp(),1,1,1000,2,1.5,"bertini",executable_string,np.array([[-1,3],[-1,3]]))
# 	process = Process(target=bertiniWalk,args=arguments)
# 	process.start()
# 	process_list.append(process)

# for process in process_list: 
# 	process.join()

