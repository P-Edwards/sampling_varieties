import numpy as np 
import tempfile
import os
import sys
from os.path import join
from shutil import copyfile
from subprocess import call

# Arguments start at 1
args = sys.argv 

path_to_template = "/home/parker/Dropbox/Sampling varieties + TDA/BertiniExample" #args[1]
path_to_workplace = tempfile.mkdtemp() #args[2]
number_of_bertini_processes = 8 #int(args[3])
number_of_functions = 1 #int(args[4])
number_of_points = 50 #int(args[5])
dimensionality = 2 #int(args[6])
walk_length = 0.5 #float(args[7])
bertini_executable = "/home/parker/Documents/BertiniLinuxMPI2_v1.5/bertini" #args[8]

space_bounds = np.array([[0,3],[0,3]]) #[] 
# for argument in range(9,len(args)-1,2):
# 	space_bounds.append([args[argument],args[argument+1]]) 
# 	space_bounds = np.array(space_bounds)


run_directory = ""

run_path = join(path_to_workplace,run_directory)
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

# Set up the seed points for the walk
pt = np.empty(dimensionality)
list_of_points = np.empty([number_of_points,dimensionality],dtype=float)
for i in range(0,number_of_points): 
	for j in range (0,dimensionality): 
		pt[j] = np.random.uniform(low=space_bounds[j][0],high=space_bounds[j][1])
	list_of_points[i] = pt

final_parameters_path = join(run_path,"final_parameters")

solutions_file_path = join(run_path,"real_finite_solutions")

for point in list_of_points: 
	final_parameters_file = open(final_parameters_path,"w")
	final_parameters_file.write(constructFinalParametersString(point,False,number_of_functions))
	final_parameters_file.close()
	call(args=["mpiexec","-np", str(number_of_bertini_processes), bertini_executable,"input_param","start_points", "> /dev/null"],cwd=run_path)
	solutions_file = open(solutions_file_path,"r")
	solutions_file.close()
	print point
