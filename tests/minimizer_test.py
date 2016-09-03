from walking import bertiniMinimizer
import numpy as np

path_to_template = "/users/mf15pe/Documents/sampling_algorithm/BertiniExample"
input_point = np.array([0.2,0.5])
number_of_bertini_processes = 8
number_of_functions = 1
dimensionality = 2
bertini_executable = "bertini"
mpi_executable = "/users/mf15pe/Documents/mpich-install/bin/mpirun"
space_bounds = np.array([[-5,5],[-5,5]])


print bertiniMinimizer(path_to_template,input_point,number_of_bertini_processes,number_of_functions,dimensionality,bertini_executable,mpi_executable,space_bounds)
