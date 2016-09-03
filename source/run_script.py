from algorithm import core_algorithm
from search_space import _splitBox,Search_Space,_matlabRectangleFormat
import search_space
from multiprocessing import Process,cpu_count,Manager
from subprocess import call
import math,sys,os
import numpy as np
import json
np.set_printoptions(threshold=np.nan)

max_number_of_points = 565
timing_switch = False

file_path = os.path.realpath(__file__)
template_stem = os.path.abspath(sys.argv[1])
parameters_path = os.path.join(template_stem,"parameters.json")
with open(parameters_path,"r") as f: 
	parameters = json.load(f)


if 'skip_cheap' not in parameters: 
	skip_cheap_segment = False
else: 
	skip_cheap_segment = True

overall_bounds = parameters["bounds"]
overall_bounds = [float(entry) for entry in overall_bounds]
dimensionality = len(overall_bounds)/2
local_template_location = os.path.join(template_stem,"cheap")
global_template_location = os.path.join(template_stem,"minimizer")
eval_template_location = os.path.join(template_stem,"evaluate")
mpi_executable_location = os.path.realpath(os.path.join(file_path,"../../mpich-install/bin/mpirun"))
bertini_string = os.path.realpath(os.path.join(file_path,"../../BertiniLinuxMPI2_v1.5/bertini"))

# This variable gets redefined by subsequent functions
output_path=""

pointcloud_to_distmat_path = os.path.abspath(os.path.join(file_path,"../../PH-roadmap/matlab/"))
dipha_script_path = os.path.realpath(os.path.join(file_path,"../../dipha/matlab/"))
matlab_path = "matlab"
dipha_executable_path = os.path.realpath(os.path.join(file_path,"../../dipha/dipha"))

number_of_min_processors = cpu_count()
number_of_cheap_processors = cpu_count()/3
number_of_functions = parameters["number_of_functions"]
density_parameter = float(parameters["density_parameter"])
density_parameter_name = ""
rounding_precision = 0.0
old_density_parameter = 0 

r_out = list() 
p_out = list()
pr_out = list() 

def single_run(): 
	global r_out,p_out,pr_out,output_path,density_parameter_name,density_parameter,old_density_parameter

	# Auto formatting for floats does not work for file naming; 
	# This fixes that situation; slightly hacky
	density_parameter_mod = int(math.log10(density_parameter))
	density_parameter_start = int(density_parameter*math.pow(10,(2-density_parameter_mod)))
	if density_parameter_mod-2 < 0: 
		density_parameter_name = str(density_parameter_start)+"e"+str(density_parameter_mod-2)
	else: 
		density_parameter_name = str(density_parameter_start)+"e+"+str(density_parameter_mod-2)
	density_parameter_name += "_2"
	sample_name = density_parameter_name + '_sample.txt'
	output_path = os.path.join(template_stem,sample_name)

	bounds = overall_bounds
	manager = Manager()
	points_queue = manager.Queue()
	rectangle_queue = manager.Queue()
	sample_rectangles_queue = manager.Queue() 


	if old_density_parameter != 0 and len(p_out) > 0: 
		core_algorithm(local_template_location,global_template_location,mpi_executable_location,bertini_string,number_of_min_processors,number_of_cheap_processors,number_of_functions,np.array(bounds),overall_bounds,density_parameter,rounding_precision,points_queue,rectangle_queue,sample_rectangles_queue,False,skip_cheap_segment,(p_out,old_density_parameter) )
	else: 
		core_algorithm(local_template_location,global_template_location,mpi_executable_location,bertini_string,number_of_min_processors,number_of_cheap_processors,number_of_functions,np.array(bounds),overall_bounds,density_parameter,rounding_precision,points_queue,rectangle_queue,sample_rectangles_queue,False,skip_cheap_segment)

	rectangle_queue.put('END')
	points_queue.put('END')
	sample_rectangles_queue.put('END')



	for points in iter(points_queue.get,'END'): 
		p_out += points
	for rectangles in iter(rectangle_queue.get,'END'):
		r_out += rectangles
	for rectangles in iter(sample_rectangles_queue.get,'END'): 
		pr_out += rectangles

	p_out = np.array(p_out)
	np.savetxt(output_path,p_out,delimiter=",")


def run_dipha(): 
	global output_path
	sample_name = density_parameter_name + '_sample.txt'
	output_path = os.path.join(template_stem,sample_name)

	with open("/dev/null","w") as devnull:
		# Save the distance matrix for dipha 
		dipha_input_stem = '%s_dipha_input' % density_parameter_name
		dipha_output_stem = '%s_dipha_output' % density_parameter_name
		dipha_input_path = os.path.join(template_stem,dipha_input_stem)
		dipha_output_path = os.path.join(template_stem,dipha_output_stem)
		call(args=[matlab_path,"-nodisplay -nosplash -nodesktop -r","\"cd %s;X=pointcloud_to_distmat('%s');cd %s;save_distance_matrix(X,'%s');exit;\"" % (pointcloud_to_distmat_path,output_path,dipha_script_path,dipha_input_path)],
			stdout=devnull.fileno(),cwd=template_stem)
	
		# Call dipha
		dipha_call_string = dipha_executable_path+" --dual --upper_dim 3 "+dipha_input_stem + " " + dipha_output_stem
		print "dipha info: ",dipha_call_string,template_stem
		call(args=dipha_call_string,cwd=template_stem,shell=True)

		# Save persistence diagram image 
		diagram_output_path = os.path.join(template_stem,'%s_diagram.jpeg' % density_parameter_name)
		call(args=[matlab_path,"-nodisplay -nosplash -nodesktop -r","\"cd %s;plot_persistence_diagram('%s');saveas(gcf,'%s');exit;\"" % (dipha_script_path,dipha_output_path,diagram_output_path)],
			stdout=devnull.fileno(),cwd=template_stem)


if timing_switch == True: 
	r_out = list() 
	p_out = list()
	pr_out = list() 
	output_path = os.path.join(template_stem,('%s' % density_parameter_name)+"_sample.txt")
	single_run()
else:
	while len(p_out) < max_number_of_points:
		r_out = list() 
		p_out = list()
		pr_out = list() 
		output_path = os.path.join(template_stem,('%s' % density_parameter_name)+"_sample.txt")
		single_run()
		if len(p_out) < max_number_of_points:
			run_dipha()
		old_density_parameter = density_parameter
		density_parameter = 0.90 * density_parameter
