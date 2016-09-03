import os,sys
sys.path.append(os.path.abspath(os.path.abspath('.')))
import search_space

sb = search_space._splitBox

bounds = [[-1.0,3.0,-2.0,3.0]]

def mrf(box): 
	rstr = "patch([%1.8f,%1.8f,%1.8f,%1.8f],[%1.8f,%1.8f,%1.8f,%1.8f],\'red\')" %(box[0],box[1],box[1],box[0],box[2],box[2],box[3],box[3]) 
	rstr += "\n hold on"
	return rstr

for i in xrange(0,5):
	new_bounds = list()
	for box in bounds: 
		new_bounds += sb(box)
	bounds = new_bounds
bounds = [mrf(box) for box in bounds]

for box_string in bounds: 
	print box_string