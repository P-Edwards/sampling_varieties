from rtree import index
import numpy as np 
from math import sqrt,pow
import sys
import operator,functools

INTERSECTS = -1
MAX_BOUNDS = float(1e6)
ROLLING_AVERAGE_LENGTH = 50
NUMBER_OF_ALLOWED_SKIPS = 49
SKIP_INTERVAL = 0.2

def sum_list(input_list): 
	if len(input_list) == 0: 
		return 0.0
	return np.sum(input_list)

def _doublePoint(point): 
	# if len(point)==2: 
	# 	return (point[0],point[1])
	doubled_point = list()
	for coordinate in point: 
		doubled_point += [coordinate]
		doubled_point += [coordinate]
	return tuple(doubled_point)

def _middleOfBox(box): 
	return_point = list()
	for i in xrange(0,len(box),2):
		midpoint = (box[i]+box[i+1])/2.0
		return_point.append(midpoint)
	return return_point

def _checkRectangles(bounding_rectangle,rectangle_to_check): 
	for coordinate in xrange(0,len(bounding_rectangle),2) :
		check_low = rectangle_to_check[coordinate]
		check_high = rectangle_to_check[coordinate+1]
		bound_low = bounding_rectangle[coordinate]
		bound_high = bounding_rectangle[coordinate+1]
		if check_low < bound_low or check_high > bound_high:		
			return False
	return True


def deinterleave(interleaved):
    dimension = len(interleaved) // 2
    def gen(): 
    	for i in xrange(dimension): 
    		yield interleaved[i]
    		yield interleaved[i+dimension]
    return list(gen())

def _matlabRectangleFormat(box): 
	box = deinterleave(box)
	if len(box)/2 == 2: 
		rstr = "patch([%1.8f,%1.8f,%1.8f,%1.8f],[%1.8f,%1.8f,%1.8f,%1.8f],\'red\')" %(box[0],box[1],box[1],box[0],box[2],box[2],box[3],box[3]) 
		rstr += "\n hold on"
	else: 
		rstr = "cubePlot([%1.8f,%1.8f,%1.8f],%1.8f,%1.8f,%1.8f,\'r\')" %(box[0],box[2],box[4],box[1]-box[0],box[3]-box[2],box[5]-box[4])
	return rstr

def _splitBox(bounds,list_of_grid_points=False):
	if not list_of_grid_points: 
		return [bounds]
	dimension = len(list_of_grid_points)
	if dimension == 1: 
		points = list_of_grid_points[0]
		return [[points[i],points[i+1]] for i in xrange(len(points)-1)]
	else: 
		output = list()
		first_dim_boxes = _splitBox(bounds,[list_of_grid_points[0]])
		second_dim_boxes = _splitBox(bounds,list_of_grid_points[1:])
		for first_box in first_dim_boxes: 
			for second_box in second_dim_boxes: 
				output.append(first_box + second_box)
		return output

def _findIntersectionPoints(first_box,second_box): 
	intersection = list() 
	for coordinate in xrange(0,len(first_box),2): 
		coordinates_to_add = [max(first_box[coordinate],second_box[coordinate]),min(first_box[coordinate+1],second_box[coordinate+1])]
		if coordinates_to_add[0] == coordinates_to_add[1]: 
			return False
		else: 
			intersection += [coordinates_to_add]
	return intersection

def _findSizeOfIntersection(first_box,second_box): 
	intersection = list() 
	for coordinate in xrange(0,len(first_box),2): 
		intersection += [max(first_box[coordinate],second_box[coordinate]),min(first_box[coordinate+1],second_box[coordinate+1])]
	return _maxLengthOfBox(intersection) 

def _maxLengthOfBox(box): 
	box_lengths = np.array([(box[i+1]-box[i]) for i in xrange(0,len(box),2)])
	return np.prod(box_lengths)

# def _maxLengthOfBox(box): 
# 	box_lengths = np.array([(box[i+1]-box[i]) for i in xrange(0,len(box),2)])
# 	# This is a dengerate case we need to check
# 	if np.amin(box_lengths) <= 0: 
# 		return 0.0
# 	return np.amax(box_lengths)

class search_box(object): 
	def __init__(self,box): 
		self.box = box
		self.measure = _maxLengthOfBox(box)


class Search_Space(object): 


	def __init__(self,density_parameter,precision_parameter,space_bounds,global_bounds,points=None,old_distance=None): 
		self.epsilon = density_parameter
		self.delta = precision_parameter
		self.space_bounds = space_bounds
		p = index.Property()
		dimension = len(space_bounds)/2
		self.dimension = dimension
		p.dimension = dimension
		self.tree = index.Index(properties=p)
		self.id_counter = 0
		self.tree.interleaved = False
		self.bad_boxes = list([search_box(space_bounds.flatten())])
		self.global_bounds = global_bounds
		self.max_length_limit = _maxLengthOfBox(self._createRadiusBox([0*i for i in xrange(0,dimension)],self.epsilon))/2.5e3
		self.current_max_length = _maxLengthOfBox(space_bounds)
		self.old_box = list([])
		self.skip_list = list()
		self.skip_radius_multiplier = pow(2.5,dimension)
		if points is not None: 
			self.bad_boxes = [search_box(self._createRadiusBox(point,old_distance)) for point in points]
			for point in points: 
				self.addPoint(point,False,True)

	def _createRadiusBox(self,point,radius): 
		radius = radius - self.delta
		range_length = radius/sqrt(self.dimension)
		box_coordinates = list()
		for coordinate in point: 
			box_coordinates.append(coordinate - range_length)
			box_coordinates.append(coordinate + range_length)
		return box_coordinates
	
	def _skipControl(self,is_skipped): 
		if len(self.skip_list) == ROLLING_AVERAGE_LENGTH: 
			del self.skip_list[0]
		self.skip_list.append(int(is_skipped))

		number_of_skips = sum_list(self.skip_list)
		if number_of_skips > NUMBER_OF_ALLOWED_SKIPS: 
			self.skip_radius_multiplier = max(0,self.skip_radius_multiplier-SKIP_INTERVAL)
			self.skip_list = list()
		return


	# The default behavior (radius=False) is for inserting sample points
	# on the variety; otherwise we're inserting exclusion region boxes	
	def addPoint(self,point,radius=False,skip_on_covered=False):
		label = (point,radius)
		no_input_radius = False
		if radius==False: 
			radius = self.epsilon
			no_input_radius = True
		if radius <= self.delta:
			return
		box = self._createRadiusBox(point,radius)
		if no_input_radius: 
			label=(point,True)
			# Prevents adding sample points outside of the bounds of the 
			# problem
			if _checkRectangles(self.global_bounds,_doublePoint(point)) == False: 
				print "Skipping out of bounds"
				return
			# Prevents adding wholly unnecessary sample points
			if skip_on_covered==True:
				intersecting_objects = list(self.tree.intersection(box,objects=True))
				intersecting_boxes = [item for item in intersecting_objects if item.object[1]==True]
				if len(intersecting_boxes) != 0: 
					intersecting_boxes = [deinterleave(item.bbox) for item in intersecting_boxes]
					box_sizes = [_findSizeOfIntersection(box,bbox) for bbox in intersecting_boxes]
					max_box_size = np.amax(box_sizes)
					if max_box_size*self.skip_radius_multiplier >= _maxLengthOfBox(box): 
						print "Skipping because too close"
						print "Number of skips: ", sum_list(self.skip_list), "Out of", len(self.skip_list)
						self._skipControl(True)
						return 
					if self.checkCover([box]) == True: 
						print "Skipping because covered"
						return
				self._skipControl(False)

		self.tree.insert(self.id_counter,box,label)
		self.id_counter += 1

		modulus = 30 + int(.0001*len(self.bad_boxes))
		if (self.id_counter % modulus == 0): 
			self.bad_boxes.sort(key=lambda box: box.measure)
			self.current_max_length = self.bad_boxes[-1].measure
		return 

	def	cull(self):
		return
		

	def _checkIntersectionStatus(self,bounds):
		# print "-----------------" 
		# print bounds 
		bounds = bounds.box
		intersecting_boxes = list(self.tree.intersection(bounds,objects=True))
		if len(intersecting_boxes)==0: 
			return False
		for item in intersecting_boxes: 
			bbox = deinterleave(item.bbox)
			if _checkRectangles(bbox,bounds)==True:
				return True

		intersecting_boxes = [deinterleave(item.bbox) for item in intersecting_boxes]
		box_sizes = [_findSizeOfIntersection(bounds,bbox) for bbox in intersecting_boxes]
		index = np.argmax(box_sizes)

		# It's possible that the rtree returned a sample/exclusion box as 
		# intersecting the query box because they share exactly a line segment
		# or smaller face in common. This is the first place we can catch and correct for that.
		if box_sizes[index] <= 0.0: 
			return False
		

		box_to_split_along = intersecting_boxes[index]
		current_box = bounds
		grid_points = _findIntersectionPoints(current_box,box_to_split_along)
		



		box_list = list() 
		for index in xrange(len(grid_points)):
			# print "Current box: ", current_box
			current_grid_points = [[current_box[i],current_box[i+1]] for i in xrange(0,len(current_box),2)]
			box_low = current_box[2*index]
			box_high = current_box[2*index+1]

			# It's possible that a given dimension has no grid points, skip it
			# if this happens; the extra list() constructor is for copying
			if len(grid_points[index]) > 0: 
				current_grid_points[index] = list(grid_points[index])

			# These need to be added back for the upcoming splitting
			if current_grid_points[index][0] != box_low: 
				current_grid_points[index] = [box_low] + current_grid_points[index]
			if current_grid_points[index][-1] != box_high:
				current_grid_points[index] += [box_high]


			list_of_split_boxes = _splitBox(current_box,current_grid_points)

			if len(grid_points[index]) > 0: 
				grid_points_to_check = grid_points[index]
			else: 
				grid_points_to_check = [box_low,box_high]
			for box_index in xrange(len(list_of_split_boxes)): 
				coordinates_to_consider = [list_of_split_boxes[box_index][2*index],list_of_split_boxes[box_index][2*index+1]]
				if coordinates_to_consider == grid_points_to_check: 
					current_box = list_of_split_boxes.pop(box_index)
					box_list += list_of_split_boxes
					break;

		box_list += [current_box]
		if len(box_list) == 1: 
			print "Grid points for unsplit box: ", grid_points
		return [search_box(box) for box in box_list]

	# The idea: Split the search space into number_of_processors equal rectangular regions 
	# In parallel, run this recursive search: 
	# (1) Does the current rectangle intersect any rectangles in the tree? If not, return False
	# (2) If so, is the current rectangle strictly contained in any of the tree rectangles? If so, 
	# return true. If not, break the current rectangle into smaller rectangles and check the 
	# smaller rectangles
	# Note: It would be nice to parallelize this in such a way that a signal from a process 
	# that finds a problem will cause the others to gracefully exit
	# Note 2: This search is breadth first- we want to stop searching at the largest 
	# "bad" rectangle, since we otherwise could waste a lot of time confirming a good recatangle
	def checkCover(self,bad_boxes=None):
		bbflag = False
		index_max = 0
		if bad_boxes is None: 
			bad_boxes = self.bad_boxes 
			bbflag = True
			print "Checking cover: ", self.space_bounds,"# of bad boxes: ", len(bad_boxes),self.id_counter
		else: 
			bad_boxes = [search_box(box) for box in bad_boxes]


		while len(bad_boxes)!=0:
			box = bad_boxes.pop() 
			intersection_status = self._checkIntersectionStatus(box) 
			if hasattr(intersection_status,"__iter__"): 
				bad_boxes = intersection_status+bad_boxes
			elif intersection_status==False:
				bad_boxes += [box]
				if bbflag==True: 
					self.bad_boxes = bad_boxes
				return _middleOfBox(box.box)

		# If we reached this point then there are no bad boxes, so 
		# the range is covered
		return True 

	def outputSamplePoints(self): 
		absolute_bounds = list()
		for i in xrange(0,self.dimension): 
			absolute_bounds.append(-MAX_BOUNDS)
			absolute_bounds.append(MAX_BOUNDS)
		all_values = self.tree.intersection(absolute_bounds,objects=True)
		all_values = [item.object[0] for item in all_values if item.object[1]==True]
		return all_values

	def outputSampleRectangles(self): 
		absolute_bounds = list()
		for i in xrange(0,self.dimension): 
			absolute_bounds.append(-MAX_BOUNDS)
			absolute_bounds.append(MAX_BOUNDS)
		all_values = self.tree.intersection(absolute_bounds,objects=True)
		all_values = [_matlabRectangleFormat(item.bbox) for item in all_values if item.object[1]==True]
		return all_values

	def outputMinPoints(self): 
		absolute_bounds = list()
		for i in xrange(0,self.dimension): 
			absolute_bounds.append(-MAX_BOUNDS)
			absolute_bounds.append(MAX_BOUNDS)
		all_values = self.tree.intersection(absolute_bounds,objects=True)
		all_values = [item.object for item in all_values if item.object[1]!=True]
		return all_values

	def outputMinPointsMatlab(self): 
		absolute_bounds = list()
		for i in xrange(0,self.dimension): 
			absolute_bounds.append(-MAX_BOUNDS)
			absolute_bounds.append(MAX_BOUNDS)
		all_values = self.tree.intersection(absolute_bounds,objects=True)
		all_values = [_matlabRectangleFormat(item.bbox) for item in all_values if item.object[1]!=True]
		return all_values

