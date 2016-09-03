import os
import sys

file_name = os.path.abspath(sys.argv[1])

output = ""
with open(file_name) as file: 
	content = file.readlines()
	for line in content: 
		for t in line.split():
			l = []
			if t == "" or len(t)==1: 
				continue
			if t[-1]=="]": 
				t = t[:-2]
			if t[0]=="[": 
				t = t[1:]
			try:
				l.append(float(t))
			except ValueError:
				pass
			for number in l: 
				output += str(number) + ","
		output = output[0:-1]
		output += "\n"

print output
