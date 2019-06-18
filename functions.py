# import required modules
import os

def getParams(classname,listname,filename,delimiter):
	# initialize `params` dictionary
	params = {}
	# find name.param in the params/ directory
	with open(os.path.join('params', classname, listname), 'r') as paramlist:
		with open(os.path.join('params', classname, filename), 'r') as paramfile:
			paramlist = paramlist.readlines()
			paramfile = paramfile.readlines() #open file with read privilege

			specs 		= [] #initialize lists
			values 		= []

			for line in paramfile: #separates columns from .param file
				parts = line.strip().split(delimiter)
				specs.append(parts[0])
				values.append(parts[1])
			for line in paramlist:
				line = line.strip().split(delimiter) #gets rid of quotation marks when comparing strings
				for spec in specs:
					if line[0] == spec:
						if spec == "wingtype" or spec == "batt_type" or spec == 'waterproof' or spec == 'VTOL':
							params[spec] = values[specs.index(spec)]
						else:
							params[spec] = float(values[specs.index(spec)])
			if len(paramlist) > len(params):
				print("~~~~~ WARNING: ", str(len(paramlist)-len(params)), " parameter(s) missing from parameters file ~~~~~")
	return params