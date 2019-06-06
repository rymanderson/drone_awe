import csv

class Drone:
	'Class used to store key drone characteristics'

	# class variables go here:
	name        = "Mavic2" # must match the .param file in the params/ directory
    wingtype    = "rotary"
	battery     = {'type':'LiPo', 'capacity':5}

	# methods go here:
	def __init__(self, name):
		self.name   = name

	def getparams(self):
		# find name.param in the params/ directory
		paramlist = open("paramlist.param","r")
		paramfile = open(self.name + ".param", "r") #open file with read privilege

		#read .txt file
		with f as inf:
     		reader = csv.reader(inf, delimiter=" ")
     		specs = list(zip(*reader))[0]

		for param in paramlist
			if "param" in specs:
				#create variable with that name
				#exec("%s = %d" % (param,param_value)) - frowned upon according to Stack Exchange
				#could instead reference dictionaries
				#need to creat list file or dictionary

			# for line in paramfile
			# 	values = line.split()
			# 	print(values[0], " - ", values[1])

		return 0

print("Successfully imported `Drone.py`")
