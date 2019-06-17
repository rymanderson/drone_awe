import csv
import os

class Drone:
	# 'Class used to store key drone characteristics'

	# class variables go here:

	# name        = "dji-Mavic2" # must match the .param file in the params/ directory
    # wingtype    = "rotary"
	# battery     = {'type':'LiPo', 'capacity':0.0, 'voltage':0.0}

	name = "dji-Mavic2"
	params = {
			'wingtype':'rotary', # some of these are default values, but otherwise the key is empty
			'TOW':None,
			'max_speed':None,
			'cruise_speed':None,
			'max_alt':None,
			'max_t':None,
			'max_t_hover':None,
			'min_temp':0,
			'max_temp':40,
			'power_rating':None,
			'batt_type':'LiPo',
			'batt_capacity':None,
			'batt_voltage':None,
			'batt_cells':None,
			'batt_energy':None,
			'batt_mass':None,
			'waterprrof':'no',
			'batt_rechargetime':None,
			'max_payload':None,
			'rotor_diam':0.2
			}

	# methods go here:
	def __init__(self, name):
		self.name   = name

	def getparams(self):
		# find name.param in the params/ directory
		with open(os.path.join('./params', 'paramlist.param'),'r') as paramlist:
			with open(os.path.join('./params', self.name + ".param"),'r') as paramfile:
				paramlist = paramlist.readlines()
				paramfile = paramfile.readlines() #open file with read privilege

				specs = [] #initialize lists
				values = []

				for line in paramfile: #separates columns from .param file
					parts = line.split()
					specs.append(parts[0])
					values.append(parts[1])
				for line in paramlist:
					line = line.strip() #gets rid of quotation marks when comparing strings
					for spec in specs:
						if line == spec:
							if spec == "wingtype" or spec == "batt_type" or spec == 'waterproof' or spec == 'VTOL':
								self.params[spec] = values[specs.index(spec)]
							else:
								self.params[spec] = float(values[specs.index(spec)])

				# for param in paramfile: #separates columns from .param file
				# 	parts = param.split()
				# 	spec = parts[0]
				# 	value = parts[1]
				# 	specs.append(spec)
				# 	values.append(value)
				# 	for line in paramlist:
				# 		line = line.strip() #gets rid of quotation marks when comparing strings
				# 		if line == spec:
				# 			if spec == "wingtype" or spec == "batt_type":
				# 				self.params[spec] = value
				# 			else:
				# 				self.params[spec] = float(value)

print("Successfully imported `Drone.py`")
