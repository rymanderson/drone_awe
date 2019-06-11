import csv

class Drone:
	# 'Class used to store key drone characteristics'

	# class variables go here:

	# name        = "Mavic2" # must match the .param file in the params/ directory
    # wingtype    = "rotary"
	# battery     = {'type':'LiPo', 'capacity':0.0, 'voltage':0.0}

	params = {
			name:"dji-Mavic2",
			wingtype:'rotary',
			battery:{'type':'LiPo',
					'capacity':3850,
					'voltage':15.4,
					'cells':4,
					'energy':59.29,
					'mass':0.297
					},
			TOW: 0.907,
			max_speed:20,
			max_alt:6000,
			max_t:31,
			max_t_hover:29,
			max_tilt:35,
			min_temp:-10,
			max_temp:40,
			power_rating:60
			}

	# methods go here:
	def __init__(self, name):
		self.name   = name

	def getparams(self):
		# find name.param in the params/ directory
		from self.name import params
        self.params[] = params[]

print("Successfully imported `Drone.py`")
