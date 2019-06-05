class Drone:
	'Class used to store key drone characteristics'

	# class variables go here:
	name        = "Mavic2"                      # must match the .param file in the params/ directory
    wingtype    = "rotary"
	battery     = {'type':'LiPo', 'capacity':5}

	# methods go here:
	def __init__(self, name):
		self.name   = name

print("Successfully imported `Drone.py`")
