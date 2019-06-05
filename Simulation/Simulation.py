class Simulation:
        'Class used to run simulations'

        # class variables go here:
        timestep = 1.0 # in seconds
	clock 	 = 0.0 # tracks the current time
	counter  = 0   # tracks the iteration number (0-indexed)

        # methods go here:
        def __init__(self, timestep):
                self.timestep = timestep

print("Successfully imported `Simulation.py`")

