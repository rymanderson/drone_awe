import matplotlib.pyplot as plt

class Plotter:
	# 'Class used to make plots of simulation results'

	# class variables go here:
	xlabel = "Time [s]"
	ylabel = "Range [m]"
	title = "Example Drone"

	# methods go here:
	def __init__(self, axistitle,xlabel,ylabel):
		self.title = axistitle
		self.xlabel = xlabel
		self.ylabel = ylabel

	def plot(self,x,y)
		plt.plot(x,y)
		plt.xlabel(self.xlabel)
		plt.ylabel(self.ylabel)
		plt.title(self.title)
		plt.show()

print("Successfully imported `Plotter.py`")
