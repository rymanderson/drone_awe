#test file for Plotter.py

# import datetime
# import matplotlib.pyplot as plt
import Plotter
# from Simulation.Plotter import Plotter

x = [0,1,2,3,4,5,6,7,8,9,10]
y = [100,95,89,82,74,65,55,44,32,19,5]
xlabel = "payload"
ylabel = "range"
title = "Test Title"

p = Plotter.Plotter(x,xlabel,y,ylabel,title)
p.plot_line()
p.plot_scatter()
