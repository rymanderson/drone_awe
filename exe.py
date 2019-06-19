import matplotlib.pyplot as plt
import numpy as np
import datetime
import csv
import os
import classes as classes
import functions as fun

simulationparams = fun.getParams('Simulation','settings_list.txt','settings.txt'," ")
print(simulationparams)
weatherlist = []
    
#instantiate drone
if simulationparams['drone'] == True:
    dronename = simulationparams['dronename']
    droneparams = fun.getParams('Drone','paramlist.param',dronename + '.param',' ')
    drone = classes.Drone(dronename,droneparams)
else:
    raise Exception('Must specify drone name')

# instantiate battery
stateofhealth = simulationparams['stateofhealth']
startstateofcharge = simulationparams['startstateofcharge']
battery = classes.Battery(drone,stateofhealth,startstateofcharge)

#get class initialization info
# rainlength = 2 #for now
raintest = simulationparams['rain']
if raintest == True:
    weatherlist.append('rain')
    dropsize = simulationparams['dropsize']
    liquidwatercontent = simulationparams['liquidwatercontent']
    # …
    rain = classes.Rain(dropsize,liquidwatercontent)

temperaturetest = simulationparams['temperature']
if temperaturetest == True:
    weatherlist.append('temperature')
    newtemperature = simulationparams['newtemperature']
    temperature = classes.Temperature(newtemperature)

windtest = simulationparams['wind']
if windtest == True:
    weatherlist.append('wind')
    speed = simulationparams['windspeed']
    direction = simulationparams['winddirection']
    wind = classes.Wind(speed,direction)

humiditytest = simulationparams['humidity']
if humiditytest == True:
    weatherlist.append('humidity')
    relativehumidity = simulationparams['relativehumidity']
    direction = classes.Humidity(relativehumidity)

icingtest = simulationparams['icing']
if icingtest == True:
    weatherlist.append('icing')
    icing = classes.Icing()

print("Weather parameters are: ")
print(str(weatherlist)[1:-1]) 

weather = classes.Weather(weatherlist)

#simulation variables
timestep = simulationparams['timestep']

power = classes.Power(drone,battery,weather)

simulationtype = simulationparams['simulationtype']
desiredresult = simulationparams['ylabel']

simulation = classes.Simulation(timestep,simulationtype,desiredresult)
simulation.run(drone,battery,weather,power)

if simulationparams['plot'] == True:
    xlabel = simulationparams['xlabel']
    ylabel = desiredresult
    axistitle = simulationparams['title']
    plotter = classes.Plotter(simulation.x,xlabel,simulation.y,ylabel,axistitle)
    plotter.plot_line()

# with open('textfile.txt','r') as fp:
#     #instantiate drone
#     drone = readline()
#     if drone == true:
#         dronename = fp.readline()
#         drone = Drone(dronename)
#     else:
#         raise Exception('Must specify drone name')

#     # instantiate battery
#     battery = Battery(drone)

#     #get class initialization info
#         rainlength = 2 #for now
#         rain = fp.readline()
#         if rain == true:
#             Dropsize = fp.readline()
#             WVC = fp.readline()
#             # …
#             rain = Rain(dropsize,WVC,...)
#         else:
#             for i = 1:rainlength:
#                 fp.readline()

#         #simulation variables
#         Timestep = fp.readline()
#         # …
#         #plot settings
#         ##time plots
#         Totalpower	true
#         Rainintensity false
#         …
#         # VVV In text file itself VVV
#         # ##time invariant plots: select two parameters from list:
#         # ### range payload endurance …
#         # #### example:
#         # #### range	payload
#         # plot range payload
#         # ^^^ 		         ^^^
#         while eof == false:
#         nextline = readline()
#         if nextline == plot:
#             plotx.push(fp.readline())
#             ploty.push(fp.readline())
#         else:
#             raise Exception('Invalid plot syntax')

# simulation = Simulation(drone,weather,power,...timestep,...)
# simulation.run()

# plotter = Plotter(inputs….)
# plotter.plot()
