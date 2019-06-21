import matplotlib.pyplot as plt
import numpy as np
import datetime
import csv
import os
import classes
import functions as fun

simulationparams    = fun.getParams('Simulation','settings_list.txt','settings.txt'," ")
weatherlist         = []
    
#instantiate drone
if simulationparams['drone'] == True:
    conversions = fun.getParams('Drone','paramlist.param','conversions.param'," ")
    dronename           = simulationparams['dronename']
    droneparams         = fun.getParams('Drone','paramlist.param',dronename + '.param',' ')
    droneconversions    = fun.getParams('Drone','paramlist.param','conversions.param',' ')
    drone               = classes.Drone(dronename,droneparams,droneconversions)
else:
    raise Exception('Must specify drone name')

# instantiate battery
stateofhealth       = simulationparams['stateofhealth']
startstateofcharge  = simulationparams['startstateofcharge']
battery             = classes.Battery(drone,stateofhealth,startstateofcharge)

# instantiate mission
missionparams       = fun.getParams('Mission','list.mission','simple.mission'," ")
mission             = classes.Mission(missionparams)

# Lines 29 - 62 can be un-commented later when Weather is ready to test
# #get class initialization info 
# raintest = simulationparams['rain']
# if raintest == True:
#     weatherlist.append('rain')
#     dropsize            = simulationparams['dropsize']
#     liquidwatercontent  = simulationparams['liquidwatercontent']
#     # …
#     rain                = classes.Rain(dropsize,liquidwatercontent)

# temperaturetest = simulationparams['temperature']
# if temperaturetest == True:
#     weatherlist.append('temperature')
#     newtemperature  = simulationparams['newtemperature']
#     temperature     = classes.Temperature(newtemperature)

# windtest = simulationparams['wind']
# if windtest == True:
#     weatherlist.append('wind')
#     speed       = simulationparams['windspeed']
#     direction   = simulationparams['winddirection']
#     wind        = classes.Wind(speed,direction)

# humiditytest = simulationparams['humidity']
# if humiditytest == True:
#     weatherlist.append('humidity')
#     relativehumidity    = simulationparams['relativehumidity']
#     direction           = classes.Humidity(relativehumidity)

# icingtest = simulationparams['icing']
# if icingtest == True:
#     weatherlist.append('icing')
#     icing   = classes.Icing()

# print("Weather parameters are: ")
# print(str(weatherlist)[1:-1]) 

weather         = classes.Weather(simulationparams['altitude'],simulationparams['temperaturesealevel'])
power           = classes.Power(drone,weather)

#simulation variables
timestep        = simulationparams['timestep'] # more relevant later
simulationtype  = simulationparams['simulationtype']
desiredresult   = simulationparams['ylabel']
xbegin          = simulationparams['xbegin']
xend            = simulationparams['xend']
xincrement      = simulationparams['xincrement']
numsimulations  = (xend - xbegin) / xincrement + 1
print("EXE.py:      Desired Result is ",desiredresult)
simulation      = classes.Simulation(timestep,simulationtype,desiredresult)
# x = [z*xincrement for z in range(2*xbegin, 2*xend+1)] #initialize x based on beginng, ending, and increment for payload
x               = np.linspace(xbegin,xend,numsimulations)
y               = []

for payload in x:
    drone.updatePayload(payload)
    drone.updateEfficiencyPropulsive(mission)
    power.updatePower(drone,weather,simulationparams['model'])
    ynext = simulation.run(drone,battery,power,weather)
    y.append(ynext)

if simulationparams['plot'] == True:
    xlabel = simulationparams['xlabel']
    ylabel = desiredresult
    axistitle = simulationparams['title']
    plotter = classes.Plotter(x,xlabel,y,ylabel,axistitle)
    plotter.plot_line()
else: 
    print('No plot functionality has been defined.')
