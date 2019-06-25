import matplotlib.pyplot as plt
import numpy as np
import datetime
import csv
import os
import classes
import functions as fun

simulationparams    = fun.getParams('Simulation','settings_list.txt','settings.txt'," ")

validation = None
if simulationparams['validation'] == True: #validation == True
    validation = True
    validationcase = simulationparams['validationcase']
    simulationparams = fun.getParams('Validation/' + validationcase,'settings_list.txt','settings.txt'," ","params/Simulation") #specifies settings_list is in separate path
else:
    validation == False

xlabel               = simulationparams['xlabel']
# ensure xlabel is an independent variable
independentvariables = [
                        "startstateofcharge",
                        "altitude",
                        "temperaturesealevel",
                        "dropsize",
                        "liquidwatercontent",
                        "newtemperature",
                        "windspeed",
                        "winddirection",
                        "relativehumidity",
                        "payload",
                        "missionspeed",
                        "model"
                        ]

if xlabel not in independentvariables:
    raise(Exception("~~~~~ ERROR: desired x variable is not independent ~~~~~"))

ylabel               = simulationparams['ylabel']

weatherlist         = []
    
#instantiate drone
if simulationparams['drone'] == True:
    conversions = fun.getParams('Drone','paramlist.param','conversions.param'," ")
    dronename           = simulationparams['dronename']
    if validation:
        droneparams         = fun.getParams('Validation/' + validationcase,'paramlist.param',dronename + ".param"," ","params/Drone")
    else:
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

weatherparams   = []
for weathertype in weatherlist:
    weatherparams = weatherparams + weathertype.params

weather         = classes.Weather(simulationparams['altitude'],simulationparams['temperaturesealevel'])
power           = classes.Power(drone,weather)

#simulation variables
timestep        = simulationparams['timestep'] # more relevant later
simulationtype  = simulationparams['simulationtype']
desiredresult   = simulationparams['ylabel']
xbegin          = simulationparams['xbegin']
xend            = simulationparams['xend']
numsimulations  = simulationparams['xnumber']

print("EXE.py:      Independent variable is ",xlabel)
print("EXE.py:      Desired Result is       ",desiredresult)

simulation      = classes.Simulation(timestep,simulationtype)#,desiredresult)
x               = np.linspace(xbegin,xend,numsimulations)
y               = []

for xvalue in x:
    # update value
    ## determine x location
    if xlabel in drone.params:
        drone.params[xlabel] = xvalue
        power.update(drone,weather,simulationparams['model'],mission)
        battery.update()
    elif xlabel in weather.params:
        weather.params[xlabel] = xvalue
        weather.update()
        power.update(drone,weather,simulationparams['model'],mission)
        battery.update()
    elif xlabel in mission.params:
        mission.params[xlabel] = xvalue
        power.update(drone,weather,simulationparams['model'],mission)
        battery.update()
    elif xlabel in simulationparams:
        simulationparams[xlabel] = xvalue
        power.update(drone,weather,simulationparams['model'],mission)
        battery.update()
    else:
        raise(Exception("~~~~~ ERROR: desired x variable not set ~~~~~"))
    
    simulation.run(drone,battery,power,weather,mission)

    if ylabel in drone.params:
        y.append(drone.params[ylabel])
    elif ylabel in simulation.params:
        y.append(simulation.params[ylabel])
    elif ylabel in weather.params:
        y.append(weather.params[ylabel])
    elif ylabel in mission.params:
        y.append(mission.params[ylabel])
    elif ylabel in power.params:
        y.append(power.params[ylabel])
    elif ylabel in simulationparams:
        y.append(simulationparams[ylabel])
    else:
        raise(Exception("~~~~~ ERROR: desired y variable not found ~~~~~"))

print("x data includes:    ",x)
print("y data includes:    ",y)

if not validation: #proceed with normal plot

    if simulationparams['plot'] == True:
        xlabel = simulationparams['xlabel']
        ylabel = desiredresult
        axistitle = simulationparams['title']
        plotter = classes.Plotter(x,xlabel,y,ylabel,axistitle)
        plotter.plot_line()
    else: 
        print('No plot functionality has been defined.')

else: # Plot validation data on top of our model
    xvalid,yvalid = fun.getXandY(validationcase,",")
    # yvalid = [x * 60.0 for x in yvalid] #only for converting from minutes to seconds until we get the conversion working before plotting

    if simulationparams['plot'] == True:
        xlabel = simulationparams['xlabel']
        ylabel = desiredresult
        axistitle = simulationparams['title'] + " Validation"
        plotter = classes.Plotter(x,xlabel,y,ylabel,axistitle)
        plotter.plot_validation(xvalid,yvalid)
    else: 
        print('No plot functionality has been defined.')
