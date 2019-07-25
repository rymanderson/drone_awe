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
    validation = False

xlabel               = simulationparams['xlabel']
# ensure xlabel is an independent variable
independentvariables = [
                        "startstateofcharge",
                        "altitude",
                        "temperature",
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
batterytechnology   = simulationparams['batterytechnology']
battery             = classes.Battery(drone,stateofhealth,startstateofcharge, batterytechnology)

# instantiate mission
missionparams       = fun.getParams('Mission','list.mission','simple.mission'," ")
mission             = classes.Mission(missionparams)

# Temperature
if xlabel == 'temperature':
    newtemperature = simulationparams['xbegin']
else:
    newtemperature  = simulationparams['temperature']
temperatureparams   = {'temperature':newtemperature}        # Ampere-hours
temperature     = classes.Temperature(temperatureparams)
weatherlist.append(temperature)

# Humidity
if xlabel == 'humidity':
    relativehumidity = simulationparams['xbegin']
else:
    relativehumidity = simulationparams['relativehumidity']
humidityparams      = {'relativehumidity':relativehumidity}
humidity            = classes.Humidity(humidityparams)
weatherlist.append(humidity)

# Rain
#     dropsize            = simulationparams['dropsize']
#     liquidwatercontent  = simulationparams['liquidwatercontent']
#     # …
#     rain                = classes.Rain(dropsize,liquidwatercontent)
#     weatherlist.append('rain')

# Wind
#     speed       = simulationparams['windspeed']
#     direction   = simulationparams['winddirection']
#     wind        = classes.Wind(speed,direction)
#     weatherlist.append(wind)

# Icing
#     weatherlist.append('icing')
#     icing   = classes.Icing()

# print("Weather parameters are: ")
# print(str(weatherlist)[1:-1]) 

# weatherparams   = []
# for weathertype in weatherlist:
#     weatherparams = weatherparams + weathertype.params

weather         = classes.Weather(simulationparams['altitude'],weatherlist)
print("Preparing to update weather:")
weather.update()
print("Weather updated.")
power           = classes.Power(drone,weather,mission)

#simulation variables
timestep        = simulationparams['timestep'] # more relevant later
simulationtype  = simulationparams['simulationtype']
desiredresult   = simulationparams['ylabel']
xbegin          = simulationparams['xbegin']
xend            = simulationparams['xend']
numsimulations  = simulationparams['xnumber']

simulation      = classes.Simulation(timestep,simulationtype)#,desiredresult)
x               = np.linspace(xbegin,xend,numsimulations)
y               = []

xplot = x
yplot = []

if "weathereffect" in simulationparams:
    weathereffect = simulationparams["weathereffect"]
    weatherbegin = simulationparams["weatherbegin"]
    weatherend = simulationparams["weatherend"]
    weathernumber = int(simulationparams["weathernumber"])
    wvector = np.linspace(weatherbegin,weatherend,weathernumber)
else:
    weathernumber = int(1)
    wvector = range(weathernumber) # only iterate once

for zvalue in wvector:
    if "weathereffect" in simulationparams:
        if weathereffect == 'temperature':
            print("weathereffect = temperature confirmed")
            weather.weatherlist[0].params["temperature"] = zvalue
        elif weathereffect == 'relativehumidity':
            weather.weatherlist[1].params["relativehummidity"] = zvalue
        else:
            raise(Exception("ERROR: weathereffect not a valid input"))
        weather.update()
        power.update(drone,weather,mission)
        battery.update()
    
    # simulation.run(drone,battery,power,weather,mission)

    # if ylabel in drone.params:
    #     y.append(drone.params[ylabel])
    # elif ylabel in simulation.params:
    #     y.append(simulation.params[ylabel])
    # elif ylabel in weather.params:
    #     y.append(weather.params[ylabel])
    # elif ylabel in mission.params:
    #     y.append(mission.params[ylabel])
    # elif ylabel in power.params:
    #     y.append(power.params[ylabel]*180.0/np.pi)
    # elif ylabel in simulationparams:
    #     y.append(simulationparams[ylabel])
    # else:
    #     raise(Exception("~~~~~ ERROR: desired y variable not found ~~~~~"))

    for xvalue in x:
        # update value
        ## determine x location
        if xlabel in drone.params:
            drone.params[xlabel] = xvalue
            power.update(drone,weather,mission)
            battery.update()
        elif xlabel in weather.params:
            if xlabel == 'temperature':
                weather.weatherlist[0].params[xlabel] = xvalue
            elif xlabel == 'relativehumidity':
                weather.weatherlist[1].params[xlabel] = xvalue
            weather.update()
            power.update(drone,weather,mission)
            battery.update()
        elif xlabel in mission.params:
            mission.params[xlabel] = xvalue
            power.update(drone,weather,mission)
            battery.update()
        elif xlabel in simulationparams:
            simulationparams[xlabel] = xvalue
            power.update(drone,weather,mission)
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

    yplot.append(y)
    y = []

print("x data includes:    ",xplot)
print("y data includes:    ",yplot)

print("")

print("EXE.py:      Independent variable is ",xlabel)
print("EXE.py:      Desired Result is       ",desiredresult)
if "weathereffect" in simulationparams:
    print("EXE.py:      Z iterator is           ",simulationparams['weathereffect'])


if not validation: #proceed with normal plot

    if simulationparams['plot'] == True:
        xlabel = simulationparams['xlabel']
        ylabel = desiredresult
        axistitle = simulationparams['title']
        plotter = classes.Plotter(xplot,xlabel,yplot,ylabel,axistitle,weathernumber)
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
        plotter = classes.Plotter(xplot,xlabel,yplot,ylabel,axistitle,weathernumber)
        plotter.plot_validation(xvalid,yvalid)
    else: 
        print('No plot functionality has been defined.')
