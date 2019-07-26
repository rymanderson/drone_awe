import classes as c
import numpy as np

# Test WeatherType class
print("WeatherType Test")
params      = {'hello':4,'there':5}
weathertype = c.WeatherType(params,{'hello':None,'there':None})
print("    Success: TRUE")

# Test Rain class
print("Rain Test")
params     = {'LWC':2.0,'dropsize':0.01,'WVC':23.0}
rain       = c.Rain(params)
print("    Success: TRUE")

# Test Temperature class
print("Temperature Test")
params     = {'temperature':300,'temperaturesealevel':308}
temperature= c.Temperature(params)
print("    Success: TRUE")

# Test Humidity class
print("Humidity Test")
params     = {'humidityrelative':50.0,'humidityabsolute':0.9}
humidity   = c.Humidity(params)
print("    Success: TRUE")

# Test Wind class
print("Wind Test")
params     = {'heading':270.0,'speednortheast':9.0}
wind       = c.Wind(params)
print("No vector or downdraft inputs: ", wind.params)
params     = {'velocityvector':[3.0,4.0,-1.0]}
wind       = c.Wind(params)
print("No heading, downdraftspeed, or speednortheast inputs: ", wind.params)
print("    Success: TRUE")

# Test Gust class
print("Gust Test")
params     = {'amplitude':2.0,'frequency':0.1}
gust       = c.Gust(params)
print("    Success: TRUE")

# Test Ice class
print("Ice Test")
params     = {'':None}
ice        = c.Ice(params)
print("    Success: TRUE")
