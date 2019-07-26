import classes
import functions

dronename = 'dji-Mavic2'
drone_params = functions.getParams('Drone','paramlist.param',dronename + '.param',' ')
conversions  = functions.getParams('Drone','paramlist.param','conversions.param',' ')
drone = classes.Drone(dronename,drone_params,conversions)
drone.getEfficiencyPropulsive()