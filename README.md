# Drone-Models
Development of drone models accounting for weather and battery characteristics.

## Instructions
This section contains instructions for how to run awesomeModelName.py.

## Directories Overview
### Drones/
This directory contains the following:
* the `Drone` class
	* class variables contain:
		* data sheet specifications of specific drone models (e.g., the Mavic Pro 4), including 
			* battery size
			* battery type
			* range under specified conditions
			* fixed wing or rotary wing
			* etc.
	* methods calculate certain characteristics based on available information
* the `params/` directory
	* contains `.param` text files containing specifications of each drone to be modeled

### Batteries/
This directory contains the following:
* the `Battery` class
	* class variables describe:
		* properties of specific batteries used to model their discharge characteristics, including
			* battery type
			* number of cells
			* low, nominal, and charged cell voltages
		* Real-time discharge characteristics for simulation, including
			* instantaneous voltage
			* instantaneous current
			* instantaneous state of charge
			* current state of health
	* methods are used to update class variables using information from the `params/` directory
* the `params/` directory
	* contains `.param` text files containing information about different battery types

### Weather/
This directory contains the following:
* scripts used to model weather effects on drones
* classes corresponding to each weather effect (e.g., rain, icing (though this might be a function of temperature and humidity- maybe it has a class that checks to see if icing is likely based on other weather characteristics), wind, etc.) that contain the following:

    * class variables characterizing the weather (e.g., )
    * methods used to modify the `Power` class
		* these probably take an instance of the `Drone` class as an argument
* the `Weather` class
	* class variables describe
		* an instance of each weather effect to be modeled
	* methods may be used to check for icing, or other miscellaneous needs
* the `params/` directory
	* contains `.param` text files containing information about different weather phenomena

### Simulation/
This directory contains the following:
* scripts used to model the range of a particular drone
* the `Power` class
	* class variables describe
		* baseline power consumption
		* an array of 'correction' objects used to modify the power consumption class variable due to weather or other effects (these could be the weather effect classes, actually)
		* total power consumption
	* class methods are used to update the total power consumption class variable
	* class variables may be constant (for speed and simplicity) or vary with time in simulation
* the `Simulation` class
	* class variables describe
		* start time
		* end time
		* timestep
		* current timestep index
		* current time
		* an instance of the `Drone` class
		* an instance of the `Battery` class
		* an instance of the `Power` class
		* an instance of the `Weather` class
	* methods are used to run and store simulation information
NOTE: the model is based on power consumption to accomodate future development. The `Power` class is designed to receive an indefinite number of modifications based on weather effects

* the `Plotter` class
	* class variables describe
		* x-axis variable
		* x-axis label
		* y-axis variable
		* y-axis label
		* other miscellaneous plotting parameters
	* methods are used to make plots
* an `awesomeModelName_importer.py` module importer script
	* this script imports the classes and definitions defined above
	* after this script is run, classes and functions can be accessed and used to start simulations
* an `awesomeModelName_exe.py` script 
	* this script is a sample simulation script containing
		* the module importer script
		* a sample simulation setup, where the following are specified:
			* the drone
			* weather
			* other variables
NOTE: `awesomeModelName_importer.py` and `awesomeModelName_exe.py` should be simple since data describing drone specifications, weather effects, battery parameters, etc. should already be contained `params/` directories

## Classes
This section contains a detailed description of each class.