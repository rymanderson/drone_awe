# Drone-Models

Development of drone models accounting for weather and battery characteristics.

## Instructions

To run this model, it is necessary to first edit the file `settings.txt` found in `params/Simulation/`. Specifics on these parameters can be found in the simulation subsection of the Parameter Files section in this document. Eventually `settings.txt` will become a GUI interface. 

Once that has been edited, simply run exe.py from the main `Drone-Models` directory and it should run everything automatically and result in a nice plot of the model (and of validation comparison data if validation is set to true).  

## Classes

This section contains a detailed description of each class, all contained in Classes.py.

* the `Drone` class

	* class variables contain:

		* data sheet specifications of specific drone models (e.g., the Mavic 2 Pro), including

			* battery size
			* battery type
			* range under specified conditions
			* fixed wing or rotary wing
			* etc.

	* methods calculate certain characteristics based on available information

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

* the `Power` class

	* class variables describe

		* baseline power consumption
		* an array of 'correction' objects used to modify the power consumption class variable due to weather or other effects (these could be the weather effect classes, actually)
		* total power consumption

	* methods are used to

		* update the total power consumption class variable
		* append `PowerCorrection` objects to `Power` objects using the `addCorrection` method
		* throw an error if `addCorrection` attempts to append a time-variant `PowerCorrection` object to a time-invariant `Power` object
		* the `PowerCorrection` method
			* adjustments to the baseline power requirements of the drone
			* whether the simulation is time-variant or time-invariant
			* methods perform miscellaneous book-keeping functions

* the `Weather` class

	* class variables describe

		* an instance of each weather effect to be modeled

	* methods may be used to check for icing, or other miscellaneous needs

* the `WeatherType` class

	* class consists of the following sub-classes:

		* the `Rain` class

			* class variables describe 

				* droplet size
				* Liquid Water Content (LWC)
				* other parameters needed to model rain

		* the `Temperature` class

			* class variables describe 

				* new temperature?

		* the `Humidity` class

			* class variables describe 

				* relative humidity

		* the `Wind` class

			* class variables describe 

				* wind speed
				* wind direction
				* amount of turbulence?
				* variation in speed and/or direction?
		
		* the `Gust` class

			* class variables describe

				* frequency
				* amplitude
		
		* the `Ice` class
			
			* Because of modeling difficulty, this class does not attempt yet to model icing effects. It may in the future be used to identify if icing conditions are present. 

* the `Mission` class

	* class variables describe
		
		* mission speed

* the `Simulation` class

  * class variables describe
	* start time
	* end time
	* timestep
	* current timestep index
	* current time
	* methods are used to run and store simulation information, including procedures to get range and endurance

NOTE: the model is based on power consumption to accomodate future development. The `Power` class is designed to receive an indefinite number of modifications based on weather effects

* the `Plotter` class

	* plots results according to labels and titles specified by the user in the `settings.txt` file. Methods can plot a line or scatter plots. The validation method plots results on top of specified validation data.

## Units

Properties and their respective units are converted within the simulation to SI units, and then converted back. Those units are:

<!-- These could probably be better organized -->

### Electricity

* Capacity: _milliamp-hours [mAh]_
* Voltage: _volts [V]_
* Current: _amperes [A]_
* Resistance: _ohms [&Omega;]_

### Mechanics

* Velocity/Speed: _meters per second [m/s]_
* Power: _watts [W]_
* Endurance or Flight time: _minutes [min]_
* Altitude: _meters [m]_
* mass: _kilograms [kg]_
	* this is often referred to as take off weight, despite actually referring to mass <!-- I feel weird quoting a weight in mass units, but I put this here because the .param files seem to use kg. Which may be fine. Let me know if you have any thoughts. :)  Yeah, that is kind of strange. -->

### Miscellaneous

* Temperature: _degrees Celcius [&deg;C]_
* Wind Resistance: _meters per second [m/s]_ 
	*refers to the maximum wind speed the drone can resist
* Battery re-charge time: _minutes [min]_

## Parameter Files

### Simulation

* The settings list file contains all the necessary simulation parameters the code will look through before running simulations. These include: 

	* validation (True/False)

		* If validation is True, the program reads in the next value, validationcase, and looks for the settings file under that directory in params/Validation/, and ignores the rest of the current settings file. 

	* validationcase
		
		* Specifies which validation case is being tested. Irrelevant if validation is False.
	
	* drone (True/False)
	* dronename
	* battery state of health
	* battery beginning state of charge
	* altitude
	* sea level temperature
	* rain (True/False)

		* droplet size
		* LWC

	* Temperature

		* New temperature

	* Humidity

		* Relative humidity

	* Wind

		* Wind speed
		* Wind direction

	* Icing
	* Timestep
	* plot (True/False)

		* xlabel
		* ylabel
		* axis title
	* simulation type 
		* simple is the only option for now
	* range of vx-values
		* xbegin
		* xend
		* xnumber
			* the simulation will loop through xnumber of the model from xbeginning to xend, according to what is put as the xlabel (also the x-variable)

NOTE: For the plotting x- and y-labels, choose from the following parameters to plot:

* range
* endurance
* payload

### Drone
There are `.param` files in the Drone directory each contain parameters specific to each drone. They are labled with the company's name or abbreviation in lowercase letters followed by a dash (-) and then the name of the drone (e.g., dji-Mavic2). New drones needing to be tested can follow similar formats, with a space (`" "`) delimiter.

### Batteries
Similar to the Drone `.param` files, battery `.param` files exist for each type of battery tested. These parameters assist in determining the discharge rate over time and amount of specific energy and power available for different types of batteries.

## functions.py - Commonly used functions
* `getparams` reads in a .txt or .csv file and outputs a dictionary with keys from a specified list and values from the specified parameter file.
* `getXandY()` reads in data from a validation case and saves the contents to lists for x and y. This function assumes the first row contains labels and ignores them.

## Testing
* test_power.py
* test_drone.py
* test_plotter.py

## Future Work

This section is to be used to record ideas for future development that cannot be immediately implemented due to time constraints.

* calculate propulsive efficiency at max range and max endurance and interpolate between the two
* go weather by weather and determine the appropriate model to be used