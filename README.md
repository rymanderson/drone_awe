# Drone-AWE: Drone Applications in Weather Environments

This repository is being developed to advance the state of the art of drone performance predictions under a wide range of weather and battery conditions with mission planning in mind. An interactive GUI is available at http://droneecon.com. This README will explain general usage as well as the underlying theory for the models used, as currently implemented. Note that *Drone-AWE* is a work in progress (see **Future Work**). This document will explain general usage for the latter, as well as the engineering theory at work behind the scenes.

## Installation

Drone-AWE may be used in two ways: using the GUI available at http://droneecon.com, or by using the source code directly. To use the source code directly, install using

```python3
pip install drone_awe
```

## Usage

### Basic Usage

Base functionality is achieved by running the following:

```python3
import drone_awe
m = drone_awe.model(args)
m.simulate()
```

where `args` is a dictionary containing simulation parameters. Note that for every parameter not specified in `args`, a default value is used. It is possible to run the default simulation by passing an empty dictionary into the `drone_awe.model()` method:

```python3
import drone_awe
m = drone_awe.model({})
m.simulate()
```

Settable dictionary keys with example values for `args` include the following:

```python3
{
    "validation":False,
    "validationcase":"DiFranco2016",
    "dronename":"dji-Mavic2",
    "batterytechnology":"near-future",
    "stateofhealth":90.0,
    "startstateofcharge":100.0,
    "rain":False,
    "dropsize":1.0,
    "liquidwatercontent":1.0,
    "temperature":15.0,
    "wind":False,
    "windspeed":10.0,
    "winddirection":0.0,
    "relativehumidity":0.0,
    "mission": {
            "missionspeed": 10.0,
            "altitude":100.0,
            "heading":0.0,
            "payload:0.0
        },
    "timestep":1,
    "plot":True,
    "xlabel":"missionspeed",
    "ylabel":"power",
    "title":"First_test",
    "simulationtype":"simple",
    "model":"abdilla",
    "xvals":[0,1,2,3,4,5],
    "weathereffect":"temperature",
    "weathervals":[10,20,30,40]
}
```

These are explained further in **Settings**. After a `drone_awe.model` object is instantiated, settings can be changed by modifying the `input` class variable as:

```python3
m.input['xlabel'] = 'payload'
m.input['validationcase'] = 'Stolaroff2018'
```

The simulation must then be re-run using:

```python3
m.simulate()
```

### Plotting

An optional boolean may be set to produce a plot on runtime by running:

```python3
m = drone_awe.model({},plot=True)
```

when a `drone_awe.model` is instantiated. Alternatively, the `drone_awe.model`'s `plot` variable may be set at any time:

```python3
m.plot = True
```

### Settings

The following subsections explain the use of various settings of a `drone_awe.model` object.

#### `'dronename'`

`'dronename'` must be set to a string that is defined in the drone database. To access a dictionary containing all supported drones, use the following methods:

```python3
droneDictionary = drone_awe.drones
```

If a drone or validation case does not exist in the drone database, it may still be used by setting a custom dictionary:

```python3
m.drone = droneDictionary
```

where `droneDictionary` is a dictionary with the same keys as the following:

```python3
{
    'id': 'dji-Mavic2',
    'wingtype': 'rotary',
    'diagonal': 0.354,
    'takeoffweight': 0.907,
    'speedmax': 20,
    'altitudemax': 6000,
    'endurancemax': 31,
    'endurancemaxspeed': 6.944,
    'endurancemaxhover': 29,
    'rangemax': 18000,
    'rangemaxspeed': 13.889,
    'temperaturemin': -10,
    'chargerpowerrating': 60,
    'batterytype': 'LiPo',
    'batterycapacity': 3850,
    'batteryvoltage': 15.4,
    'batterycells': 4,
    'batteryenergy': 59.29,
    'batterymass': 0.297,
    'waterproof': 'no',
    'windspeedmax': 10.8,
    'batteryrechargetime': 90,
    'rotorquantity': 4,
    'rotordiameter': 0.2,
    'cruisespeed': 6.94,
    'payload': 0.0,
    'length': 0.322,
    'width': 0.242,
    'height': 0.084
}
```

Note that not all parameters need be specified, but if a simulation is run that requires unspecified parameters, the model will not run.

#### `'validation'` and `'validationcase'`

To view validation cases for `drone_awe` models, set `'validation'=True` and `'validationcase'` to a string contained in the validation case database. To access a dictionary containing all supported validation cases, run the following:

```python3
validationDictionary = drone_awe.validationdata()
```

Additionally, `validationCaseDictionary` is a dictionary with the following format:

```python3
{
    'id': 'Stolaroff2018',
    'xvalid': [1.358974359,1.7179487179,1.6923076923,1.7692307692,1.7692307692,1.8205128205,1.8717948718,1.8974358974,1.9230769231,1.9743589744,2.0512820513,2.1025641026,2.1025641026,2.1794871795,2.2820512821,2.3076923077,2.3333333333,2.4358974359,2.5384615385,2.5384615385,2.6666666667,2.7948717949,2.9487179487,3.1025641026,3.2564102564,3.358974359,3.4615384615,3.641025641,3.7435897436,3.8717948718,4.0512820513,4.2564102564,4.358974359,4.641025641,4.8974358974,5.8461538462,6,6.2820512821,6.5384615385,7.1538461538,7.4871794872,7.8205128205,8.4871794872,9.1794871795,9.7692307692,10.2564102564,10.6923076923,11.1282051282,11.5897435897],
    'yvalid': [133.7510729614,149.2375075128,142.6477385276,144.5897859156,140.9450186657,141.7583868756,142.3874173587,140.4340604922,139.3656195241,139.5335348046,140.2321151941,139.9589946754,136.1689649626,138.3618186589,141.3711557508,139.4001574523,136.7492021569,136.5131929807,137.712167001,133.136399421,135.7755034665,136.58240428,141.8694500174,143.7278953027,145.8751386173,144.4503136349,143.5492800366,144.4966689523,140.8944307591,139.7885398414,140.0842115956,140.6905892611,140.7593265104,140.4412051028,148.3297356325,161.0751453894,160.7738527567,158.8079335653,157.0439596719,161.2248774665,165.1381601781,164.7032531681,177.1408013138,175.6982333173,190.095233258,197.4485952036,209.0736216573,203.7684942987,215.7551870381],
    'drone': {
        'validationcase':'Stolaroff2018',
        'wingtype': 'rotary',
        'rotorquantity': 4,
        'takeoffweight': 1.3,
        'batterytype': 'LiPo',
        'batteryvoltage': 11.1,
        'batterymass': 0.262,
        'props': '10x4.7',
        'endurancemaxhover': 16,
        'payloadmax': 0.4,
        'batterycapacity': 5500,
        'payload': 0.0,
        'rotordiameter': 0.254,
        'batterycells': 3,
        'length': 0.280,
        'width': 0.140,
        'height': 0.100,
        'waterproof': 'no'
        },
    'settings': {
        'dronename': 'drone',
        'stateofhealth': 100.0,
        'startstateofcharge': 100.0,
        'altitude': 100.0,
        'temperaturesealevel': 15.0,
        'rain': False,
        'dropsize': 0.0,
        'liquidwatercontent': 1.0,
        'temperature': 15.0,
        'wind': False,
        'windspeed': 0.0,
        'winddirection': 0.0,
        'relativehumidity': 85.0,
        'icing': False,
        "mission": {
            "missionspeed": 10.0
        },
        'timestep': 1,
        'xlabel': 'missionspeed',
        'ylabel': 'alpha',
        'title': 'Stolaroff',
        'simulationtype': 'simple',
        'model': 'abdilla',
        'xvals': [0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0],
        'validation': False,
        'validationcase': 'Stolaroff2018',
        'batterytechnology': 'current'
        }
}
```

Note that the `'altitude'` setting of validation cases is not always known and is set by default to `'0'`.

#### `'L/D'` and `'propulsiveefficiency'`

When simulating a fixed-wing drone, the current model requires input parameters for the lift-to-drag ratio (`'L/D'`) and the efficiency of the propulsive system (`'propulsiveefficiency'`). Ballpark estimates for small UAVs are on the order of 10 and 30%, respectively, although these values heavily depend on the specific geometries and propulsive systems of the drones.

#### Battery parameters

#### `'batterytechnology'`
`'batterytechnology'` must be set to one of the following options: `'current'`, `'near-future'`, or `'far-future'`. `'current'` will use the current battery capacity as specified in `drone.params`. `'near-future'` projects a future capacity in five years based on an increase of 3.5% of battery capacity each year. `'far-future'` projects battery capacity for lithium-air batteries, about ten times the capacity of current LiPo batteries.

#### `'stateofhealth'` and `'stateofcharge'`
`'stateofhealth'` refers to the amount of capacity contained within the battery relative to its initial capacity at first use. 100% would be a new battery. LiPo batteries (most commonly used battery for drones) are typically retired when their state of health falls to 80-85%. At lower states of health, the battery capacity is decreased and cycles through a full charge more quickly.

`'stateofcharge'` refers to the current capacity level of the battery, where 100% is fully charged (regardless of how much the battery has been used before), and 0% is a dead battery.

#### Weather parameters

#### `'dropsize'`, `'liquidwatercontent'` and `'rainfallrate'`
These parameters specifiy rain characteristics. If there is no rain in the simulation, set these to 0.0. 

`'dropsize'` refers to the diameter of a raindrop (assuming a spherical shape), with units of meters. `'liquidwatercontent'` refers to the amount of water in a given volume of air, with units of kg/m^3. `'rainfallrate'` is the rate of rainfall in units of mm/hr, which translates to litres per cubic meter per hour.

Only one of `'liquidatercontent'` and `'rainfallrate'` needs to be specified. 

#### `'temperature'`
Units of temperature are given in degrees Celcius.

#### `'relativehumidity'`

`'relativehumidity'` refers to water content in the air, ranging from 0% to 100%.

#### Mission parameters
`'mission'` is a dictionary with the following keys: 

* `'missionspeed'`&mdash; the cruise velocity of the drone, with units of m/s.

* `'altitude'`&mdash; in units of meters.

* `'heading'`&mdash; the heading angle of the drone, in units of ___________.

* `'payload'`&mdash; the mass of any exra payload attached to the drone (e.g., camera), in units of kg.

#### Plotting parameters
`'xlabel'` and `'ylabel'` not only specify axis labels of the plot, but also tell the simulation what to solve for (`'ylabel'`) and what parameter to loop through (`'xlabel'`). `'ylabel'` can be `range`, `endurance`, or `power`. `'xlabel'` can be any variable from the following list:

```
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
"missionspeed"
```

`'xvals'` is a list containing all the values of x the simulation will loop through and produce a single data point for. 

`'zlabel'` can be specified to any of the variables used for the `'xlabel'` parameter (see list above) to loop through. This will produce a curve for each of the values in `'zvals'` and plot them on the same figure. This is useful for looping through a weather effect, such as temperature or relative humidity, to see their effects on a plot of range vs payload.

`'title'` is the title displayed on the top of the plot, followed by the current date. Currently, the title needs to be one word. A title of multiple words can have each word separated by an underscore en lieu of a space.

## Theory

`Drone AWE` is intended to predict performance parameters of rotary and fixed-wing drones based on readily-available specifications. To accomplish this, the power requirements for a given flight maneuver is calculated and used to predict battery drain behavior. Then, parameters such as range and endurance can be calculated. Comprehensive state data is collected at each step of a simulation, providing a versatile data set as output for study.

### Power

#### Rotary Drones

For rotary drones, power consumption is predicted in two steps. First, model parameters are calibrated based on known specifications, like maximum hover time and rotor diameter. Then, momentum theory is used to predict the power consumption.

#### Calibration

#### Momentum Theory

According to rotor momentum theory, five equations govern the flight of a rotorcraft:

1. $T = \sqrt{W^2 + D^2}$
2. $tan(\alpha) = \frac{D}{W}$
3. $D = \frac{1}{2} \rho V_\infty^2 C_D A_{\bot}$
4. $v_i = \frac{T}{2 A_{rotor} N_{rotor} \rho \sqrt{V_\infty^2 cos^2(\alpha) + (V_\infty sin(\alpha) + v_i^2)}}$
5. $A_{\bot} = A_{front} cos(\alpha) + A_{top} sin(\alpha)$

A modified set of these equations were applied to rotary drones by Stolaroff, Samaras, O'Neill et. al. in [1]. `Drone AWE` predicts the power consumption of a rotary drone by solving this system. GEKKO, a software package introduced in [2], is used to solve for:

* $T$ - thrust per rotor
* $\alpha$ - angle of attack
* $D$ - drag
* $v_i$ - rotor induced velocity
* $A_{\bot}$ - planform area of the drone perpendicular to its velocity

given the following:

* $W$ - effective weight
* $V_\infty$ - drone speed
* $C_D$ - drag coefficient
* $A_{rotor}$ - rotor area
* $N_{rotor}$ - number of rotors
* $\rho$ - air density
* $A_{front}$ - frontal drone area
* $A_{top}$ - drone top area

Note that the model currently assumes that the drag coefficient does not change with velocity and/or angle of attack.

The power requirement for the hover case is calculated according to another equation from [1]:

* $P_{hover} = \frac{W^{3/2}}{\eta \sqrt{2N_{rotor} A_{rotor} \rho}}$

where $\eta$ is the propulsive efficiency calculated as:

* $\eta = \frac{P_{hover}}{P_{hover, actual}}

where $P_{hover, actual}$ is predicted using:

* $P_{hover, actual} = $

This $\eta$ is used to predict power consumption for the non-hover case as:

* $P = \eta T(V_\infty sin(\alpha) + v_i)$

#### Fixed-Wing Drones

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

	* note that "takeoff weight" is measured in mass units <!-- I feel weird quoting a weight in mass units, but I put this here because the .param files seem to use kg. Which may be fine. Let me know if you have any thoughts. :)  Yeah, that is kind of strange. -->

### Miscellaneous

* Temperature: _degrees Celcius [&deg;C]_
* Wind Resistance: _meters per second [m/s]_

	*refers to the maximum wind speed rating for the drone

* Battery re-charge time: _minutes [min]_

## functions.py &mdash; Commonly used functions

* `getparams` reads in a .txt or .csv file and outputs a dictionary with keys from a specified list and values from the specified parameter file.
* `getXandY()` reads in data from a validation case and saves the contents to lists for x and y. This function assumes the first row contains labels and ignores them.
* `interpolate()` does a simple linear interpolation with inputs of 2 x-values, 2 y-values, and the x-value of the interpolated value.

## Testing
* test_power.py
* test_drone.py
* test_plotter.py

## Future Work

This section is to be used to record ideas for future development that cannot be immediately implemented due to time constraints.

* calculate propulsive efficiency at max range and max endurance and interpolate between the two
* go weather by weather and determine the appropriate model to be used

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
        * droplet size (rain)
        * Liquid Water Content (LWC) (rain)
        * rainfall rate (rain)

    * the `__updateRain` method

		* For all drone types, this calculates the momentum imparted to the drone from falling rain droplets based on their size and liquid water content. 

    * the `__getWebernumber` method
        * Calculates weber number based on rain density, velocity, diameter, and frequency. This is used to obtain momentum in the __updateRain method.

    * the `__getSurfaceTension` method
        * empirically interpolates surfaces tension based on current temperature

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

## References

1. Beal, L.D.R., Hill, D., Martin, R.A., and Hedengren, J. D., GEKKO Optimization Suite, Processes, Volume 6, Number 8, 2018, doi: 10.3390/pr6080106.
2. Stolaroff, J. K., Samaras, C., O’Neill, E. R., Lubers, A., Mitchell, A. S., & Ceperley, D. (2018). Energy use and life cycle greenhouse gas emissions of drones for commercial package delivery. Nature Communications, 9(1), 1–13. https://doi.org/10.1038/s41467-017-02411-5
3. 