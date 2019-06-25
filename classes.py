import datetime
import matplotlib.pyplot as plt
import numpy as np

# insert classes here


class Drone:
    'Class used to store key drone characteristics'

    name = None
    params = None
    conversions = None
    correctunits = False

    def __init__(self, name, params, conversions):
        self.name = name
        self.params = params
        # this input is obtained using functions.getParams() in exe.py
        self.conversions = conversions
        self.__convertUnits()

    def __convertUnits(self):
        if not self.correctunits:
            for spec in self.conversions:
                if spec in self.params:
                    self.params[spec] = self.params[spec] * \
                        self.conversions[spec]
                else:
                    print("Tried to convert:", spec, "- skipping for now.")
            self.correctunits = True
        else:
            print("~~~~~ WARNING: problem converting units: skipping ~~~~~")


print("Successfully imported `Drone` class")


class Battery:
    'Class used to track battery characteristics and performance during simulation'

    params = {
                'capacity':None,        # Ampere-hours
                'soc':None,             # state of charge (in percent nominal capacity)
                'startsoc':None,        # state of charge at simulation start
                'soh':None,             # state of health (actual capacity divided by ideal capacity)
                'batterytype':None,     # possible values include LiPo, Li-ion, NiCd, NiMH, SLA
                'voltage':None,         # Volts; this is the instantaneous voltage
                'voltagemean':None,     # Volts; average voltage used for time-invariant simulations
                'voltagecharged':None,  # Volts
                'voltagedead':None,     # Volts
                'current':None          # Amperes; this is the instantaneous current
            }

    # constructor
    # default value for soh is based on the assumption that batteries are retired at a soh of 80%
    def __init__(self, drone, soh=90.0, startsoc=100.0):
        # import parameters from drone object
        self.batterytype = drone.params['batterytype']
        self.voltagemean = drone.params['batteryvoltage']
        self.capacity = drone.params['batterycapacity']

        # update parameters
        self.update()

        # estimate list lengths for prior memory allocation

    def update(self):
        print("still working on Battery.update method")

    # time-variant methods go here:
    def updateLoad(self, power):
        self.current = power.params['power'] / self.voltage

    def discharge(self, power, timestep):
        self.soc = self.soc - power.params['power'] * timestep
        self.current = power.params['power'] / self.voltage


print("Successfully imported `Battery` class")


class Power:
    'Class used to predict the drone\'s power requirement'

    # class variables go here:
    params = {
        'efficiencypropulsive': None,
        'power':None
    }

    # methods go here:
    def __init__(self, drone, weather):
        pass

    def update(self, drone, weather, model, mission):
        self.__updateEfficiencyPropulsive(drone, mission)
        if model == 'dandrea':
            self.__getPowerDandrea(drone, weather)
        elif model == 'abdilla':
            self.__getPowerAbdilla(drone, weather)
        else:
            # raise Exception(f"~~~~~ ERROR: model { model } not available ~~~~~") #
            raise Exception("~~~~~ ERROR: model '" +
                            model + "' not available ~~~~~")

    # super simple estimate for power from D'Andrea `Can Drones Deliver` *****Doesn't work well*******
    def __getPowerDandrea(self, drone, weather):
        powerelectronics = 0.1          # kW, estimate from paper
        L_D = 3.0          # quick estimate for initial functionality TODO: Change this to something more scientific
        self.params['power'] = (drone.params['takeoffweight'] + drone.params['payload']) * drone.params['endurancemaxspeed'] / (
            370.0 * self.params['efficiencypropulsive'] * L_D) - powerelectronics

    # slightly more complicated estimate for power
    def __getPowerAbdilla(self, drone, weather):
        # print("power:           takeoffweight is        ",drone.params['takeoffweight'])
        # print("weather:         gravitationconstant is  ",weather.params['gravitationconstant'])
        # print("power:           etapropulsive is        ",self.params['efficiencypropulsive'])
        # print("drone:           payload is              ",drone.params['payload'])
        # print("drone:           rotorquantity is        ",drone.params['rotorquantity'])
        # print("weather:         airdensity is           ",weather.params['airdensity'])
        # print("")
        self.power = (drone.params['takeoffweight']/weather.params['gravitationconstant'] + drone.params['payload'])**1.5 / \
                     (self.params['efficiencypropulsive'] * drone.params['rotordiameter'] / 2.0) * \
            weather.params['gravitationconstant']**1.5 / \
            np.sqrt(2 * drone.params['rotorquantity'] *
                    weather.params['airdensity'] * np.pi) + \
            0.0  # Camera power consumption estimate

    def __updateEfficiencyPropulsive(self, drone, mission):
        # default value:
        if 'endurancemax' not in drone.params or 'endurancemaxspeed' not in drone.params:
            self.params['efficiencypropulsive'] = 0.4
        else:
            # get efficiency at max endurance conditions
            etamaxendurance = self.__getEfficiencyPropulsive(
                drone, drone.params['endurancemax'])
            vmaxendurance = drone.params['endurancemaxspeed']
            # get efficiency at max range conditions
            etamaxrange = self.__getEfficiencyPropulsive(
                drone, drone.params['rangemax']/drone.params['rangemaxspeed'])
            vmaxrange = drone.params['rangemaxspeed']
            # interpolate for the current velocity
            self.params['efficiencypropulsive'] = etamaxendurance - (vmaxendurance - mission.params['missionspeed']) / \
                (vmaxendurance - vmaxrange) * (etamaxendurance - etamaxrange)
            # print("Drone:           etamaxendurance is ",etamaxendurance)
            # print("Drone:           etamaxrange is     ",etamaxrange)
            # print("Mission:         mission speed is   ",mission.params['speed'])
            # print("Drone:           etamission is      ",self.params['efficiencypropulsive'])
            # print("Drone:           vmaxendurance is   ",vmaxendurance)
            # print("Drone:           vmaxrange is       ",vmaxrange)
            # print("")
            # print("")

    def __getEfficiencyPropulsive(self, drone, endurance):
        # Analyzing a single propeller
        thrust          = drone.params['takeoffweight']/drone.params['rotorquantity']
        batteryenergy   = drone.params['batteryenergy'] / \
            drone.params['rotorquantity']
        rotorarea       = drone.params['rotordiameter']**2/4*np.pi
        airdensity      = 1.225  # assuming air density was equal to 1.225 kg/m3 during drone testing

        poweractual     = batteryenergy/endurance
        powerideal      = thrust * np.sqrt(thrust/(2*rotorarea*airdensity))
        efficiency      = powerideal/poweractual
        # prediction based on momentum theory for hover case
        # slides from https://fenix.tecnico.ulisboa.pt/downloadFile/282093452028191/3-Momentum%20Theory%20in%20hover.pdf

        return efficiency


print("Successfully imported `Power` class")


class Weather:
    '''
    Class describing ambient weather conditions and is used to predict the drone\'s power requirement
    '''

    # miscellaneous weather
    # a list of weather objects (e.g., rain, icing, humidity, etc.)
    weatherlist = []

    # ambient air quality
    params = {
        'airdensity': 1.225,
        'gravitationconstant': 9.807
    }

    # airdensity      = 1.225 #kg/m^3
    # gravity         = 9.806655 #m/s^2
    # altitude = 100 #m
    # temperature_sl = 15
    # temperature = 15
    # pressure_sl = 1
    # pressure = 101300 #Pa
    # rain = 'False'
    # ice = 'False'
    # wind = 0
    # humidity = 0

    # need to figure out how to access current or at least average weather
    # data based on lattitude, longitude, altitude?
    # then calculate pressure, temp, and density based on altitude (see
    # eqns 1.21 - 1.22 in Dr. Ning's book
    # or do above with average values of temperature/pressure/density for a given location
    # or none of the above - we vary temperature and assume either constant pressure or density (likely pressure)

    # How about this: we assume that we know the pressure/density/temperature at ground level.

    # methods go here:
    def __init__(self, altitude, temperaturesealevel):  # keeping it simple to begin with
        pass
        # self.altitude = altitude
        # self.temperaturesealevel = temperaturesealevel
        # self.temperature = self.temperaturesealevel - 71.5 + 2*np.log(1 + np.exp(35.75 - 3.25*self.altitude) + np.exp(-3 + 0.0003 * self.altitude**3))
        # self.pressure = self.pressure_sl * np.exp(-0.118 * self.altitude - (0.0015*self.altitude**2) / (1 - 0.018*self.altitude + 0.0011 * self.altitude**2))

    def update(self):
        # update independent parameters
        for weathertype in self.weatherlist:
            for param in weathertype.params:
                self.params[param] = weathertype.params[param]
        # update dependent parameters
        for weathertype in self.weatherlist: 
            pass

    def getStandardAtmosphere(self, altitude):
        '''This function currently assumes STP conditions at sea level and should probably be adjusted to use ground level conditions as a baseline'''

        # set sea level parameters
        temperaturesealevel             = 288.15               # Kelvin
        pressuresealevel                = 1.01325e5            # Pascals
        gravitationconstantsealevel     = 9.80665              # m/s2
        specificheatratio               = 1.4
        airgasconstant                  = 287.053              # J/(kg-K)
        S                               = 110.4                # K
        beta                            = 1.458e-6             # kg/(smK^1/2)

        # compute parameters at altitude
        altitude            = altitude / 1000.0    # convert to kilometers
        pressure            = pressuresealevel *\
                                np.exp(-0.118 * altitude - \
                                    (0.0015 * altitude**2) / \
                                    (1 - 0.018 * altitude + \
                                        0.0011 * altitude**2) \
                                    )
        temperature         = temperaturesealevel - 71.5 + 2 * \
            np.log(1 + np.exp(35.75 - 3.25 * altitude) +
                   np.exp(-3 + 0.0003 * altitude**3))
        airdensity          = pressure / (airgasconstant * temperature)
        speedsound          = np.sqrt(specificheatratio * airgasconstant * temperature)
        dynamicviscocity    = beta * temperature**1.5 / (temperature + S)

        return (pressure, temperature, airdensity, speedsound, dynamicviscocity)


print("Successfully imported `Weather` class")


class WeatherType:
    params = {}

    def __init__(self, params, paramnames):
        'constructor'
        for name in paramnames:
            if not(name in params):
                print("~~~~~ WARNING: parameter ", name,
                      " missing in WeatherType declaration ~~~~~")
                print(
                    "~~~~~        creating as None type ~~~~~")
                params[name] = None
        if len(paramnames) < len(params):
            print("~~~~~ WARNING: "+str(-len(paramnames)+len(params)) +
                  " extra parameters in WeatherType declaration ~~~~~")
        self.params = params


class Rain(WeatherType):
    '''
    Class used to define rain characteristics. Params variables include:
    * variable (parameter) [units] {example}
    * LWC (liquid water content) [kg/m3]
    * dropsize [m]
    * WVC (water vapor content) []

    Rain is expected to affect:
    * lift
    * drag
    * weight
    * impart downward momentum
    '''

    def __init__(self, params):
        paramnames = ['LWC', 'dropsize', 'WVC']
        WeatherType.__init__(self, params, paramnames)


print("Successfully imported `Rain` class")


class Temperature(WeatherType):
    '''
    Class used to define temperature characteristics. Params variables include:
    * variable (parameter) [units] {example}
    * temperature [K]
    * temperaturesealevel [K]

    Temperature is expected to affect:
    * air density
    * battery performance
    '''

    def __init__(self, params):
        paramnames = ['temperature', 'temperaturesealevel']
        WeatherType.__init__(self, params, paramnames)


print("Successfully imported `Temperature` class")


class Humidity(WeatherType):
    '''
    Class used to define humidity characteristics. Params variables include:
    * variable (parameter) [units] {example}
    * humidityrelative (relative humidity) [%]
    * humidityabsolute (absolute humidity) []

    Humidity is expected to affect:
    * air density
    * condensation
    '''

    def __init__(self, params):
        paramnames = ['humidityrelative', 'humidityabsolute']
        WeatherType.__init__(self, params, paramnames)


print("Successfully imported `Humidity` class")


class Wind(WeatherType):
    '''
    Class used to define wind characteristics. Params variables include:
    * variable (parameter) [units] {example}
    * velocityvector (wind direction and magnitude) [<north,east,down> m/s]
    * heading [degrees East of North]
    * downdraftspeed (wind velocity component downwards) [m/s]
    * speednortheast (wind speed neglecting downdraft component) [m/s]

    Wind is expected to affect:
    * No downdraft: relative airspeed and therefore effective velocity
    * With downdraft: ?? effective climb/descent ??
    '''

    def __init__(self, params):
        paramnames = ['velocityvector', 'heading',
                      'downdraftspeed', 'speednortheast']
        WeatherType.__init__(self, params, paramnames)
        if self.params['velocityvector'] == None:
            if not(self.params['heading'] == None) and not(self.params['speednortheast'] == None):
                if self.params['downdraftspeed'] == None:
                    # set downdraftspeed to 0 by default
                    self.params['downdraftspeed'] = 0.0
                velocityvector = [self.params['speednortheast'],
                                  0.0, self.params['downdraftspeed']]
                theta = self.params['heading'] * \
                    np.pi/180.0  # convert to radians
                rotationmatrix = np.array(
                    [[np.cos(theta), -np.sin(theta), 0.0],
                     [np.sin(theta), np.cos(theta), 0.0],
                     [0.0, 0.0, 0.0]]
                )
                self.params['velocityvector'] = np.dot(
                    rotationmatrix, velocityvector).tolist()
            else:
                raise(
                    Exception("~~~~~ ERROR: could not construct wind velocity vector ~~~~~"))
        if self.params['heading'] == None and not(self.params['velocityvector'] == None):
            self.params['heading'] = np.arctan2(
                self.params['velocityvector'][1], self.params['velocityvector'][0])
            if self.params['heading'] < 0.0:
                self.params['heading'] = self.params['heading'] + 2*np.pi
        if self.params['speednortheast'] == None and not(self.params['velocityvector'] == None):
            self.params['speednortheast'] = np.hypot(
                self.params['velocityvector'][1], self.params['velocityvector'][0])
        if self.params['downdraftspeed'] == None and not(self.params['velocityvector'] == None):
            self.params['downdraftspeed'] = self.params['velocityvector'][2]


print("Successfully imported `Wind` class")


class Gust(WeatherType):
    '''
    Class used to define gust characteristics. Params variables include:
    * variable (parameter) [units] {example}
    * amplitude (gust amplitude) [m/s]
    * frequency (gust frequency) [Hz]

    Gusts are expected to affect:
    * control effort
    * ?? overall minimal effect on power consumption ??
    '''

    def __init__(self, params):
        paramnames = ['amplitude', 'frequency']
        WeatherType.__init__(self, params, paramnames)


print("Successfully imported `Gust` class")


class Ice(WeatherType):
    '''
    Class used to define icing characteristics. Params variables include:
    * variable (parameter) [units] {example}

    Ice is expected to affect:
    * lift
    * drag
    NOTE: Ice is very dangerous and difficult to model
    '''

    def __init__(self, params):
        paramnames = ['']
        WeatherType.__init__(self, params, paramnames)
        print("Too complicated. Stop now while you can!")  # haha seriously!


print("Successfully imported `Ice` class")


class Mission:
    'Class used to define the mission, including flight trajectory, maneuvers, etc.'

    params = {}

    def __init__(self, params):
        self.params = params


class Simulation:
    'Class used to run simulations'

    # class variables go here:
    params = {
        # setup parameters
        'simulationtype': None,
        #'desiredresult': None,
        'timestep': None,        # in seconds
        'clock': 0.0,            # tracks the current time
        'counter': 0,            # tracks the iteration number (0-indexed)
        # results
        'range':None,
        'endurance':None
    }

    # methods go here:
    def __init__(self, timestep, simulationtype): #, desiredresult):
        # get parameters from `settings.txt`
        self.params['simulationtype'] = simulationtype
        # self.params['desiredresult'] = desiredresult
        # print("Desired Result is        ", self.params['desiredresult'])

    def run(self, drone, battery, power, weather, mission):
        if self.params['simulationtype'] == 'simple':
            # insert simple model here
            self.__runSimpleModel(drone, battery, power, weather, mission)
        # we'll need to define a list of these terms in the README
        elif self.params['simulationtype'] == 'complicated':
            # insert another model here
            pass
        else:
            # insert another model here
            raise(Exception("~~~~~ ERROR: simulation model not available ~~~~~"))

    def __runSimpleModel(self, drone, battery, power, weather, mission):
        self.params['endurance']    = battery.capacity * battery.voltagemean / \
                                      power.params['power']  # simple endurance model
        cruisespeed                 = mission.params["missionspeed"]
        self.params['range']        = self.params['endurance'] * cruisespeed

    def __model2(self, drone, battery, power, weather):
        print('Model 2 is still in development')

    def __model3(self, drone, battery, power, weather):
        print('Model 3 is still in development')

    # # time-invariant methods go here (indicated by suffix *_ti):
    # def getEndurance_ti(self, power):
    #         # check that objects are time-invariant
    #         if power.timevariant == True:
    #                 print("")
    #                 sys.exit("ERR: attempted to run a time-invariant simulation with a time-variant `Power` object")
    #         else:
    #                 # calculate endurance
    #                 return self.capacity * self.voltage_mean / power.params['power']


print("Successfully imported `Simulation` class")


class Plotter:
    'Class used to make plots of simulation results'

    # class variables go here:
    fig_num = 1

    # methods go here:
    # def __init__(self,x,xlabel: string,y,ylabel,axistitle) -> None:
    def __init__(self, x, xlabel, y, ylabel, axistitle):
        d = datetime.datetime.today()
        self.title = axistitle + " (" + d.strftime("%b-%d-%Y") + ")"
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.x = x
        self.y = y
        # plt.plot(self.x,self.y)
        # plt.xlabel(self.xlabel)
        # plt.ylabel(self.ylabel)
        # plt.title(self.title)
        # plt.show()

    def plot_line(self):
        fig = plt.figure(Plotter.fig_num)
        fig.patch.set_facecolor('w')
        plt.plot(self.x, self.y)
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        plt.title(self.title)
        Plotter.fig_num += 1
        fig.show()
        input()

    def plot_scatter(self):
        fig = plt.figure(Plotter.fig_num)
        fig.patch.set_facecolor('w')
        plt.plot(self.x, self.y, 'ro')
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        plt.title(self.title)
        Plotter.fig_num += 1
        fig.show()
        input()

    def plot_validation(self, xvalid, yvalid):
        fig = plt.figure(Plotter.fig_num)
        fig.patch.set_facecolor('w')
        plt.plot(self.x, self.y)
        plt.plot(xvalid, yvalid, 'ko')
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        plt.title(self.title)
        Plotter.fig_num += 1
        fig.show()
        input()

print("Successfully imported `Plotter` class")
