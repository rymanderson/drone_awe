import datetime
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve
import functions as fun
import gekko

# insert classes here


class Drone:
    'Class used to store key drone characteristics'

    name            = None
    params          = None
    conversions     = None
    correctunits    = False

    def __init__(self, name, params, conversions):
        self.name                   = name
        self.params                 = params
        # this input is obtained using functions.getParams() in exe.py
        self.conversions            = conversions
        self.__convertUnits()
        # extrapolate useful parameters
        if 'width' in self.params and 'height' in self.params:
            self.params['frontalarea']  = self.params['width'] * self.params['height']
        else:
            self.params['frontalarea']  = 1.0
            print("Drone.__init__:  'width' and 'height' not found; 'frontalarea' set to 1.0")
            # NOTE: if width and height are not available AND a test run is not available to calibrate C_D, method will not work
        if 'width' in self.params and 'length' in self.params:
            self.params['toparea'] = self.params['width'] * self.params['length']
        else:
            self.params['toparea'] = 8.0
            print("Drone.__init__:  'width' and 'length' not found; 'frontalarea' set to 1.0")
        if 'rangemax' in self.params and 'rangemaxspeed' in self.params:
            self.params['endurancemaxrange'] = self.params['rangemax'] / self.params['rangemaxspeed']
        if 'rotordiameter' in self.params:
            self.params['rotorarea'] = self.params['rotordiameter']**2/4*np.pi # area per rotor
        else:
            raise(Exception("~~~~~ ERROR: rotor diameter not found ~~~~~"))
        if self.params['wingtype'] == 'fixed':
            self.params['spanefficiency'] = 0.8 #estimate from Dr. Ning's book (he lists 0.7-0.9). If we want to we could decrease this further based on fuselage diameter, but maybe that's requiring too much detail
        if 'numbatteries' in self.params:
            if self.params['numbatteriesconnection'] == 'parallel':
                self.params['batterycapacity'] *= self.params['numbatteries'] #increase capacity in parallel
            elif self.params['numbatteriesconnection'] == 'series':
                self.params['batteryvoltage'] *= self.params['nummbatteries'] #increase voltage in series
            else:
                 raise(Exception("~~~~~ ERROR: incorrect battery connection parameter applied ~~~~~"))
        # print parameters for debugging
        self.__printParameters()

    def __convertUnits(self):
        if not self.correctunits:
            for spec in self.conversions:
                if spec in self.params:
                    self.params[spec] = self.params[spec] * \
                        self.conversions[spec]
                else:
                    print("~~~~~ WARNING: drone parameter ", spec, " not found- skipping ~~~~~")
            self.correctunits = True
        else:
            print("~~~~~ WARNING: problem converting units: skipping ~~~~~")

    def __printParameters(self):
        print("DRONEIMPORT:     takeoffweight is: ",self.params['takeoffweight'])

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
        self.current = power.params['power'] / self.params['voltage']

    def discharge(self, power, timestep):
        self.soc = self.soc - power.params['power'] * timestep
        self.current = power.params['power'] / self.params['voltage']


print("Successfully imported `Battery` class")


class Power:
    'Class used to predict the drone\'s power requirement'

    # class variables go here:
    params = {
        'efficiencypropulsive': None,
        'power':None,
        'model':None,
        'dragcoefficient':None,
        'alpha':None
    }

    # methods go here:
    def __init__(self, drone, weather, mission):
        # initial propulsive efficiency values

        self.update(drone, weather, mission)


    def update(self, drone, weather, mission):
        self.__setParameters(drone,weather)
        self.__getDragCoefficient(drone)
        self.__updateEfficiencyPropulsive(drone, weather, mission)
        self.__getAlphaVelocity(drone,weather,mission)
        self.__getDrag(drone,weather,mission)
        self.__getThrust(drone,weather,mission)
        
        if drone.params['wingtype'] == 'rotary':
            # self.__getPowerAbdilla(drone, weather, mission)
            self.__getPowerMomentum(drone, weather, mission)
        elif drone.params['wingtype'] == 'fixed':
            self.__getPowerTraub(drone, weather, mission)
        else:
            # raise Exception(f"~~~~~ ERROR: model { model } not available ~~~~~") #
            raise Exception("~~~~~ ERROR: model '" +
                            self.params['model'] + "' not available ~~~~~")
        # print("power.update(): power is            ",self.params['power'])
        # print("power.update(): drag coefficient is ",self.params['dragcoefficient'])

    # def update(self, drone, weather, mission):
    #     self.__getDragCoefficient(drone)
    #     self.__updateEfficiencyPropulsive(drone, mission)
    #     if self.params['model'] == 'dandrea':
    #         self.__getPowerDandrea(drone, weather)
    #     elif self.params['model'] == 'abdilla':
    #         self.__getPowerAbdilla(drone, weather, mission)
    #     else:
    #         # raise Exception(f"~~~~~ ERROR: model { model } not available ~~~~~") #
    #         raise Exception("~~~~~ ERROR: model '" +
    #                         self.params['model'] + "' not available ~~~~~")

    # # super simple estimate for power from D'Andrea `Can Drones Deliver` *****Doesn't work well*******
    # def __getPowerDandrea(self, drone, weather):
    #     powerelectronics = 0.1          # kW, estimate from paper
    #     L_D = 3.0          # quick estimate for initial functionality TODO: Change this to something more scientific
    #     self.params['power'] = (drone.params['takeoffweight'] + drone.params['payload']) * drone.params['endurancemaxspeed'] / (
    #         370.0 * self.params['efficiencypropulsive'] * L_D) - powerelectronics

    # slightly more complicated estimate for power
    def __getPowerAbdilla(self, drone, weather, mission):
        self.params['power'] = self.params['thrust']**1.5 / \
                               self.params['efficiencypropulsive'] / \
                                np.sqrt(2 * drone.params['rotorquantity'] * \
                                weather.params['airdensity'] * drone.params['rotorarea']) + \
                                0.0  # Camera power consumption estimate

        self.__printParameters(drone,weather,mission)

    def __getPowerMomentum(self, drone, weather, mission):
        self.__getBladeProfilePower()
        self.params['power'] = self.params['thrust']*self.params['velocityinduced'] + \
                                                self.params['drag']*mission.params['missionspeed'] + \
                                                self.params['bladeprofilepower']
        self.__printParameters(drone,weather,mission)

    def __getPowerTraub(self, drone, weather, mission): #fixed-wing power model
        density = weather.params['airdensity']
        cruisespeed = mission.params['missionspeed']
        if 'wingarea' in drone.params:
            wingarea = drone.params['wingarea']
        elif 'wingspan' in drone.params and 'chord' in drone.params:
            wingarea = drone.params['wingspan'] * drone.params['chord']
        else:
            raise(Exception("~~~~~ ERROR: wing area needed to calculate power ~~~~~"))
        if 'wingspan' not in drone.params:
            raise(Exception("~~~~~ ERROR: wing span needed to calculate power ~~~~~"))
        else:
            span = drone.params['wingspan']
        self.__getDragCoefficient(drone)
        dragcoefficient = self.params['dragcoefficient']
        weight = drone.params['takeoffweight'] * weather.params['gravitationconstant']
        spanefficiency = drone.params['spanefficiency']
        k = 1 / (np.pi*span**2 / wingarea * spanefficiency)

        self.params['power'] = 0.5 * density * cruisespeed**3 * wingarea * dragcoefficient \
                                + 2*weight**2*k \
                                / (density*cruisespeed*wingarea)

    def __updateEfficiencyPropulsive(self, drone, weather, mission):
        # default value:
        if 'endurancemax' not in drone.params or 'endurancemaxspeed' not in drone.params:
            self.params['efficiencypropulsive'] = 0.5
        else:
            # get efficiency at max endurance conditions
            etamaxendurance = self.__getEfficiencyPropulsive(
                drone, weather, drone.params['endurancemax'])
            vmaxendurance = drone.params['endurancemaxspeed']
            # get efficiency at max range conditions
            etamaxrange = self.__getEfficiencyPropulsive(
                drone, weather, drone.params['rangemax']/drone.params['rangemaxspeed'])
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
        
        # print("Drone:       efficiencypropulsive is ",self.params['efficiencypropulsive'])

    def __getEfficiencyPropulsive(self, drone, weather, endurance):
        # Analyzing a single propeller
        thrust          = drone.params['takeoffweight'] * weather.params['gravitationconstant'] / drone.params['rotorquantity']
        batteryenergy   = drone.params['batteryenergy'] / \
                          drone.params['rotorquantity']
        rotorarea       = drone.params['rotorarea']
        airdensity      = 1.225  # assuming air density was equal to 1.225 kg/m3 during drone testing

        poweractual     = batteryenergy/endurance
        powerideal      = thrust * np.sqrt(thrust/(2*rotorarea*airdensity))
        efficiency      = powerideal/poweractual
        # prediction based on momentum theory for hover case
        # slides from https://fenix.tecnico.ulisboa.pt/downloadFile/282093452028191/3-Momentum%20Theory%20in%20hover.pdf
        # prediction based on momentum theory for forward flight
        # slides from https://fenix.tecnico.ulisboa.pt/downloadFile/845043405440064/6-Momentum%20Theory%20in%20forward%20flight.pdf
        # see also: http://www.aerospaceweb.org/design/helicopter/momentum.shtml

        return efficiency

    def __getAlphaVelocity(self,drone,weather,mission):

        velocityinfinity    = mission.params['missionspeed']
        totalweight         = drone.params['totalweight']
        airdensity          = weather.params['airdensity']
        rotorarea           = drone.params['rotorarea']
        rotorquantity       = drone.params['rotorquantity']
        dragcoefficient     = self.params['dragcoefficient']
        frontalarea         = drone.params['frontalarea']
        toparea             = drone.params['toparea']
        pi                  = np.pi

        # def momentumTheoryEquations(variables):
        #     alpha, velocityinduced  = variables
            
        #     return1     = velocityinduced**4 + velocityinduced**3 * (2*velocityinfinity*np.sin(alpha*np.pi/180.0)) + \
        #                 velocityinduced**2 * velocityinfinity**2 - (totalweight/(2*airdensity*rotorarea*rotorquantity))**2
        #     return2     = totalweight/(rotorquantity*np.cos(alpha*np.pi/180.0)) - \
        #                 airdensity*velocityinfinity**2*dragcoefficient / (2*rotorquantity*np.sin(alpha*np.pi/180.0)) * \
        #                 toparea*np.sin(alpha*np.pi/180.0) + frontalarea*np.cos(alpha*np.pi/180.0)

        #     return (return1,return2)
        
        # (alpha, velocityinduced)        = fsolve(momentumTheoryEquations,(0.0,0.0))

        m               = gekko.GEKKO(remote=False)             # create GEKKO model
        alpha           = m.Var(value=0.0)      # define new variable, initial value=0
        velocityinduced = m.Var(value=self.params['velocityinducedhover'])      # define new variable, initial value

        # print values for debugging
        print("GEKKO:       airdensity*velocityinfinity**2*pi*toparea = ",airdensity*velocityinfinity**2*pi*0.006)

        qty             = (totalweight/rotorquantity + 0.00000000/rotorquantity)**2

        m.Equations([ \


# airdensity*velocityinfinity**2*pi*alpha*0.006


            2*airdensity*rotorarea*velocityinduced*m.sqrt(velocityinfinity**2 + 2*velocityinfinity*velocityinduced*m.sin(alpha) + velocityinduced**2) == \
            m.sqrt( \
            (totalweight/rotorquantity + 0/rotorquantity)**2 + \
            # qty + \
            (1/2*airdensity*velocityinfinity**2*dragcoefficient * (toparea*m.sin(-alpha) + frontalarea*m.cos(alpha)))**2 ), \
            velocityinduced == totalweight/rotorquantity/2.0/airdensity/rotorarea/m.sqrt((velocityinfinity*m.cos(alpha))**2 + (velocityinfinity*m.sin(alpha) + velocityinduced)**2)

            # velocityinduced**4 + velocityinduced**3 * (2*velocityinfinity*m.sin(alpha*np.pi/180.0)) + \
            # velocityinduced**2 * velocityinfinity**2 - (totalweight/(2*airdensity*rotorarea*rotorquantity))**2==0.0, \
            # totalweight/(rotorquantity*m.cos(alpha*np.pi/180.0)) - \
            # airdensity*velocityinfinity**2*dragcoefficient / (2*rotorquantity*m.sin(alpha*np.pi/180.0)) * \
            # toparea*m.sin(alpha*np.pi/180.0) + frontalarea*m.cos(alpha*np.pi/180.0) == 0.0 \

            ]) # equations
        m.solve(disp=False)     # solve

        alpha                           = alpha.value[0]
        # ensure alpha is between -pi and pi
        if alpha > np.pi or alpha < -np.pi:
            print("NOTE: alpha was evaluated as ",alpha,"; adjusting")
            while alpha > np.pi or alpha < -np.pi:
                if alpha > np.pi:
                    alpha = alpha - 2*np.pi
                elif alpha < -np.pi:
                    alpha = alpha + 2*np.pi
        self.params['alpha']            = alpha            
        self.params['velocityinduced']  = velocityinduced.value[0]

    def __getDrag(self,drone,weather,mission):
        drag = 0.5 * self.params['dragcoefficient'] * weather.params['airdensity'] * mission.params['missionspeed']**2 * \
               (drone.params['frontalarea'] * np.cos(self.params['alpha']*np.pi/180.0) + \
                drone.params['toparea'] * np.sin(self.params['alpha']*np.pi/180.0))
        self.params['drag'] = drag

    def __getThrust(self,drone,weather,mission):
        totalweight = drone.params['totalweight']
        drag        = self.params['drag']
        thrust      = np.sqrt(totalweight**2 + drag**2)

        self.params['thrust'] = thrust

    def __getDragCoefficient(self,drone):
        self.params['dragcoefficient'] = 0.07
    
    def __getBladeProfilePower(self):
        self.params['bladeprofilepower'] = 145.0

    def __setParameters(self,drone,weather):
        drone.params['totalweight']             = (drone.params['takeoffweight'] + drone.params['payload']) * \
                                                  weather.params['gravitationconstant']
        self.params['area']                     = np.pi * drone.params['rotordiameter']**2/4
        self.params['velocityinducedhover']     = np.sqrt(drone.params['totalweight'] / \
                                                  drone.params['rotorquantity'] / 2 / \
                                                  self.params['area']/weather.params['airdensity'])

    def __printParameters(self,drone,weather,mission):
        velocityinfinity    = mission.params['missionspeed']
        totalweight         = drone.params['totalweight']
        airdensity          = weather.params['airdensity']
        rotorarea           = drone.params['rotorarea']
        rotorquantity       = drone.params['rotorquantity']
        dragcoefficient     = self.params['dragcoefficient']
        frontalarea         = drone.params['frontalarea']
        toparea             = drone.params['toparea']
        velocityinduced     = self.params['velocityinduced']
        alpha               = self.params['alpha']

        thrust2 = 2*airdensity*rotorarea*velocityinduced * \
                np.sqrt(velocityinfinity**2 + 2*velocityinfinity*velocityinduced*np.sin(alpha) + velocityinduced**2)
        thrust3 = np.sqrt((totalweight/rotorquantity)**2 + (1/2*airdensity*velocityinfinity**2 * \
                dragcoefficient * (toparea*np.sin(alpha) + frontalarea*np.cos(alpha)))**2)
        velocityinduced2 = totalweight/rotorquantity/2.0/airdensity/rotorarea/np.sqrt((velocityinfinity*np.cos(alpha))**2 + \
                           (velocityinfinity*np.sin(alpha) + velocityinduced)**2)

        print("----- MISSION SPEED = ",mission.params['missionspeed']," -----")
        print("")
        print("power:           alpha is                ",self.params['alpha']*180/np.pi)
        print("power:           thrust is               ",self.params['thrust'])
        print("power:           thrust2 is              ",thrust2)
        print("power:           thrust3 is              ",thrust3)
        print("power:           drag is                 ",self.params['drag'])
        print("power:           velocityinduced is      ",self.params['velocityinduced'])
        print("power:           velocityinduced2 is     ",velocityinduced2)
        print("power:           velocityinducedhover is ",self.params['velocityinducedhover'])
        print("power:           etapropulsive is        ",self.params['efficiencypropulsive'])
        print("")
        print("drone:           totalweight is          ",drone.params['totalweight'])
        print("drone:           takeoffweight is        ",drone.params['takeoffweight'])
        print("drone:           payload is              ",drone.params['payload'])
        print("drone:           rotorquantity is        ",drone.params['rotorquantity'])
        print("")
        print("weather:         airdensity is           ",weather.params['airdensity'])
        print("weather:         gravitationconstant is  ",weather.params['gravitationconstant'])
        print("")  
        print("")

print("Successfully imported `Power` class")


class Weather:
    '''
    Class describing ambient weather conditions and is used to predict the drone\'s power requirement
    '''

    # miscellaneous weather
    # a list of weather objects (e.g., rain, icing, humidity, etc.)
    weatherlist = None

    # ambient air quality
    params = {
        'airdensity': 1.225,
        'airdensitysealevel': 1.225,
        'gravitationconstant': 9.807,
        'temperaturesealevel': 288.15, #15 deg C
        'temperature': None,
        'humidity': None,
        'altitude':None
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
    def __init__(self, altitude, weatherlist):  # keeping it simple to begin with
        self.weatherlist = weatherlist
        self.altitude = altitude
        # self.temperature = temperaturesealevel
        # self.temperature = self.temperaturesealevel - 71.5 + 2*np.log(1 + np.exp(35.75 - 3.25*self.altitude) + np.exp(-3 + 0.0003 * self.altitude**3))
        # self.pressure = self.pressure_sl * np.exp(-0.118 * self.altitude - (0.0015*self.altitude**2) / (1 - 0.018*self.altitude + 0.0011 * self.altitude**2))

    def update(self):
        # update independent parameters
        # print("Weatherlist:", self.weatherlist)
        for weathertype in self.weatherlist:
            for param in weathertype.params:
                self.params[param] = weathertype.params[param]
        # update dependent parameters
        densityfactor = []
        for weathertype in self.weatherlist: 
            densityfactor.append(weathertype.updateDensity(self))
        airdensity = self.params['airdensitysealevel']
        for x in densityfactor:
            airdensity *= x
        self.params['airdensity'] = airdensity 
        # since weatherlist does not contain strings, we may need to 
        # include a method for density in each weather class. For those 
        # that do nothing to density, we simply return densityfactor = 1.0

        # print('Need to add more here for wind, rain, etc.')
        #affect efficiency, missionspeed (wind), lift, weight, etc.

    def getStandardAtmosphere(self, altitude):
        '''This function currently assumes STP conditions at sea level and should probably be adjusted to use ground level conditions as a baseline'''

        # set sea level parameters
        temperaturesealevel             = self.params['temperaturesealvel'] # Kelvin
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

    Temperature is expected to affect:
    * air density
    * battery performance
    '''

    def __init__(self, params):
        paramnames = ['temperature']
        WeatherType.__init__(self, params, paramnames)
    
    def updateDensity(self,weather):
        newtemperature = weather.weatherlist[0].params['temperature'] +273.15 #convert to kelvin
        oldtemperature = weather.params['temperaturesealevel'] #convert to kelvin
        # print("old temperature is",oldtemperature)
        # print("new temperature is",newtemperature)
        # print("Self.params['temperature'] is",self.params['temperature'])
        olddensity = weather.params['airdensity']
        # print("old density is",olddensity)
        newdensity = olddensity * ((oldtemperature-newtemperature)/oldtemperature + 1) #inverse relationship
        temperatureeffect = newdensity / olddensity
        return temperatureeffect


print("Successfully imported `Temperature` class")


class Humidity(WeatherType):
    '''
    Class used to define humidity characteristics. Params variables include:
    * variable (parameter) [units] {example}
    * relativehumidity (relative humidity) [%]
    * absolutehumidity (absolute humidity) []

    Humidity is expected to affect:
    * air density
    * condensation
    '''

    def __init__(self, params):
        paramnames = ['relativehumidity']    #, 'absolutehumidity'] - most of the time it is only given in relative terms
        WeatherType.__init__(self, params, paramnames)

    def updateDensity(self,weather): # 2D linear interpolation based on Yue (2017)
        relativehumiditylist = [0.0,25.0,50.0,70.0,90.0] #percentage
        temperaturelist = [15.0,20.0,25.0,30.0,35.0] #may need to convert this to kelvin

        #humidity effects only
        humidityarray = np.array([[1.0, 1.0, 1.0, 1.0, 1.0],\
                                    [0.99837, 0.99751, 0.99662, 0.99656, 0.99389],\
                                    [0.99673, 0.99585, 0.99409, 0.99055, 0.98866],\
                                    [0.99592, 0.99419, 0.99155, 0.98883, 0.98517],\
                                    [0.99429, 0.99252, 0.98902, 0.98625, 0.98080]]) # increasing humditiy in rows, temp in columns


        # temperature and humidity effects included
        # humidityarray = np.array([[1.0, 0.98286, 0.96653, 0.95020, 0.93551],\
        #                             [0.99837, 0.98041, 0.96327, 0.94694, 0.92980],\
        #                             [0.99673, 0.97878, 0.96082, 0.94122, 0.92490],\
        #                             [0.99592, 0.97714, 0.95837, 0.93959, 0.92163],\
        #                             [0.99429, 0.97551, 0.95592, 0.93714, 0.91755]]) # increasing humditiy in rows, temp in columns

        relativehumidity = weather.weatherlist[1].params['relativehumidity']
        # temperature = weather.params['temperaturesealevel']
        # temperature = 288.15 - 273.15 #convert to celcius - need to update this if temperature is altered
        temperature = weather.weatherlist[0].params['temperature'] #already to celcius
        # print("Temperature is",temperature)
        # determine 4 closest points before doing a 2D interpolation

        if temperature < temperaturelist[0]:
            temperature = temperaturelist[0]
            print("WARNING: Temperature is below interpolation bounds. Humidy effects are assumed to be as if the temperature is",temperaturelist[0],"degrees C")
        elif temperature > temperaturelist[-1]:
            temperature = temperaturelist[-1]
            print("WARNING: Temperature is below interpolation bounds. Humidy effects are assumed to be as if the temperature is",temperaturelist[0],"degrees C")
        
        if relativehumidity < 0.0 or relativehumidity > 100.0:
            raise(Exception("~~~~~ ERROR: relative humidity is outside of available bounds ~~~~~"))
        elif relativehumidity > relativehumiditylist[-1]:
            relativehumidity = relativehumiditylist[-1]
            print("WARNING: Humidity is above interpolation bounds. Humidy effects are assumed to be as if the temperature is",humiditylist[-1],"%")

        counter = 0
        for temp in temperaturelist:
            if temperature == temp:
                idT1 = counter
                idT2 = counter
                break
            elif temperature < temp:
                idT1 = counter - 1
                idT2 = counter
                break
            else:
                counter += 1
        counter = 0
        for hum in relativehumiditylist:
            if relativehumidity == hum:
                idH1 = counter
                idH2 = counter
                break
            elif relativehumidity < hum:
                idH1 = counter - 1
                idH2 = counter
                break
            else:
                counter += 1

        #interpolate temperature at each humidity index
        x1 = temperaturelist[idT1]
        x2 = temperaturelist[idT2]
        y1 = humidityarray[idH1,idT1]
        y2 = humidityarray[idH1,idT2]
        x = temperature
        tempid1 = fun.interpolate(x1,x2,y1,y2,x)
        
        x1 = temperaturelist[idT1]
        x2 = temperaturelist[idT2]
        y1 = humidityarray[idH2,idT1]
        y2 = humidityarray[idH2,idT2]
        x = temperature
        tempid2 = fun.interpolate(x1,x2,y1,y2,x)

        #now interpolate based on humidity
        x1 = relativehumiditylist[idH1]
        x2 = relativehumiditylist[idH2]
        y1 = tempid1
        y2 = tempid2
        x = relativehumidity
        humidityeffect = fun.interpolate(x1,x2,y1,y2,x)
        return humidityeffect


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
