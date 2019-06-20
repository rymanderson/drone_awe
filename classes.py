import datetime
import matplotlib.pyplot as plt
import numpy as np

#insert classes here

class Drone:
        'Class used to store key drone characteristics'

        name            = None
        params          = None
        conversions     = None
        correctunits    = False

        payload         = 0.0

        def __init__(self, name, params, conversions):
                self.name       = name
                self.params     = params
                self.conversions= conversions # this input is obtained using functions.getParams() in exe.py
                self.__convertUnits()
                self.__getEfficiencyPropulsive()
        
        def updatePayload(self,payload)
                self.payload = payload

        def __convertUnits(self):
                if not self.correctunits:
                        for spec in self.conversions:
                                self.params[spec] = self.params[spec] * self.conversions[spec]
                        self.correctunits    = True
                else:
                        print("~~~~~ WARNING: problem converting units: skipping ~~~~~")

        def __getEfficiencyPropulsive(self):
                thrust                  = self.params['takeoffweight']
                batteryenergy           = self.params['batteryenergy']
                endurancemaxhover       = self.params['endurancemaxhover']
                rotorarea               = self.params['rotordiameter']**2/4*np.pi
                airdensity              = 1.225                                                # assuming air density was equal to 1.225 kg/m3 during drone testing
                poweractual             = batteryenergy/endurancemaxhover
                powerideal              = thrust * np.sqrt(thrust/(2*rotorarea*airdensity))
                efficiency              = powerideal/poweractual
                # prediction based on momentum theory for hover case
                # slides from https://fenix.tecnico.ulisboa.pt/downloadFile/282093452028191/3-Momentum%20Theory%20in%20hover.pdf

                self.params['efficiencypropulsive'] = efficiency

print("Successfully imported `Drone` class")


class Battery:
        'Class used to track battery characteristics and performance during simulation'

        # class variables go here:
        capacity        = None          # Ampere-hours
        soc             = None          # state of charge (in percent nominal capacity)
        startsoc        = None          # state of charge at simulation start
        soh             = None          # state of health (actual capacity divided by ideal capacity)
        startsoh        = None          # state of health at simulation start (for running many simulations with the same drone)
        batterytype     = None          # possible values include LiPo, Li-ion, NiCd, NiMH, SLA
        voltage         = None          # Volts; this is the instantaneous voltage
        voltagemean     = None          # Volts; this is the average voltage used for time-invariant simulations
        voltagecharged  = None          # Volts
        voltagedead     = None          # Volts
        current         = None          # Amperes; this is the instantaneous current

        # constructor
        def __init__(self, drone, soh = 90.0, startsoc = 100.0): # default value for soh is based on the assumption that batteries are retired at a soh of 80%
                # import parameters from drone object
                self.batterytype = drone.params['batterytype']
                self.voltagemean = drone.params['batteryvoltage']
                self.capacity    = drone.params['batterycapacity']

                # update parameters
                self.updateBattery()

                # estimate list lengths for prior memory allocation

        def updateBattery(self):
                print("still working on Battery.updateBattery method")

        # time-variant methods go here:
        def updateLoad(self, power):
                self.current    = power.power / self.voltage

        def discharge(self, power, timestep):
                self.soc        = self.soc - power.power * timestep
                self.current    = power.power / self.voltage
                

print("Successfully imported `Battery` class")


class Power:
        'Class used to predict the drone\'s power requirement'

        # class variables go here:
        power           = None      # watts

        # methods go here:
        def __init__(self,params,drone,weather):
                # self.density    = weather.density
                # self.mass       = drone.params[TOW]
                # if drone.params['wingtype'] == 'rotary':
                #         self.Nrotors    = drone.params[num_rotors]
                #         self.L_D        = 3 #may need to change
                #         Power.getPowerRotor1()
                # elif drone.params['wingtype'] == 'fixed':
                #         self.L_D        = 10 #may need to change
                #         self.power      = 0
                pass # what parameters does the power class really need to store? I don't like double-storing, but it might increase speed, in which case it could be worth it

        def updatePower(self,drone,weather,model):
                if model == 'dandrea':
                        __getPowerDandrea(drone,weather)
                elif model == 'abdilla':
                        __getPowerAbdilla(drone,weather)

        def __getPowerDandrea(self,drone,weather): #super simple estimate for power from D'Andrea `Can Drones Deliver`
                powerelectronics        = 0.1          # kW, estimate from paper
                L_D                     = 3.0          # quick estimate for initial functionality
                self.power              = (drone.param['takeoffweight'] + drone.payload) * drone.param['endurancemaxspeed'] / (370.0 * drone.params['efficiencypropulsive'] * L_D) - powerelectronics

        def __getPowerAbdilla(self,drone,weather): #slightly more complicated estimate for power
                self.power = (drone.param['takeoffweight'] + drone.payload) / (drone.params['efficiencypropulsive'] * drone.params['rotordiameter'] / 2.0) * np.sqrt(weather.params['gravitationconstant']**3 / (2 * drone.params['rotorquantity'] * weather.params['airdensity'] * np.pi))

print("Successfully imported `Power` class")

class Weather:
        'Class describing ambient weather conditions and is used to predict the drone\'s power requirement'

        # miscellaneous weather   
        weatherlist     = [] # a list of weather objects (e.g., rain, icing, humidity, etc.)

        # ambient air quality
        params          = {
                        'airdensity':1.225,
                        'gravitationconstant':9.807
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

        #need to figure out how to access current or at least average weather
                #data based on lattitude, longitude, altitude?
        #then calculate pressure, temp, and density based on altitude (see
                # eqns 1.21 - 1.22 in Dr. Ning's book
        #or do above with average values of temperature/pressure/density for a given location
        #or none of the above - we vary temperature and assume either constant pressure or density (likely pressure)

        # methods go here:
        def __init__(self, altitude,temperature_sl): # keeping it simple to begin with
                self.altitude = altitude       
                self.temperature_sl = temperature_sl 
                self.temperature = self.temperature_sl - 71.5 + 2*np.log(1 + np.exp(35.75 - 3.25*self.altitude) + np.exp(-3 + 0.0003 * self.altitude**3))
                self.pressure = self.pressure_sl * np.exp(-0.118 * self.altitude - (0.0015*self.altitude**2) / (1 - 0.018*self.altitude + 0.0011 * self.altitude**2))

        def calculateAtmosphere(self,altitude)
                temperaturesealevel             = 288.15               # Kelvin
                pressuresealevel                = 1.01325e5            # Pascals
                gravitationconstantsealevel     = 9.80665              # m/s2
                specificheatratio               = 1.4
                airgasconstant                  = 287.053              # J/(kg-K)
S = 110.4; %K
beta = 1.458e-6; %kg/(smK^1/2)

p = psl * exp( -0.118 * h - (0.0015 * h.^2)./(1 - 0.018 * h + 0.0011 * h.^2) );

T = Tsl - 71.5 + 2 * log( 1 + exp(35.75 - 3.25 * h) + exp(-3 + 0.0003 * h.^3) );

rho = p ./ ( R .* T );

a = sqrt( gamma .* R .* T );

mu = beta * T.^(3/2) ./ ( T + S );

        def update_temp(self,new_temp):
                test = PowerCorrection('temp',new_temp)

        def update_wind(self,velocity,heading):
                wind_vars = [velocity, heading]
                test = PowerCorrection('temp',wind_vars)

        def update_rain(self,LWC):
                test = PowerCorrection('temp',LWC)

        def update_humidity(self,rel_hum):
                test = PowerCorrection('temp',rel_hum)

        def update_icing(self):
                test = PowerCorrection('temp',new_temp) 

print("Successfully imported `Weather` class")

# Add other weather classes: rain, temp, humid, wind, etc.

class Rain:
        'Class used to define rain characteristics'
        #class variables go here:
        LWC             = None
        dropsize        = None
        WVC             = None

print("Successfully imported `Rain` class")


class Temperature:
        'Class used to define temperature characteristics'
        #class variables go here:
        temperature = None

        def __init__(self,temp):
                self.temperature = temp

print("Successfully imported `Temperature` class")


class Humidity:
        'Class used to define humidity characteristics'
        #class variables go here:
    
        def __init__(self,rel_hum):
                self.rel_hum = rel_hum

print("Successfully imported `Humidity` class")


class Wind:
        'Class used to define wind characteristics'
        #class variables go here:
        speed   = None
        heading = 0 #degrees from north maybe? We could simply use a vector in North-east-down directions. Do we want to worry about up/downdrafts?

        def __init__(self,speed,heading):
                self.speed = speed
                self.direction = direction

print("Successfully imported `Wind` class")


class Ice:
        'Class used to define icing characteristics'
        #class variables go here:
        def __init__(self):
                print("Too complicated. Stop now while you can!") # Hah! I laugh in the face of impending complexity and assured coding madness.

print("Successfully imported `Rain` class")



class Simulation:
        'Class used to run simulations'

        # class variables go here:
        params = {
                'timestep':None,        # in seconds
                'clock':0.0,            # tracks the current time
                'counter':0,            # tracks the iteration number (0-indexed)
                }

        #Notes from call with Ryan:
                # Simulation init sets up variables
                # run.Simulation would have if statements based on a simulationtype variable that would then call other class methods. It would also have class instances as inputs, which would then be inputs to the called class methods.

        # methods go here:
        def __init__(self,timestep,simulationtype,desiredresult):
                # get parameters from `settings.txt`
                print("still working on simulation class constructor")
                self.simulationtype = simulationtype
                self.desiredresult = desiredresult

        def run(self,drone,battery,power,weather):
                if self.simulationtype == 'simple':
                        #insert simple model here
                        self.runSimpleModel(drone,battery,power,weather)
                elif self.simulationmodel == 'complicated': # we'll need to define a list of these terms in the README
                        #insert another model here
                        pass 
                else:
                        #insert another model here
                        pass

        def runSimpleModel(self,drone,battery,power,weather):
                if self.desiredresult == 'Endurance' or 'endurance':
                        return battery.capacity * battery.voltagemean / power.power # simple endurance model
                elif self.desiredresult == 'Range' or 'range':
                        return drone.params['cruisespeed'] * battery.capacity * battery.voltagemean / power.power # multiplies endurance by cruise speed to get range

        def model2(self,drone, battery,power,weather):
                print('Model 2 is still in development')
        
        def model3(self,drone, battery,power,weather):
                print('Model 3 is still in development')

        # # time-invariant methods go here (indicated by suffix *_ti):
        # def getEndurance_ti(self, power):
        #         # check that objects are time-invariant
        #         if power.timevariant == True:
        #                 print("")
        #                 sys.exit("ERR: attempted to run a time-invariant simulation with a time-variant `Power` object")
        #         else:
        #                 # calculate endurance
        #                 return self.capacity * self.voltage_mean / power.power


print("Successfully imported `Simulation` class")


class Plotter:
	'Class used to make plots of simulation results'

	# class variables go here:
	fig_num = 1

	# methods go here:
	def __init__(self,x,xlabel,y,ylabel,axistitle):
		d = datetime.datetime.today()
		self.title      = axistitle + " (" + d.strftime("%b-%d-%Y") + ")"
		self.xlabel     = xlabel
		self.ylabel     = ylabel
		self.x          = x
		self.y          = y
		# plt.plot(self.x,self.y)
		# plt.xlabel(self.xlabel)
		# plt.ylabel(self.ylabel)
		# plt.title(self.title)
		# plt.show()

	def plot_line(self):
		fig = plt.figure(Plotter.fig_num)
		fig.patch.set_facecolor('w')
		plt.plot(self.x,self.y)
		plt.xlabel(self.xlabel)
		plt.ylabel(self.ylabel)
		plt.title(self.title)
		Plotter.fig_num += 1
		fig.show()
		input()

	def plot_scatter(self):
		fig = plt.figure(Plotter.fig_num)
		fig.patch.set_facecolor('w')
		plt.plot(self.x,self.y,'ro')
		plt.xlabel(self.xlabel)
		plt.ylabel(self.ylabel)
		plt.title(self.title)
		Plotter.fig_num += 1
		fig.show()
		input()

print("Successfully imported `Plotter` class")
