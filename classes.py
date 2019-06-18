import datetime
import matplotlib.pyplot as plt

#insert classes here

class Drone:
        'Class used to store key drone characteristics'

        name    = None
        params  = None

        def __init__(self, name, params):
                self.name       = name
                self.params     = params

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
                self.batterytype = drone.params['batt_type']
                self.voltagemean = drone.params['batt_voltage']
                self.capacity    = drone.params['batt_capacity']

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
        basepower           = 0.0       # watts
        mass                = 1         # kg
        propulsive_eta      = 0.8       # propulsive efficiency
        radius              = None      # m
        g                   = 9.81      # m/s2
        Nr                  = 4         # number of rotors
        density             = None      # kg/m3
        power               = None      # watts

    # methods go here:
        def __init__(self,name,drone,weather):
                self.name = name
                # drone = Drone(name)
                drone.getParams()
                # self.Op1 = Weather.Weather(altitude)
                self.density = weather.density
                self.mass = drone.params[TOW]
                if drone.params['wingtype'] == 'rotary':
                        self.Nr = drone.params[num_rotors]
                        self.L_D = 3 #may need to change
                        Power.getpower_rotor1()
                elif drone.params['wingtype'] == 'fixed':
                        self.L_D = 10 #may need to change
                        self.power = 0

        def getpower_rotor1(self): #super simple estimate for power
                power_elec = 0.1 #kW, estimate from paper
                self.power = self.mass * self.drone.param[cruise_speed] / (370 * self.propulsive_eta * self.L_D) - power_elec

        def getpower_rotor2(self): #slightly more complicated estimate for power
                self.power = self.mass * 1 / (self.propulsive_eta * self.drone.params[prop_diameter]/2.0) * np.sqrt(self.g**3 / (2 * self.Nr * self.density * np.pi))


print("Successfully imported `Power` class")


class Weather:
        'Class describing ambient weather conditions and is used to predict the drone\'s power requirement'

        # miscellaneous weather   
        weatherlist     = [] # a list of weather objects (e.g., rain, icing, humidity, etc.)

        # ambient air quality
        density = 1.225 #kg/m^3
        gravity = 9.806655 #m/s^2
        altitude = 100 #m
        temperature_sl = 15
        temperature = 15
        pressure_sl = 1
        pressure = 101300 #Pa
        rain = 'False'
        ice = 'False'
        wind = 0
        humidity = 0

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
        LWC = None
        dropsize = None
        WVC = None

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
        speed = None
        heading = 0 #degrees from north maybe?

        def __init__(self,speed,heading):
                self.speed = speed
                self.heading = heading

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

        # methods go here:
        def __init__(self, drone, weather, power):
                # get parameters from `settings.txt`
                print("still working on simulation class constructor")

        # time-invariant methods go here (indicated by suffix *_ti):
        def getEndurance_ti(self, power):
                # check that objects are time-invariant
                if power.timevariant == True:
                        print("")
                        sys.exit("ERR: attempted to run a time-invariant simulation with a time-variant `Power` object")
                else:
                        # calculate endurance
                        return self.capacity * self.voltage_mean / power.power


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