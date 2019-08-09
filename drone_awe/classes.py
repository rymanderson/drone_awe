import datetime
import matplotlib.pyplot as plt
import numpy as np
from .functions import getParams, getXandY, interpolate
from .validationdata import validationdatabase
from .drones import drones
from .drones import conversions
from scipy.optimize import fsolve
import gekko

# insert classes here

class Drone:
    '''
    Class used to store key drone characteristics.
    Class instances are used only to instantiate other classes;
    i.e., Drone class objects are never referenced during a single simulation step
    '''

    name            = None
    params          = None
    conversions     = None
    correctunits    = None
    debug           = False
    log = {
        'SUCCESS': [],
        'WARNING': [],
        'NOTE': [],
        'ERROR': []
    }

    def __init__(self, name, params, conversions, correctunits, debug=False):
        self.debug = debug
        self.name                   = name
        self.params                 = params
        self.params['takeoffweightoriginal'] = self.params['takeoffweight']
        # this input is obtained using functions.getParams() in exe.py
        self.correctunits           = correctunits
        self.conversions            = conversions
        self.__convertUnits()
        # extrapolate useful parameters
        self.params['frontalarea']  = self.params['width'] * self.params['height']
        self.params['toparea']      = self.params['width'] * self.params['length']
        if 'rangemax' in self.params and 'rangemaxspeed' in self.params:
            self.params['endurancemaxrange'] = self.params['rangemax'] / self.params['rangemaxspeed']
        if 'endurancemaxhover' not in self.params and 'endurancemax' in self.params:
            self.params['endurancemaxhover'] = self.params['endurancemax'] * 0.9
            self.__updateLog('WARNING','`endurancemaxhover` not found; assuming 90%% of `endurancemax`.')
            if self.debug:
                print("Drone.__init__:  'endurancemaxhover' not found; assuming 90%% of 'endurancemax'")
        if 'rotordiameter' not in self.params:
            self.__updateLog('ERROR','Rotor diameter not found.')
            raise(Exception("DRONE: ~~~~~ ERROR: rotor diameter not found"))
        if 'rotorarea' not in self.params:
            self.params['rotorarea'] = np.pi*self.params['rotordiameter']**2/4
        if self.params['wingtype'] == 'fixed':
            self.params['spanefficiency'] = 0.8 #estimate from Dr. Ning's book (he lists 0.7-0.9). If we want to we could decrease this further based on fuselage diameter, but maybe that's requiring too much detail
            print("DRONE: ????? Should we be pre-defining a spanwise efficiency?")
        if 'energyfull' not in self.params['battery']:
            self.params['battery']['energyfull'] = self.params['battery']['batteryvoltage']*self.params['battery']['batterycapacity']
        if 'numbatteries' in self.params['battery']:
            if self.params['battery']['numbatteriesconnection'] == 'parallel':
                self.params['battery']['batterycapacity'] *= self.params['numbatteries'] #increase capacity in parallel
            elif self.params['battery']['numbatteriesconnection'] == 'series':
                self.params['battery']['batteryvoltage'] *= self.params['nummbatteries'] #increase voltage in series
            else:
                self.__updateLog('ERROR','Incorrect battery connection parameter applied.')
                raise(Exception("DRONE: ~~~~~ ERROR: incorrect battery connection parameter applied"))
        # print parameters for debugging
        # self.__printParameters()

    def __convertUnits(self):
        if self.debug:
            print("DRONE: ----- drone.correctunits = ",self.correctunits)
            print("DRONE: ----- drone.params contains: ",self.params)
        if not self.correctunits:
            self.__updateLog('NOTE','Converting units.')
            if self.debug:
                print("DRONE: ===== converting units")
            for spec in self.conversions:
                if spec in self.params:
                    self.params[spec] = self.params[spec] * \
                        self.conversions[spec]
                elif spec in self.params['battery']:
                    self.params['battery'][spec] = self.params['battery'][spec] * \
                        self.conversions[spec]
                else:
                    self.__updateLog('NOTE','No units to convert.')
                    if self.debug:
                        print("DRONE: ===== no units to convert")
            self.correctunits = True
            if self.debug:
                print("DRONE: 22222 drone.correctunits = ",self.correctunits)
        else:
            self.__updateLog('WARNING','Problem converting units: skipping.')
            if self.debug:
                print("~~~~~ WARNING: problem converting units: skipping ~~~~~")

    def __printParameters(self):
        print("DRONEIMPORT:     takeoffweight is: ",self.params['takeoffweight'])

    # TYLER: I commented out this class assuming that the power class would handle this. Let me know if you had something else in mind.
    def update(self,weather):
        if "dropletforce" in weather.params:
            self.params['takeoffweight'] = self.params['takeoffweightoriginal'] + \
                (weather.params['dropletforce'] / weather.params['gravitationconstant'])
            self.__updateLog('NOTE','Takeoff weight adjusted due to rain.')
            if self.debug:
                print("Takeoff weight adjusted due to rain.")
    
    def __updateLog(self,entrytype,message):
            self.log[entrytype].append(message)
            # entrytype can be `SUCCESS`, `WARNING`, or `ERROR`

class Battery:
    'Class used to track battery characteristics and performance during simulation'

    params = {
                'batterycapacity': None,           # Ampere-hours
                'stateofcharge': 100,             # state of charge (in percent nominal capacity)
                'startstateofcharge': 100,        # state of charge at simulation start
                'stateofhealth': 100,             # state of health (actual capacity divided by ideal capacity)
                'batterytype': None,     # possible values include LiPo, Li-ion, NiCd, NiMH, SLA
                'voltage': None,         # Volts; this is the instantaneous voltage
                'voltagemean': None,     # Volts; average voltage used for time-invariant simulations
                'voltagecharged': None,  # Volts
                'voltagedead': None,     # Volts
                'current': None,         # Amperes; this is the instantaneous current
                'batterytechnology': 'current' # current, near-future, or far-future
            }

    debug = False
    log = {
        'SUCCESS': [],
        'WARNING': [],
        'NOTE': [],
        'ERROR': []
    }

    # constructor
    # default value for soh is based on the assumption that batteries are retired at a soh of 80%
    def __init__(self, drone, stateofhealth=100.0, startstateofcharge=100.0, batterytechnology='current', debug=False):
        # update parameters
        self.debug = debug
        self.update(drone)

        # estimate list lengths for prior memory allocation

    def __defineCapacity(self):
        if self.debug:
            print("defineCapacity:  batterytechnology is: ",self.params['batterytechnology'])
        if self.params['batterytechnology'] == 'current':
            self.params['batterycapacity'] = self.params['batterycapacity']
        elif self.params['batterytechnology'] == 'near-future':
            self.__updateLog('NOTE','Assuming a LiPo battery capacity increase of 3.5% per year for 5 years.')
            if self.debug:
                print("Assuming a LiPo battery capacity increase of 3.5% per year for 5 years.")
            self.params['batterycapacity'] = self.params['batterycapacity'] * (1.035**5) # capacity increases by 3-4% each year, this assumes 3.5% for 5 years. It may be an option to let the user specify a timeline, but past 5 years I don't know if that growth in capacity is sustainable.
        elif self.params['batterytechnology'] == 'far-future': #Li-air batteries
            self.__updateLog('NOTE','Assuming Lithium-air batteries with a capacity of ____ mAh.')
            if self.debug:
                print("Assuming Lithium-air batteries with a capacity of ____ mAh.")
            self.params['batterycapacity'] = self.params['batterycapacity'] * 10 # estimate for now - ten times the capacity
        else:
            self.__updateLog('ERROR','Incompatible battery technology input')
            raise(Exception("ERROR: Incompatible battery technology input."))

    def update(self,drone):
        # import required parameters from drone object
        self.params['batterycapacity'] = drone.params['battery']['batterycapacity']
        self.params['voltagemean'] = drone.params['battery']['batteryvoltage']
        # import optional parameters from drone object
        for param in ['batterytechnology','batterytype','batterycells']:
            if param in drone.params['battery']:
                self.params[param] = drone.params['battery'][param]

        #update capacity based on future technology if needed
        self.__defineCapacity()
        
    # time-variant methods go here:
    def updateLoad(self, power):
        self.current = power.params['power'] / self.params['voltage']

    def discharge(self, power, timestep):
        self.soc = self.soc - power.params['power'] * timestep
        self.current = power.params['power'] / self.params['voltage']

    def __updateLog(self,entrytype,message):
        self.log[entrytype].append(message)
        # entrytype can be `SUCCESS`, `WARNING`, or `ERROR`


class Power:
    'Class used to predict the drone\'s power requirement'

    # class variables go here:
    params = {
        'efficiencypropulsive': None,
        'power': None,
        'model': None,
        'dragcoefficient': 1.0,
        'alpha': 0.0,
        'model': None,
        'powerModel': None
    }

    log = {
        'SUCCESS': [],
        'NOTE': [],
        'WARNING': [],
        'ERROR' : []
    }

    # methods go here:
    def __init__(self, drone, weather, mission, debug=False):
        self.debug = debug
        self.__initializeParameters(drone,weather,mission)         
        self.__setDragCoefficient()                # hopefully constant across different rotary drones
        # self.__setupModel(mission.params['missionspeed'],weather.params['airdensity'],drone.params['rotorarea'],drone.params['rotorquantity'],drone.params['frontalarea'],drone.params['toparea'])

        if drone.params['wingtype'] == 'rotary':
            # self.__getPowerAbdilla(drone, weather, mission)
            self.params['powerModel'] = self.__setPowerMomentum
        elif drone.params['wingtype'] == 'fixed':
            self.params['powerModel'] = self.__setPowerTraub
        else:
            self.__updateLog('ERROR','No power model defined for ' + drone.params['wingtype'])
            if self.debug:
                print("POWER: ~~~~~ ERROR: no power model defined for drone.params['wingtype']")

    def __initializeParameters(self,drone,weather,mission):
        if self.debug:
            print("POWER: ===== mission.params = ",mission.params)
        self.__updateWeight(drone.params['takeoffweight'],weather.params['gravitationconstant'],weather.params['weightadjustment'],mission.params['payload'])
        self.__updateVelocityInducedHover(self.params['effectiveweight'],drone.params['rotorquantity'],drone.params['rotorarea'],weather.params['airdensity'])
        self.params['velocityinduced']          = self.params['velocityinducedhover']
        self.params['thrust']                   = self.params['effectiveweight']
        self.params['alpha']                    = 0.0
        self.params['alpha_gekko']              = 0.0
        self.params['drag']                     = 0.0
        self.params['dragarea']                 = drone.params['frontalarea']
        self.__setDragCoefficient()


        # frontalarea         = drone.params['frontalarea']
        # toparea             = drone.params['toparea']

    def __printParameters(self,drone,weather,mission):
        velocityinfinity    = mission.params['missionspeed']
        effectiveweight     = self.params['effectiveweight']
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
        thrust3 = np.sqrt((effectiveweight/rotorquantity)**2 + (1/2*airdensity*velocityinfinity**2 * \
                dragcoefficient * (toparea*np.sin(alpha) + frontalarea*np.cos(alpha)))**2)
        velocityinduced2 = effectiveweight/rotorquantity/2.0/airdensity/rotorarea/np.sqrt((velocityinfinity*np.cos(alpha))**2 + \
                           (velocityinfinity*np.sin(alpha) + velocityinduced)**2)

        # print("----- MISSION SPEED = ",mission.params['missionspeed']," -----")
        # print("")
        # print("power:           alpha is                ",self.params['alpha']*180/np.pi)
        # print("power:           thrust is               ",self.params['thrust'])
        # print("power:           thrust2 is              ",thrust2)
        # print("power:           thrust3 is              ",thrust3)
        # print("power:           drag is                 ",self.params['drag'])
        # print("power:           velocityinduced is      ",self.params['velocityinduced'])
        # print("power:           velocityinduced2 is     ",velocityinduced2)
        # print("power:           velocityinducedhover is ",self.params['velocityinducedhover'])
        # print("power:           etapropulsive is        ",self.params['efficiencypropulsive'])
        # print("")
        # print("drone:           effectiveweight is          ",drone.params['effectiveweight'])
        # print("drone:           takeoffweight is        ",drone.params['takeoffweight'])
        # print("drone:           payload is              ",drone.params['payload'])
        # print("drone:           rotorquantity is        ",drone.params['rotorquantity'])
        # print("")
        # print("weather:         airdensity is           ",weather.params['airdensity'])
        # print("weather:         gravitationconstant is  ",weather.params['gravitationconstant'])
        # print("weather:         temperature is          ",weather.weatherlist[0].params['temperature'])
        # print("weather:         relativehumidity is     ",weather.weatherlist[1].params['relativehumidity'])
        # print("")  
        # print("")

    def update(self, drone, weather, mission):
        self.__updateWeight(drone.params['takeoffweight'],weather.params['gravitationconstant'],weather.params['weightadjustment'],mission.params['payload'])
        self.__updateVelocityInducedHover(self.params['effectiveweight'],drone.params['rotorquantity'],drone.params['rotorarea'],weather.params['airdensity'])
        self.params['powerModel'](drone, weather, mission)

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

    def __updateWeight(self,takeoffweight,g,weightadjustment,payload):
        self.params['effectiveweight']      = takeoffweight*g + \
                                                weightadjustment + \
                                                payload*g
    
    def __updateVelocityInducedHover(self,effectiveweight,rotorquantity,rotorarea,airdensity):
        self.params['velocityinducedhover'] = np.sqrt(effectiveweight / \
                                                rotorquantity / 2 / \
                                                rotorarea / airdensity)

    # def __getPowerAbdilla(self, drone, weather, mission):
    #     self.__updateEfficiencyPropulsive(drone, weather, mission)
    #     self.params['power'] = self.params['thrust']**1.5 / \
    #                            self.params['efficiencypropulsive'] / \
    #                             np.sqrt(2 * drone.params['rotorquantity'] * \
    #                             weather.params['airdensity'] * drone.params['rotorarea']) + \
    #                             0.0  # Camera power consumption estimate

        # self.__printParameters(drone,weather,mission)

    def __setPowerMomentum(self, drone, weather, mission):
        self.__setupModel(mission.params['missionspeed'],weather.params['airdensity'],drone.params['rotorarea'],drone.params['rotorquantity'],drone.params['frontalarea'],drone.params['toparea'])
        # self.__solveModel()
        # self.__getDrag(drone,weather,mission)
        # self.__getThrust(drone,weather,mission)
        self.__setBladeProfilePower()
        self.params['efficiencypropulsive'] = self.__getEfficiencyPropulsive(drone.params['endurancemaxhover'],drone.params['takeoffweight'],9.81,drone.params['rotorquantity'],drone.params['rotorarea'],drone.params['battery']['energyfull'],1.225)
        self.params['power']                = (self.params['thrust']*self.params['velocityinduced'] + \
                                                self.params['drag']*mission.params['missionspeed'] + \
                                                self.params['bladeprofilepower'])/self.params['efficiencypropulsive']
        if self.debug:
            self.__printParameters(drone,weather,mission)

    def __setPowerTraub(self, drone, weather, mission): #fixed-wing power model
        density = weather.params['airdensity']
        cruisespeed = mission.params['missionspeed']
        if 'wingarea' in drone.params:
            wingarea = drone.params['wingarea']
        elif 'wingspan' in drone.params and 'chord' in drone.params:
            wingarea = drone.params['wingspan'] * drone.params['chord']
        else:
            self.__updateLog('ERROR','Wing area needed to calculate power')
            raise(Exception("~~~~~ ERROR: wing area needed to calculate power ~~~~~"))
        if 'wingspan' not in drone.params:
            self.__updateLog('ERROR','Wingspan needed to calculate power')
            raise(Exception("~~~~~ ERROR: wingspan needed to calculate power ~~~~~"))
        else:
            span = drone.params['wingspan']
        self.__setDragCoefficient()
        dragcoefficient  = self.params['dragcoefficient']
        weight           = drone.params['takeoffweight'] * weather.params['gravitationconstant']
        spanefficiency   = drone.params['spanefficiency']
        k                = 1 / (np.pi*span**2 / wingarea * spanefficiency)

        self.params['power'] = 0.5 * density * cruisespeed**3 * wingarea * dragcoefficient \
                                + 2*weight**2*k \
                                / (density*cruisespeed*wingarea)

    def getPowerFixedWing(self, drone, weather, mission): #from Dr. Nin's textbook
        
        LD = 10 #placeholder, this could be a user input

        if 'LDadjustment' in weather.params: #this probably can be its own method
            L = drone.params['effectivetakeoffweight'] * weather.params['gravitationconstant']
            D = L / LD
            CLfactor = weather.params['LDadjustment'][0]
            CDfactor = weather.params['LDadjustment'][1]
            L *= CLfactor
            D += CDfactor
            LD = L/D

        propulsiveefficiency = .4 #placeholder

        self.params['power'] = weather.params['gravitationconstant'] * drone.params['effectivetakeoffweight'] * mission.params['missionspeed'] / (propulsiveefficiency * LD)


    # def __updateEfficiencyPropulsive(self, drone, weather, mission):
    #     # default value:
    #     if 'endurancemax' not in drone.params or 'endurancemaxspeed' not in drone.params:
    #         self.params['efficiencypropulsive'] = 0.5
    #     else:
    #         # get efficiency at max endurance conditions
    #         etamaxendurance = self.__getEfficiencyPropulsive(
    #             drone, weather, drone.params['endurancemax'])
    #         vmaxendurance = drone.params['endurancemaxspeed']
    #         # get efficiency at max range conditions
    #         etamaxrange = self.__getEfficiencyPropulsive(
    #             drone, weather, drone.params['rangemax']/drone.params['rangemaxspeed'])
    #         vmaxrange = drone.params['rangemaxspeed']
    #         # interpolate for the current velocity
    #         self.params['efficiencypropulsive'] = etamaxendurance - (vmaxendurance - mission.params['missionspeed']) / \
    #             (vmaxendurance - vmaxrange) * (etamaxendurance - etamaxrange)

    def __getEfficiencyPropulsive(self, testmass, g, rotorquantity, rotorarea, airdensity, batteryenergy, endurance):
        # Analyzing a single propeller
        thrust          = testmass * g / rotorquantity
        batteryenergy   = batteryenergy / \
                          rotorquantity
        poweractual     = batteryenergy/endurance
        powerideal      = thrust**1.5 / np.sqrt(2*rotorarea*airdensity)

        # print('gEP: poweractual is      ',poweractual)
        # print('gEP: powerideal is      ',powerideal)

        efficiencypropulsive = powerideal/poweractual

        return efficiencypropulsive

        # prediction based on momentum theory for hover case
        # slides from https://fenix.tecnico.ulisboa.pt/downloadFile/282093452028191/3-Momentum%20Theory%20in%20hover.pdf
        # prediction based on momentum theory for forward flight
        # slides from https://fenix.tecnico.ulisboa.pt/downloadFile/845043405440064/6-Momentum%20Theory%20in%20forward%20flight.pdf
        # see also: http://www.aerospaceweb.org/design/helicopter/momentum.shtml

    def __setupModel(self,velocityinfinity,airdensity,rotorarea,rotorquantity,frontalarea,toparea): #
        if self.debug:
            self.__print5('__solveModel BEFORE ')

            print('POWER: ----- frontalarea = ',frontalarea)
            print('POWER: ----- toparea       ',toparea)
            print('POWER: ----- velocityinf   ',velocityinfinity)
            print('POWER: ----- airdensity    ',airdensity)
            print('POWER: ----- rotorarea     ',rotorarea)
            print('POWER: ----- rotorquantity ',rotorquantity)
            print('')
            print('')

        effectiveweight   = self.params['effectiveweight']
        dragcoefficient   = self.params['dragcoefficient']

        m                 = gekko.GEKKO(remote=False)             # create GEKKO modelz
        m.alpha           = m.Var(value=self.params['alpha_gekko'])      # define new variable, initial value=0
        m.velocityinduced = m.Var(value=self.params['velocityinduced'])      # define new variable, initial value
        m.thrust          = m.Var(value=self.params['thrust'])
        # dragcoefficient = m.Var(value=1.0)
        m.drag            = m.Var(value=self.params['drag'])
        m.dragarea        = m.Var(value=self.params['dragarea'])
        # print values for debugging
        # print("GEKKO:       airdensity*velocityinfinity**2*pi*toparea = ",airdensity*velocityinfinity**2*pi*0.006)

        qty             = (effectiveweight/rotorquantity + 0.00000000/rotorquantity)**2

        m.Equations([ \
            m.thrust          == m.sqrt(effectiveweight**2 + m.drag**2), \
            m.tan(m.alpha)    == m.drag / effectiveweight, \
            m.drag            == 0.5*airdensity*velocityinfinity**2 * dragcoefficient * m.dragarea, \
            m.velocityinduced == m.thrust / (2*rotorarea*rotorquantity*airdensity*m.sqrt((velocityinfinity*m.cos(m.alpha))**2 + (velocityinfinity*m.sin(m.alpha) + m.velocityinduced)**2)), \
            m.dragarea        == (frontalarea * m.cos(m.alpha) + toparea * m.sin(m.alpha)), \
            m.dragarea        > 0.0, \
            m.alpha           < 1.5708, \
            m.alpha           > 0.0
            ])
        
        self.params['model'] = m
    
    # def __warmModel(self,drone,weather,mission):
    #     velocityinfinity    = mission.params['missionspeed']
    #     effectiveweight         = drone.params['effectiveweight']
    #     airdensity          = weather.params['airdensity']
    #     rotorarea           = drone.params['rotorarea']
    #     rotorquantity       = drone.params['rotorquantity']
    #     dragcoefficient     = self.params['dragcoefficient']
    #     frontalarea         = drone.params['frontalarea']
    #     toparea             = drone.params['toparea']
    #     dynamicpressure     = 0.5*airdensity*velocityinfinity**2

    #     self.params['model'].alpha.value           = self.params['alpha']
    #     self.params['model'].velocityinduced.value = self.params['velocityinduced']
    #     self.params['model'].thrust.value          = self.params['thrust']
    #     # dragcoefficient = m.Var(value=1.0)
    #     self.params['model'].drag.value            = self.params['drag']
    #     self.params['model'].area.value            = self.params['area']

    #     self.__print5('__warmModel')

    #     print(" *** +/= : EQUATIONs = ",self.params['model'].Equations, "before")

    #     self.params['model'].Equations = [ \
    #         self.params['model'].thrust          == self.params['model'].sqrt(effectiveweight**2 + self.params['model'].drag**2), \
    #         self.params['model'].tan(self.params['model'].alpha)    == self.params['model'].drag / effectiveweight, \
    #         self.params['model'].drag            == dynamicpressure * dragcoefficient * self.params['model'].area, \
    #         self.params['model'].velocityinduced == self.params['model'].thrust / (2*rotorarea*rotorquantity*airdensity*self.params['model'].sqrt(( \
    #                                                 velocityinfinity*self.params['model'].cos(self.params['model'].alpha))**2 + (velocityinfinity* \
    #                                                 self.params['model'].sin(self.params['model'].alpha) + self.params['model'].velocityinduced)**2)), \
    #         self.params['model'].area            == (frontalarea * self.params['model'].cos(self.params['model'].alpha) + toparea * self.params['model'].sin(self.params['model'].alpha)) \
    #         ]

    #     print(" *** >/> : EQUATIONs = ",self.params['model'].Equations, "after")

    # def __solveModel(self):
        m.options.SOLVER=3
        m.solve(disp=False)     # solve

        self.params['alpha_gekko']      = m.alpha.value[0]
        alpha                           = m.alpha.value[0]
        # ensure alpha is between -pi and pi
        if alpha > np.pi or alpha < -np.pi:
            print("NOTE: alpha was evaluated as ",alpha,"; adjusting")
            while alpha > np.pi or alpha < -np.pi:
                if alpha > np.pi:
                    alpha = alpha - 2*np.pi
                elif alpha < -np.pi:
                    alpha = alpha + 2*np.pi
        self.params['alpha']            = alpha        
        self.params['velocityinduced']  = m.velocityinduced.value[0]
        self.params['thrust']           = m.thrust.value[0]
        self.params['drag']             = m.drag.value[0]
        self.params['dragarea']         = m.dragarea.value[0]

        if self.debug:
            self.__print5('__solveModel AFTER')
            print('POWER: ----- frontalarea = ',frontalarea)
            print('POWER: ----- toparea       ',toparea)
            print('POWER: ----- velocityinf   ',velocityinfinity)
            print('POWER: ----- airdensity    ',airdensity)
            print('POWER: ----- rotorarea     ',rotorarea)
            print('POWER: ----- rotorquantity ',rotorquantity)
            print('POWER: ----- power         ',self.params['power'])
            print('')
            print('')

    def __print5(self,functionname):
        print('')
        print('===== ',functionname,' =====')
        print('     > thrust    = ',self.params['thrust'])
        print('     > alpha     = ',self.params['alpha'])
        print('     > velocityi = ',self.params['velocityinduced'])
        print('     > drag      = ',self.params['drag'])
        print('     > area      = ',self.params['dragarea'])
        print('     > weight    = ',self.params['effectiveweight'])
        print('     > cd        = ',self.params['dragcoefficient'])
        print('     > ')
        print('')
        pass

    # def __getDrag(self,drone,weather,mission):
    #     drag = 0.5 * self.params['dragcoefficient'] * weather.params['airdensity'] * mission.params['missionspeed']**2 * \
    #            (drone.params['frontalarea'] * np.cos(self.params['alpha']*np.pi/180.0) + \
    #             drone.params['toparea'] * np.sin(self.params['alpha']*np.pi/180.0))
    #     self.params['drag'] = drag

    # def __getThrust(self,drone,weather,mission):
    #     effectiveweight = drone.params['effectiveweight']
    #     drag        = self.params['drag']
    #     thrust      = np.sqrt(effectiveweight**2 + drag**2)

    #     self.params['thrust'] = thrust

    def __setDragCoefficient(self):
        self.params['dragcoefficient'] = 2.0
    
    def __setBladeProfilePower(self):
        self.params['bladeprofilepower'] = 50.0

    def __updateLog(self,entrytype,message):
        self.log[entrytype].append(message)
        # entrytype can be `SUCCESS`, `WARNING`, or `ERROR`, `NOTE`

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
        'altitude': None,
        'weightadjustment': 0.0
    }

    log     = {
        'SUCCESS': [],
        'WARNING': [],
        'ERROR': [],
        'NOTE': []
    }
    debug   = False

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
    def __init__(self, altitude, weatherlist, debug=False):  # keeping it simple to begin with
        self.weatherlist    = weatherlist
        self.altitude       = altitude
        self.debug          = debug
        # self.temperature = temperaturesealevel
        # self.temperature = self.temperaturesealevel - 71.5 + 2*np.log(1 + np.exp(35.75 - 3.25*self.altitude) + np.exp(-3 + 0.0003 * self.altitude**3))
        # self.pressure = self.pressure_sl * np.exp(-0.118 * self.altitude - (0.0015*self.altitude**2) / (1 - 0.018*self.altitude + 0.0011 * self.altitude**2))

        # set self.params        
        for weathertype in weatherlist:
            for param in weathertype.params:
                self.params[param] = weathertype.params[param]

    def update(self,drone):
        # update independent parameters
        # print("Weatherlist:", self.weatherlist)
        
        # update weather effects:
        densityfactor = 1.0
        densityfactor *= self.__updateDensityTemperature()
        densityfactor *= self.__updateDensityHumidity()
        self.params['weightadjustment'] = self.__updateRain(drone.params['toparea'])
        if drone.params['wingtype'] == 'fixed':
            self.params['LDadjustment'] = self.updateLD() #list with [CLfactor, CD factor]
        self.params['airdensity'] = self.params['airdensitysealevel'] * densityfactor #for the case that we update weather multiple times, we need to not compound weather effects
        self.__warning()

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


    # Rain Methods:

    def __updateRain(self, dronearea):
        
        diameter = self.params['dropsize']
        liquidwatercontent = self.params['liquidwatercontent'] / 1000 #convert to kg/m3
        if self.debug:
            print("Rain.update is in process.")
            print("diameter =",diameter)
            print("liquidwatercontent =",liquidwatercontent)

        if diameter == 0 or liquidwatercontent == 0:
            dropletforce = 0.0
        else:
            
            # single droplet
            diameter = self.params['dropsize']
            radius = diameter/2.0
            dropletvolume = 4.0/3.0 * np.pi * radius**3 #sphere estimate
            airdensity = self.params['airdensity'] #kg/m^3
            waterdensity = 1000 #kg/m^3
            dropletmass = dropletvolume * waterdensity #kg
            g = self.params['gravitationconstant']
            Cd = 0.5 #sphere estimate
            area = np.pi * radius**2
            terminalvelocity = np.sqrt(2*dropletmass*g / (Cd*airdensity*area)) 
            # terminalvelocity2 = 9.58*(1-np.exp(-(diameter/1.77)**1.147)) #(from Cao et al.)
            # print("terminalvelocity1",terminalvelocity)
            # print("terminalvelocity2 =",terminalvelocity2)
            
            #determine raindrop incidence using rate - Method 1
            # rate = self.params['rainfallrate']
            # rate *= 0.001 #convert from L/m3/h to m3/m2/h
            # rate /= 3600 #convert from m3/m3/h to m3/m2/s
            # dronearea = drone.params['area']
            # rate *= dronearea #m3 per drone area per second
            # numdrops = rate / dropletvolume #number of drops per drone area per second

            #Method 2 - using LWC
            numdrops = liquidwatercontent / dropletmass * dronearea * terminalvelocity #numsrops hitting drone per second

            surfacetension = self.__getSurfaceTension()
            webernumber = self.__getWeberNumber(terminalvelocity, diameter, surfacetension)

            upperthreshold = 18.0**2 * diameter * (waterdensity/surfacetension)**0.5 * terminalvelocity**0.25 * numdrops**0.75
            # print("rain upper threshold weber number is:",upperthreshold)

            if webernumber < 5:
                #raindrop sticks
                C = 0.0
            elif webernumber < 10:
                #raindrop bounces back
                C = 1.0
            elif webernumber < upperthreshold:
                #raindrop spreads
                C = 0.0
            else:
                #raindrop splashes
                C = 2/np.pi #estimate, (integral of sin(x) from 0 to pi) / pi
            #calculate change in velocity
            delta_velocity = terminalvelocity * (C+1)
            
            dropletforce = numdrops * (dropletmass*delta_velocity)

            if self.debug:
                print("")
                print("Rain Parameters:")
                print("dropletmass is:", dropletmass)
                print("dropletvolume is:", dropletvolume)
                print("velocity is:", terminalvelocity)
                print("webernumber is:", webernumber)
                print("numdrops is:", numdrops)
                print("dropletforce is:", dropletforce)
                print("")

        return dropletforce

    def __getWeberNumber(self, velocity, diameter, surfacetension):
        density = 1000 #kg/m^3
        webernumber = density * velocity**2 * diameter / surfacetension

        return webernumber

    def __getSurfaceTension(self):
        Tlist = [0, 5, 10, 20, 30, 40, 50] #deg C
        sigmalist = [7.56, 7.49, 7.42, 7.28, 7.12, 6.96, 6.79] #N/m
        temperature = self.params['temperature']

        if temperature < Tlist[0]:
            surfacetension = sigmalist[0]
            self.__updateLog('WARNING','Temperature is too low to predict surface tension of raindrops. Assuming a surface tension at 0 degrees C.')
            if self.debug:
                print("temperature is too low to predict surface tension of raindrops. Assuming a surface tension at 0 degrees C.")
        elif temperature > Tlist[-1]:
            surfacetension = sigmalist[-1]
            self.__updateLog('WARNING','temperature is too high to predict surface tension of raindrops. Assuming a surface tension at 50 degrees C.')
            if self.debug:
                print("temperature is too high to predict surface tension of raindrops. Assuming a surface tension at 50 degrees C.")
        else:
            counter = 0
            for temp in Tlist:
                if temperature == temp:
                    id1 = counter
                    id2 = counter
                    break
                elif temperature < temp:
                    id1 = counter - 1
                    id2 = counter
                    break
                else:
                    counter += 1
            #interpolate
            x1 = Tlist[id1]
            x2 = Tlist[id2]
            y1 = sigmalist[id1]
            y2 = sigmalist[id2]
            x = temperature
            surfacetension = interpolate(x1,x2,y1,y2,x)

        return surfacetension

    def updateLD(self):
            CLfactor = 0.94 #94% of original value
            CDfactor = 0.01 #added on to original value 
            #if CD needs to be percentage, increasing by 20% is about average, although less acurate than an addition. 
            return [CLfactor, CDfactor]

    # Icing Methods:

    def __warning(self):

        if self.params['temperature'] < 5:
            if self.params['relativehumidity'] <= 0:
                self.__updateLog('WARNING','Cold temperatures may result in dangerous icing conditions. Avoid routes near clouds and/or rain.')
                if self.debug:
                    print("WARNING: Cold temperatures may result in dangerous icing conditions. Avoid routes near clouds and/or rain.")
            elif self.params['relativehumidity'] < 30:
                self.__updateLog('WARNING','Cold temperatures and moderate humidity could result in dangerous icing conditions. Be cautious when flying.')
                if self.debug:
                    print("WARNING: Cold temperatures and moderate humidity could result in dangerous icing conditions. Be cautious when flying.")
            else:
                self.__updateLog('WARNING','Cold temperatures and high humidity will likely result in dangerous icing conditions. Take extreme caution if flying.')
                if self.debug:
                    print("WARNING: Cold temperatures and high humidity will likely result in dangerous icing conditions. Take extreme caution if flying.")

    # Temperature Methods:

    def __updateDensityTemperature(self):
        newtemperature = self.params['temperature'] +273.15 #convert to kelvin
        oldtemperature = self.params['temperaturesealevel'] #convert to kelvin
        # print("old temperature is",oldtemperature)
        # print("new temperature is",newtemperature)
        # print("Self.params['temperature'] is",self.params['temperature'])
        olddensity = self.params['airdensity']
        # print("old density is",olddensity)
        newdensity = olddensity * ((oldtemperature-newtemperature)/oldtemperature + 1) #inverse relationship
        temperatureeffect = newdensity / olddensity

        return temperatureeffect
    
    # Humidity Methods:

    def __updateDensityHumidity(self): # 2D linear interpolation based on Yue (2017)
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

        relativehumidity = self.params['relativehumidity']
        # temperature = weather.params['temperaturesealevel']
        # temperature = 288.15 - 273.15 #convert to celcius - need to update this if temperature is altered
        temperature = self.params['temperature'] #already to celcius
        # print("Temperature is",temperature)
        # determine 4 closest points before doing a 2D interpolation

        if temperature < temperaturelist[0]:
            temperature = temperaturelist[0]
            self.__updateLog('WARNING','Temperature is below interpolation bounds. Humidity effects are assumed to be as if the temperature is'+str(temperaturelist[0])+'degrees C')
            if self.debug:
                print("WARNING: Temperature is below interpolation bounds. Humidity effects are assumed to be as if the temperature is",temperaturelist[0],"degrees C")
        elif temperature > temperaturelist[-1]:
            temperature = temperaturelist[-1]
            self.__updateLog('WARNING','Temperature is below interpolation bounds. Humidity effects are assumed to be as if the temperature is'+str(temperaturelist[0])+'degrees C')
            if self.debug:
                print("WARNING: Temperature is below interpolation bounds. Humidity effects are assumed to be as if the temperature is",temperaturelist[0],"degrees C")
        
        if relativehumidity < 0.0 or relativehumidity > 100.0:
            self.__updateLog('ERROR','Relative humidity is outside of available bounds')
            raise(Exception("~~~~~ ERROR: relative humidity is outside of available bounds ~~~~~"))
        elif relativehumidity > relativehumiditylist[-1]:
            relativehumidity = relativehumiditylist[-1]
            self.__updateLog('WARNING','Humidity is above interpolation bounds. Humidity effects are assumed to be as if the temperature is'+str(relativehumiditylist[-1])+"%")
            print("WARNING: Humidity is above interpolation bounds. Humidity effects are assumed to be as if the temperature is",relativehumiditylist[-1],"%")

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
        tempid1 = interpolate(x1,x2,y1,y2,x)
        
        x1 = temperaturelist[idT1]
        x2 = temperaturelist[idT2]
        y1 = humidityarray[idH2,idT1]
        y2 = humidityarray[idH2,idT2]
        x = temperature
        tempid2 = interpolate(x1,x2,y1,y2,x)

        #now interpolate based on humidity
        x1 = relativehumiditylist[idH1]
        x2 = relativehumiditylist[idH2]
        y1 = tempid1
        y2 = tempid2
        x = relativehumidity
        humidityeffect = interpolate(x1,x2,y1,y2,x)

        return humidityeffect

    def __updateLog(self,entrytype,message):
        self.log[entrytype].append(message)
        # entrytype can be `SUCCESS`, `WARNING`, or `ERROR`


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
    * LWC (liquid water content) [g/m3]
    * dropsize (refers to drpolet diameter) [m]
    * rainfallrate [mm/h] - 1 mm of rain is equivalent to 1 L/m^2

    Rain is expected to affect:
    * lift
    * drag
    * weight
    * impart downward momentum
    '''

    def __init__(self, params):
        if 'liquidwatercontent' not in params:
            if 'rainfallrate' not in params:
                raise(Exception('ERROR: either liquidwatercontent (LWC) or rainfall rate must be specified.'))
            R = params['rainfallrate']
            if R < 0.5: #drizzle
                params['liquidwatercontent'] = 30000*np.pi * 10**-3 / (5.7**4) * R**0.84 #see Cao, pp. 89
            elif R < 4.0: #moderate
                params['liquidwatercontent'] = 7000*np.pi * 10**-3 / (4.1**4) * R**0.84
            elif R > 4.0: #heavy
                params['liquidwatercontent'] = 1400*np.pi * 10**-3 / (3**4) * R**0.84
        elif 'rainfallrate' not in params:
            liquidwatercontent = params['liquidwatercontent']
            if liquidwatercontent < 0.05: #drizzle
                params['rainfallrate'] = 1/((30000*np.pi*10**-3 / 5.7**4*liquidwatercontent)**(1/0.84))
            elif liquidwatercontent < 0.25: #moderate
                params['rainfallrate'] = 1/((7000*np.pi*10**-3 / 4.1**4*liquidwatercontent)**(1/0.84))
            elif liquidwatercontent > 0.25: #heavy
                params['rainfallrate'] = 1/((1400*np.pi*10**-3 / 3.0**4*liquidwatercontent)**(1/0.84))
        paramnames = ['dropsize', 'rainfallrate', 'liquidwatercontent']

        WeatherType.__init__(self, params, paramnames)

    # def update(self, weather, drone):
        
    #     diameter = weather.weatherlist[2].params['dropsize']
    #     liquidwatercontent = weather.weatherlist[2].params['dropsize'] / 1000 #convert to kg/m3
    #     # print("Rain.update is in process.")
    #     print("diameter =",diameter)
    #     print("liquidwatercontent =",liquidwatercontent)

    #     if diameter == 0 or liquidwatercontent == 0:
    #         dropletforce = 0.0
    #     else:
            
    #         # single droplet
    #         diameter = weather.weatherlist[2].params['dropsize']
    #         radius = diameter/2.0
    #         dropletvolume = 4.0/3.0 * np.pi * radius**3 #sphere estimate
    #         airdensity = weather.params['airdensity'] #kg/m^3
    #         waterdensity = 1000 #kg/m^3
    #         dropletmass = dropletvolume * waterdensity #kg
    #         g = weather.params['gravitationconstant']
    #         Cd = 0.5 #sphere estimate
    #         area = np.pi * radius**2
    #         terminalvelocity = np.sqrt(2*dropletmass*g / (Cd*airdensity*area)) 
    #         # terminalvelocity2 = 9.58*(1-np.exp(-(diameter/1.77)**1.147)) #(from Cao et al.)
    #         # print("terminalvelocity1",terminalvelocity)
    #         # print("terminalvelocity2 =",terminalvelocity2)
            
    #         #determine raindrop incidence using rate - Method 1
    #         # rate = self.params['rainfallrate']
    #         # rate *= 0.001 #convert from L/m3/h to m3/m2/h
    #         # rate /= 3600 #convert from m3/m3/h to m3/m2/s
    #         # dronearea = drone.params['area']
    #         # rate *= dronearea #m3 per drone area per second
    #         # numdrops = rate / dropletvolume #number of drops per drone area per second

    #         #Method 2 - using LWC
    #         dronearea = drone.params['toparea']
    #         numdrops = liquidwatercontent / dropletmass * dronearea * terminalvelocity #numsrops hitting drone per second

    #         surfacetension = self.getSurfaceTension(weather)
    #         webernumber = self.getWeberNumber(weather, terminalvelocity, diameter, surfacetension)

    #         upperthreshold = 18.0**2 * diameter * (waterdensity/surfacetension)**0.5 * terminalvelocity**0.25 * numdrops**0.75
    #         # print("rain upper threshold weber number is:",upperthreshold)

    #         if webernumber < 5:
    #             #raindrop sticks
    #             C = 0.0
    #         elif webernumber < 10:
    #             #raindrop bounces back
    #             C = 1.0
    #         elif webernumber < upperthreshold:
    #             #raindrop spreads
    #             C = 0.0
    #         else:
    #             #raindrop splashes
    #             C = 2/np.pi #estimate, (integral of sin(x) from 0 to pi) / pi
    #         #calculate change in velocity
    #         delta_velocity = terminalvelocity * (C+1)
            
    #         dropletforce = numdrops * (dropletmass*delta_velocity)

    #         print("")
    #         print("Rain Parameters:")
    #         print("dropletmass is:", dropletmass)
    #         print("dropletvolume is:", dropletvolume)
    #         print("velocity is:", terminalvelocity)
    #         print("webernumber is:", webernumber)
    #         print("numdrops is:", numdrops)
    #         print("dropletforce is:", dropletforce)
    #         print("")

    #     return dropletforce

    #     # if drone.params['wingtype'] == 'fixed':
    #     #     #TODO loss of lift, increase in drag

    # def getWeberNumber(self, weather, velocity, diameter, surfacetension):
    #     density = 1000 #kg/m^3
    #     webernumber = density * velocity**2 * diameter / surfacetension

    #     return webernumber

    # def getSurfaceTension(self, weather):
    #     Tlist = [0, 5, 10, 20, 30, 40, 50] #deg C
    #     sigmalist = [7.56, 7.49, 7.42, 7.28, 7.12, 6.96, 6.79] #N/m
    #     temperature = weather.weatherlist[0].params['temperature']

    #     if temperature < Tlist[0]:
    #         surfacetension = sigmalist[0]
    #         print("temperature is too low to predict surface tension of raindrops. Assuming a surface tension at 0 degrees C.")
    #     elif temperature > Tlist[-1]:
    #         surfacetension = sigmalist[-1]
    #         print("temperature is too high to predict surface tension of raindrops. Assuming a surface tension at 50 degrees C.")
    #     else:
    #         counter = 0
    #         for temp in Tlist:
    #             if temperature == temp:
    #                 id1 = counter
    #                 id2 = counter
    #                 break
    #             elif temperature < temp:
    #                 id1 = counter - 1
    #                 id2 = counter
    #                 break
    #             else:
    #                 counter += 1
    #         #interpolate
    #         x1 = Tlist[id1]
    #         x2 = Tlist[id2]
    #         y1 = sigmalist[id1]
    #         y2 = sigmalist[id2]
    #         x = temperature
    #         surfacetension = interpolate(x1,x2,y1,y2,x)

    #     return surfacetension

    # def updateLD(self):
    #     CLfactor = 0.94 #94% of original value
    #     CDfactor = 0.01 #added on to original value 
    #     #if CD needs to be percentage, increasing by 20% is about average, although less acurate than an addition. 
    #     return [CLfactor, CDfactor]


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
    
    # def updateDensity(self,weather):
    #     newtemperature = weather.weatherlist[0].params['temperature'] +273.15 #convert to kelvin
    #     oldtemperature = weather.params['temperaturesealevel'] #convert to kelvin
    #     # print("old temperature is",oldtemperature)
    #     # print("new temperature is",newtemperature)
    #     # print("Self.params['temperature'] is",self.params['temperature'])
    #     olddensity = weather.params['airdensity']
    #     # print("old density is",olddensity)
    #     newdensity = olddensity * ((oldtemperature-newtemperature)/oldtemperature + 1) #inverse relationship
    #     temperatureeffect = newdensity / olddensity
    #     return temperatureeffect


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

    # def updateDensity(self,weather): # 2D linear interpolation based on Yue (2017)
    #     relativehumiditylist = [0.0,25.0,50.0,70.0,90.0] #percentage
    #     temperaturelist = [15.0,20.0,25.0,30.0,35.0] #may need to convert this to kelvin

    #     #humidity effects only
    #     humidityarray = np.array([[1.0, 1.0, 1.0, 1.0, 1.0],\
    #                                 [0.99837, 0.99751, 0.99662, 0.99656, 0.99389],\
    #                                 [0.99673, 0.99585, 0.99409, 0.99055, 0.98866],\
    #                                 [0.99592, 0.99419, 0.99155, 0.98883, 0.98517],\
    #                                 [0.99429, 0.99252, 0.98902, 0.98625, 0.98080]]) # increasing humditiy in rows, temp in columns


    #     # temperature and humidity effects included
    #     # humidityarray = np.array([[1.0, 0.98286, 0.96653, 0.95020, 0.93551],\
    #     #                             [0.99837, 0.98041, 0.96327, 0.94694, 0.92980],\
    #     #                             [0.99673, 0.97878, 0.96082, 0.94122, 0.92490],\
    #     #                             [0.99592, 0.97714, 0.95837, 0.93959, 0.92163],\
    #     #                             [0.99429, 0.97551, 0.95592, 0.93714, 0.91755]]) # increasing humditiy in rows, temp in columns

    #     relativehumidity = weather.weatherlist[1].params['relativehumidity']
    #     # temperature = weather.params['temperaturesealevel']
    #     # temperature = 288.15 - 273.15 #convert to celcius - need to update this if temperature is altered
    #     temperature = weather.weatherlist[0].params['temperature'] #already to celcius
    #     # print("Temperature is",temperature)
    #     # determine 4 closest points before doing a 2D interpolation

    #     if temperature < temperaturelist[0]:
    #         temperature = temperaturelist[0]
    #         print("WARNING: Temperature is below interpolation bounds. Humidity effects are assumed to be as if the temperature is",temperaturelist[0],"degrees C")
    #     elif temperature > temperaturelist[-1]:
    #         temperature = temperaturelist[-1]
    #         print("WARNING: Temperature is below interpolation bounds. Humidity effects are assumed to be as if the temperature is",temperaturelist[0],"degrees C")
        
    #     if relativehumidity < 0.0 or relativehumidity > 100.0:
    #         raise(Exception("~~~~~ ERROR: relative humidity is outside of available bounds ~~~~~"))
    #     elif relativehumidity > relativehumiditylist[-1]:
    #         relativehumidity = relativehumiditylist[-1]
    #         print("WARNING: Humidity is above interpolation bounds. Humidity effects are assumed to be as if the temperature is",relativehumiditylist[-1],"%")

    #     counter = 0
    #     for temp in temperaturelist:
    #         if temperature == temp:
    #             idT1 = counter
    #             idT2 = counter
    #             break
    #         elif temperature < temp:
    #             idT1 = counter - 1
    #             idT2 = counter
    #             break
    #         else:
    #             counter += 1
    #     counter = 0
    #     for hum in relativehumiditylist:
    #         if relativehumidity == hum:
    #             idH1 = counter
    #             idH2 = counter
    #             break
    #         elif relativehumidity < hum:
    #             idH1 = counter - 1
    #             idH2 = counter
    #             break
    #         else:
    #             counter += 1

    #     #interpolate temperature at each humidity index
    #     x1 = temperaturelist[idT1]
    #     x2 = temperaturelist[idT2]
    #     y1 = humidityarray[idH1,idT1]
    #     y2 = humidityarray[idH1,idT2]
    #     x = temperature
    #     tempid1 = interpolate(x1,x2,y1,y2,x)
        
    #     x1 = temperaturelist[idT1]
    #     x2 = temperaturelist[idT2]
    #     y1 = humidityarray[idH2,idT1]
    #     y2 = humidityarray[idH2,idT2]
    #     x = temperature
    #     tempid2 = interpolate(x1,x2,y1,y2,x)

    #     #now interpolate based on humidity
    #     x1 = relativehumiditylist[idH1]
    #     x2 = relativehumiditylist[idH2]
    #     y1 = tempid1
    #     y2 = tempid2
    #     x = relativehumidity
    #     humidityeffect = interpolate(x1,x2,y1,y2,x)
    #     return humidityeffect


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


class Ice(WeatherType):
    '''
    Class used to define icing characteristics. Params variables include:
    * variable (parameter) [units] {example}

    Ice is expected to affect:
    * lift
    * drag
    NOTE: Ice is very dangerous and difficult to model
    '''

    def __init__(self):
        pass
    
    # def warning(self,weather):
    #     for weatherclass in weather.weatherlist:
    #         if 'temperature' in weatherclass.params:
    #             temperature = weatherclass.params['temperature']
    #         elif 'relativehumidty' in weatherclass.params:
    #             relativehumidity = weatherclass.params['relativehumidity']
        
    #     if temperature < 5:
    #         if relativehumidity <= 0:
    #             print("WARNING: Cold temperatures may result in dangerous icing conditions. Avoid routes near clouds and/or rain.")
    #         elif relativehumidity < 30:
    #             print("WARNING: Cold temperatures and moderate humidity could result in dangerous icing conditions. Be cautious when flying.")
    #         else:
    #             print("WARNING: Cold temperatures and high humidity will likely result in dangerous icing conditions. Take extreme caution if flying.")


class Mission:
    'Class used to define the mission, including flight trajectory, maneuvers, etc.'

    params = {
        'payload': 0.0,
        'missionspeed': 10.0,
        'altitude': 5000
    }

    log     = {
        'SUCCESS': [],
        'WARNING': [],
        'ERROR': [],
        'NOTE': []
    }
    debug = False

    def __init__(self, params, debug=False):
        self.params = params
        self.debug  = debug

    def update(self):
        pass

    def __updateLog(self,entrytype,message):
        self.log[entrytype].append(message)
        # entrytype can be `SUCCESS`, `WARNING`, or `ERROR`


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

    log         = {
        'SUCCESS': [],
        'WARNING': [],
        'ERROR': [],
        'NOTE': []
    }

    debug = False

    # methods go here:
    def __init__(self, timestep, simulationtype, debug=False): #, desiredresult):
        # get parameters from `settings.txt`
        self.params['simulationtype'] = simulationtype
        self.debug = debug

    def run(self, drone, battery, power, weather, mission):
        if self.params['simulationtype'] == 'simple':
            # insert simple model here
            self.__runSimpleModel(drone, battery, power, weather, mission)
        # we'll need to define a list of these terms in the README
        elif self.params['simulationtype'] == 'complicated':
            # insert another model here
            pass
        else:
            self.__updateLog('ERROR','Simulation model not available.')
            raise(Exception("~~~~~ ERROR: simulation model not available ~~~~~"))

    def __updateLog(self,entrytype,message):
        self.log[entrytype].append(message)
        # entrytype can be `SUCCESS`, `WARNING`, or `ERROR`

class Plotter:
    'Class used to make plots of simulation results'

    # class variables go here:
    fig_num = 1

    # methods go here:
    # def __init__(self,x,xlabel: string,y,ylabel,axistitle) -> None:
    def __init__(self, x, xlabel, y, ylabel, axistitle, numplots):
        d = datetime.datetime.today()
        self.title = axistitle + " (" + d.strftime("%b-%d-%Y") + ")"
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.x = x
        self.y = y
        self.numplots = numplots
        # plt.plot(self.x,self.y)
        # plt.xlabel(self.xlabel)
        # plt.ylabel(self.ylabel)
        # plt.title(self.title)
        # plt.show()

    def plot_line(self):
        fig = plt.figure(Plotter.fig_num)
        fig.patch.set_facecolor('w')
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        plt.title(self.title)
        for i in range(self.numplots):
            plt.plot(self.x, self.y[i])
        Plotter.fig_num += 1
        plt.show()
        # input()

    def plot_scatter(self):
        fig = plt.figure(Plotter.fig_num)
        fig.patch.set_facecolor('w')
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        plt.title(self.title)
        for i in range(self.numplots):
            plt.plot(self.x, self.y[i], 'ro')
        Plotter.fig_num += 1
        plt.show()
        # input()

    def plot_validation(self, xvalid, yvalid):
        fig = plt.figure(Plotter.fig_num)
        fig.patch.set_facecolor('w')
        plt.xlabel(self.xlabel)
        plt.ylabel(self.ylabel)
        plt.title(self.title)
        plt.plot(xvalid, yvalid, 'ko')
        for i in range(self.numplots):
            plt.plot(self.x, self.y[i])
        Plotter.fig_num += 1
        plt.show()
        # input()


class model:
    '''
    Drone Applications in Weather Environments
    Top-level class used by pip3 package to run simulations using python dictionary inputs.
    `model.simulate()` returns `model.output` if `model.verbose=True`
    `model.simulate()` produces a plot if `model.plot=True`
    '''

    input       = {}
    output      = {}
    params      = {}
    classes     = {
        'drone': None,
        'weather': None,
        'power': None,
        'simulation': None,
        'mission': None,
        'battery': None,
    }
    plotter     = None
    correctunits= False

    verbose     = False # returns self.output after simulation
    plot        = False # toggles plotting feature
    debug       = {
        'drone': False,
        'weather': False,
        'battery': False,
        'mission': False,
        'model': False,
        'power': False
    }  # toggles debugging logs

    log         = {
        'SUCCESS': [],
        'WARNING': [],
        'ERROR': [],
        'NOTE': []
    }

    outputlog   = {
        'log': [],
        'xlabel': None,
        'zlabel': None
    }

    def __init__(self,input,verbose=False,plot=False,debug={}):
        self.input      = input
        self.verbose    = verbose
        self.plot       = plot
        for entry in debug:
            self.debug[entry] = debug[entry]

    def __resetParams(self):
        self.params = {
            "validation":True,
            "validationcase":"DiFranco2016",
            "drone":None,
            "dronename":"dji-Mavic2",
            "L/D":10.0, #we can change these default values vvv
            "propulsiveefficiency":0.35,
            "batterytechnology":"near-future",
            "stateofhealth":90.0,
            "startstateofcharge":100.0,
            "altitude": 5000,
            "dropsize":0.0,
            "liquidwatercontent":1.0, #one of this and rainfallrate needs to be specified for rain, but not both
            "rainfallrate":None,
            "temperature":15.0,
            "wind":False,
            "windspeed":10.0,
            "winddirection":0.0,
            "downdraft":0.0,
            "relativehumidity":0.0,
            "icing":False,
            "mission": {
                    "missionspeed":10.0,
                    "altitude":100.0,
                    "heading":0.0,
                    "payload":0.0
                },
            "timestep":1,
            # "plot":True,
            "xlabel":"missionspeed",
            "ylabel":"power",
            "title":"First_test",
            "simulationtype":"simple",
            "model":"abdilla",
            # "xbegin":0,
            # "xend":1,
            # "xnumber":5,
            "xvals":[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15],
            "zlabel":"temperature",
            # "weatherbegin":10,
            # "weatherend":40,
            # "weathernumber":3
            "zvals":[0,10,20,30,40] # must contain only unique elements
        }

        for key in self.input:
            self.params[key] = self.input[key]

    def __setupValidation(self):
        validationdata = getParams(validationdatabase,self.params['validationcase'])
        self.params             = validationdata['settings']
        self.params['xvalid']   = validationdata['xvalid']
        self.params['yvalid']   = validationdata['yvalid']
        self.classes['drone']              = Drone(validationdata['settings']['dronename'],validationdata['drone'],conversions,self.correctunits,self.debug['drone'])
        self.correctunits       = True

    def __setupSimulation(self):
        droneparams             = getParams(drones,self.params['dronename'])
        self.classes['drone']   = Drone(self.params['dronename'],droneparams,conversions,self.correctunits,debug=self.debug['drone'])
        self.correctunits       = True

    def __prepareSimulation(self):
        self.__resetParams()
        if self.params['validation'] == True:
            # validation = True
            # validationcase = self.params['validationcase']
            # validationdata = getParams(validationdatabase,validationcase) #specifies settings_list is in separate path
            # self.params = validationdata['settings']
            self.__setupValidation()
        else:
            self.__setupSimulation()

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
                                "missionspeed"
                                ]
        if self.params['xlabel'] not in independentvariables:
            self.__updateLog('ERROR','Desired x variable is not independent.')
            raise(Exception("~~~~~ ERROR: desired x variable is not independent ~~~~~"))

        weatherlist         = []

        # instantiate battery
        stateofhealth       = self.params['stateofhealth']
        startstateofcharge  = self.params['startstateofcharge']
        batterytechnology   = self.params['batterytechnology']
        self.classes['battery']        = Battery(self.classes['drone'],stateofhealth,startstateofcharge, batterytechnology,debug=self.debug['battery'])

        # instantiate mission
        missionparams       = self.params['mission']
        self.classes['mission']        = Mission(missionparams,debug=self.debug['mission'])

        # Temperature
        if self.params['xlabel'] == 'temperature':
            newtemperature  = self.params['xvals'][0]
        else:
            newtemperature  = self.params['temperature']
        temperatureparams   = {'temperature':newtemperature}        # Ampere-hours
        temperature     = Temperature(temperatureparams)
        weatherlist.append(temperature)

        # Humidity
        if self.params['xlabel'] == 'humidity':
            relativehumidity = self.params['xvals'][0]
        else:
            relativehumidity = self.params['relativehumidity']
        humidityparams      = {'relativehumidity':relativehumidity}
        humidity            = Humidity(humidityparams)
        weatherlist.append(humidity)

        # Rain
        dropsize            = self.params['dropsize']
        if self.debug['model']:
            print("MODEL: ----- Dropsize = ",dropsize)
        if self.params['liquidwatercontent'] != None:
            liquidwatercontent  = self.params['liquidwatercontent']
            if self.debug['model']:
                print("MODEL: ----- LWC = ",liquidwatercontent)
            rainparams = {'dropsize':dropsize, 'liquidwatercontent':liquidwatercontent}
        elif self.params['rainfallrate'] != None:
            rainfallrate = self.params['rainfallrate']
            rainparams = {'dropsize':dropsize, 'rainfallrate':rainfallrate}
        # …
        rain                = Rain(rainparams)
        weatherlist.append(rain)

        # Wind
        #     speed       = self.params['windspeed']
        #     direction   = self.params['winddirection']
        #     wind        = Wind(speed,direction)
        #     weatherlist.append(wind)

        # Icing
        #     weatherlist.append('icing')
        #     icing   = Icing()

        # print("Weather parameters are: ")
        # print(str(weatherlist)[1:-1]) 

        # weatherparams   = []
        # for weatherclass in weatherlist:
        #     weatherparams = weatherparams + weatherclass.params

        self.classes['weather']        = Weather(self.classes['mission'].params['altitude'],weatherlist,debug=self.debug['weather'])
        if self.debug['model']:
            print("Preparing to update weather:")
        self.classes['weather'].update(self.classes['drone'])
        if self.debug['model']:
            print("Weather updated.")
        self.classes['drone'].update(self.classes['weather'])
        if self.debug['model']:
            print("MODEL: ----- Preparing power class:")
            print("MODEL: ----- drone is",self.classes['drone'])
            print("MODEL: ----- weather is",self.classes['weather'])
            print("MODEL: ----- mission is",self.classes['mission'])
        self.classes['power']          = Power(self.classes['drone'],self.classes['weather'],self.classes['mission'], debug=self.debug['power'])
        if self.debug['model']:
            print("MODEL: ----- power class instantiated.")
        self.classes['simulation']     = Simulation(self.params['timestep'],self.params['simulationtype'])
        
        if 'zlabel' not in self.params:
                    self.params['zlabel']   = None
                    self.params['zvals']    = [0]

    def simulate(self):
        # self.__clearLogs()
        self.__prepareSimulation()
        x               = self.params['xvals']
        y               = []
        yplot = []

        for zvalue in self.params['zvals']:
            if self.params['zlabel']:
                for myclass in self.classes:
                    if self.params['zlabel'] in self.classes[myclass].params:
                        self.classes[myclass].params[self.params['zlabel']] = zvalue
                        foundz = True
                        break
                if foundz == False:
                    self.__updateLog('ERROR','Manipulated z variable not set.')
                    raise(Exception("MODEL: ~~~~~ ERROR: manipulated z variable not set"))
                else:
                    self.__updateLog('SUCCESS','Manipulated z variable set.')
                    if self.debug['model']:
                        print("MODEL: ===== SUCCESS: manipulated z variable set")
                        print("MODEL: ----- missionspeed = ",self.classes['mission'].params['missionspeed'])

            # i = 0
            # for weatherclass in self.classes['weather'].weatherlist:
            #     if weathereffect in weatherclass.params:
            #         self.classes['weather'].weatherlist[i].params[weathereffect] = zvalue
            #         # self.classes['weather'].update(self.classes['drone'])
            #         # self.classes['drone'].update(self.classes['weather'])
            #         # self.classes['power'].update(self.classes['drone'],self.classes['weather'],self.classes['mission'])
            #         # self.classes['battery'].update()
            #     i += 1

        # if weathereffect == 'undefined':
        #     pass
        # elif weathereffect == 'temperature':
        #     # print("weathereffect = temperature confirmed")
        #     self.weather.weatherlist[0].params["temperature"] = zvalue
        #     self.weather.update(self.classes['drone']) #splitting up so as to only update drone class if rain occurs (which happens after weather.update)
        # elif weathereffect == 'relativehumidity':
        #     self.weather.weatherlist[1].params["relativehummidity"] = zvalue
        #     self.weather.update(self.classes['drone'])
        # elif weathereffect == 'dropsize' or weathereffect == 'liquidwatercontent':
        #     self.weather.weatherlist[2].params[weathereffect] = zvalue
        #     self.weather.update(self.classes['drone'])
        #     self.classes['drone'].update(self.weather)
        # else:
        #     raise(Exception("~~~~~ ERROR: weathereffect not a valid input ~~~~~"))
        # self.classes['power'].update(self.classes['drone'],self.weather,self.classes['mission'])
        # self.classes['battery'].update()
        
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
        # elif ylabel in self.params:
        #     y.append(self.params[ylabel])
        # else:
        #     raise(Exception("~~~~~ ERROR: desired y variable not found ~~~~~"))

            foundx = False
            for xvalue in x:
                for myclass in self.classes:
                    if self.params['xlabel'] in self.classes[myclass].params:
                        self.classes[myclass].params[self.params['xlabel']] = xvalue
                        foundx = True
                        break
            
                if foundx == False:
                    raise(Exception("MODEL: ~~~~~ ERROR: manipulated x variable not set"))
                else:
                    if self.debug['model']:
                        print("MODEL: ===== SUCCESS: manipulated x variable set")
                        print("MODEL: ----- missionspeed = ",self.classes['mission'].params['missionspeed'])

                # update value
                ## determine x location
                # if xlabel in self.classes['drone'].params:
                #     self.classes['drone'].params['xlabel'] = xvalue
                #     self.classes['power'].update(self.classes['drone'],self.weather,self.classes['mission'])
                #     self.classes['battery'].update()
                # elif xlabel in self.classes['weather'].params:
                #     i = 0
                #     for weatherclass in self.classes['weather'].weatherlist:
                #         if weathereffect in weatherclass.params:
                #             self.classes['weather'].weatherlist[i].params[weathereffect] = xvalue
                #             self.classes['weather'].update(self.classes['drone'])
                #             self.classes['drone'].update(self.classes['weather'])
                #             self.classes['power'].update(self.classes['drone'],self.classes['weather'],self.classes['mission'])
                #             self.classes['battery'].update()
                #         i += 1
                # elif xlabel in self.classes['mission'].params:
                #     self.classes['mission'].params[xlabel] = xvalue
                #     self.classes['power'].update(self.classes['drone'],self.classes['weather'],self.classes['mission'])
                #     self.classes['battery'].update()
                # elif xlabel in self.params:
                #     self.params[xlabel] = xvalue
                #     self.classes['power'].update(self.classes['drone'],self.classes['weather'],self.classes['mission'])
                #     self.classes['battery'].update()

                self.classes['battery'].update(self.classes['drone'])
                self.classes['weather'].update(self.classes['drone'])
                self.classes['mission'].update()
                self.classes['power'].update(self.classes['drone'],self.classes['weather'],self.classes['mission'])
                if self.debug['model']:
                    print("========================= ABOUT TO RUN ========================")
                # self.classes['simulation'].run(self.classes['drone'],self.classes['battery'],self.classes['power'],self.classes['weather'],self.classes['mission'])
                self.__runSimpleModel()
                if self.debug['model']:
                    print("=========================   FINISHED   ========================")

                sety = False
                for myclass in self.classes:
                    if self.params['ylabel'] in self.classes[myclass].params:
                        y.append(self.classes[myclass].params[self.params['ylabel']])
                        sety = True
                        break
                if sety == False:
                    self.__updateLog('ERROR','Dependent variable not set.')
                    raise(Exception("MODEL: ~~~~~ ERROR: dependent variable not set"))
                else:
                    if self.debug['model']:
                        print("MODEL: ===== SUCCESS: dependent variable set")
                
                self.__updateOutput([self.classes['drone'],self.classes['battery'],self.classes['power'],self.classes['weather'],self.classes['mission'],self.classes['simulation']],self.params['zvals'].index(zvalue))
                # print('MODEL: ----- xlabel is ',self.params['xlabel'])
                self.__updateOutputLog([self.classes['drone'],self.classes['battery'],self.classes['power'],self.classes['weather'],self.classes['mission'],self.classes['simulation']],xvalue,zvalue)

            yplot.append(y)
            y = []

        # print("x data includes:    ",xplot)
        # print("y data includes:    ",yplot)

        # print("")

        # print("EXE.py:      Independent variable is ",xlabel)
        # print("EXE.py:      Desired Result is       ",self.params['ylabel'])

        if self.params['zlabel']:   
            if self.debug['model']:
                print("EXE.py:      Z iterator is           ",self.params['zlabel'])
            self.output['zvals'] = self.params['zvals']
        else:
            self.output['zvals'] = "undefined"

        self.params['yplot'] = yplot
        if self.plot == True:
            if self.params['validation'] == True:
                self.validationPlot()
            else:
                self.simulationPlot()
        
        if self.verbose:
            return self.output
        else:
            pass

    def __runSimpleModel(self):

        self.params['endurance']    = self.classes['battery'].params['batterycapacity'] * self.classes['battery'].params['voltagemean'] / \
                                      self.classes['power'].params['power']  # simple endurance model
        self.params['range']        = self.params['endurance'] * self.classes['mission'].params["missionspeed"]

    def printWeatherClasses(self,name): # method for debugging
        print("")
        print("//////===     ",name,"     ===\\\\\\")
        print("//////===     Weather Classes     ===\\\\\\")
        print("///  drone:      ",self.classes['drone'].params)
        print("///  power:      ",self.classes['power'].params)
        print("///  weather:    ",self.classes['weather'].weatherlist)
        print("///  mission:    ",self.classes['mission'].params)
        print("///  simulation: ",self.classes['simulation'].params)
        print("///  plotter:    ",self.plotter.params)
        print("///  battery:    ",self.classes['battery'].params)
        print("//////////////////////////////////////////")

    def simulationPlot(self):
        xlabel      = self.params['xlabel']
        ylabel      = self.params['ylabel']
        yplot       = self.params['yplot']
        axistitle   = self.params['title']
        plotter     = Plotter(self.params['xvals'],xlabel,yplot,ylabel,axistitle,len(self.params['zvals']))
        plotter.plot_line()

    def validationPlot(self):
        xvalid = self.params['xvalid']
        yvalid = self.params['yvalid']

        xplot  = self.params['xvals']
        xlabel = self.params['xlabel']
        yplot  = self.params['yplot']
        ylabel = self.params['ylabel']
        axistitle = self.params['title'] + " Validation"
        if 'weathervals' not in self.params:
            self.params['weathervals'] = [0]
        plotter = Plotter(xplot,xlabel,yplot,ylabel,axistitle,len(self.params['weathervals']))
        plotter.plot_validation(xvalid,yvalid)
    
    def setDefaultParams(self):
        pass

    def __updateOutput(self,classes,zindex):
        for myclass in classes:
            for param in myclass.params:
                if param not in self.output and param != "model" and param != "powerModel":
                    self.output[param] = []
                    for index in range(len(self.params['zvals'])):
                        self.output[param].append([])
                if param != "model" and param != "powerModel":
                    self.output[param][zindex].append(myclass.params[param])
    
    def __updateLog(self,entrytype,message):
        self.log[entrytype].append(message)
        # entrytype can be `SUCCESS`, `WARNING`, or `ERROR`

    def __updateOutputLog(self,classes,xvalue,zvalue):
        #TODO: add function to set self.outputlog['xlabel'] and self.outputlog['zlabel']
        self.outputlog['log'].append({})
        self.outputlog['log'][-1]['xvalue'] = xvalue
        self.outputlog['log'][-1]['zvalue'] = zvalue

        for myclass in classes+[self]:
            for entrytype in myclass.log:
                if entrytype not in self.outputlog['log'][-1]:
                    self.outputlog['log'][-1][entrytype] = []
                self.outputlog['log'][-1][entrytype] += myclass.log[entrytype]