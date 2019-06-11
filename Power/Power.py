import sys
sys.path.append('../')
import Drone

class Power:
        # 'Class used to predict the drone\'s power requirement'

        # class variables go here:
        basepower = 0.0   # watts
        mass = 1 #kg
        propulsive_eta = 0.8 #propulsive efficiency
        radius = .1 #m
        g = 9.806655 #m/s2
        Nr = 4 #number of rotors
        density = 1.225 #kg/m3

        # methods go here:
        def __init__(self,name,altitude):
                self.name = name
                self.drone = Drone.Drone(name)
                self.drone.getParams()
                self.Op1 = Weather.Weather(altitude)
                self.density = Op1.density
                self.mass = self.drone.params[TOW]
                if drone.params['wingtype'] == 'rotary':
                        self.Nr = self.drone.params[num_rotors]
                        self.L_D = 3
                        Power.getpower_rotor1()
                elif drone.params['wingtype'] == 'fixed':
                        self.L_D = 10
                        self.power = 0

        def getpower_rotor1(self): #super simple estimate for power
                power_elec = 0.1 #kW, estimate from paper
                self.power = self.mass * self.drone.param[cruise_speed] / (370 * self.propulsive_eta * self.L_D) - power_elec

        def getpower_rotor2(self): #slightly more complicated estimate for power
                self.power = self.mass * 1 / (self.propulsive_eta * self.drone.params[prop_diameter]/2.0) * np.sqrt(self.g**3 / (2 * self.Nr * self.density * np.pi))


print("Successfully imported `Power.py`")
