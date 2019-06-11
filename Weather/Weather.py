import numpy as np

class Weather:
        # 'Class used to predict the drone\'s power requirement'

        # class variables go here:
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
                #data based on lattitude, longitude
        #then calculate pressure, temp, and density based on altitude (see
                # eqns 1.21 - 1.22 in Dr. Ning's book

        # methods go here:
        def __init__(self, altitude): # keeping it simple to begin with
                self.altitude = altitude        
                self.temperature = self.temperature_sl - 71.5 + 2*np.log(1 + np.exp(35.75 - 3.25*self.altitude) + np.exp(-3 + 0.0003 * self.altitude**3))
                self.pressure = self.pressure_sl * np.exp(-0.118 * self.altitude - (0.0015*self.altitude**2) / (1 - 0.018*self.altitude + 0.0011 * self.altitude**2))


        def gettemperature(self):
                return 0

        print("Successfully imported `Weather.py`")
