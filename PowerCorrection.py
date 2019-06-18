class PowerCorrection:
        'Class used to make corrections to the drone\'s power requirement based on maneuvers, weather, payload, etc.'

        # class variables go here:
        addedpower      = 0.0           # watts

        # methods go here:
        def __init__(self, sink, value):
            # sink is an object describing the reason for an increase/decrease in the power requirement (e.g., an instance of the Rain.py class)
            if sink == 'altitude':
                PowerCorrection.altitude(value)
            elif sink == 'payload':
                PowerCorrection.payload(value)
            elif sink == 'rain':
                PowerCorrection.rain(value)
            elif sink == 'temp':
                PowerCorrection.temp(value)
            elif sink == 'wind':
                PowerCorrection.wind(value)
            elif sink == 'humidity':
                PowerCorrection.humidity(value)
            elif sink == 'ice':
                PowerCorrection.ice()
            else:
                print("Error: Value for 'sink' is not recognized.")
        
        def altitude(self,new_alt):
            pass

        def rain(self,LWC):
            pass

        def temp(self,new_temp):
            del_temp = new_temp / old_temp
            del_rho = 1 / del_temp
            # del_rho = new_rho/old_rho
            del_power = 1 / sqrt(del_rho)
            new_power = old_power * del_power

        def wind(self,wind_velocity):
            #need to input heading
            #calculate wind in direction of heading (trig or dot product)
            #adjust velocity based on that - this method doesn't affect power consumption, only range
                #or we adjust power consumption to maintain the same velocity?
            new_velocity = old_velocity + wind_velocity #defined as in the same direction as drone heading

        def payload(self,new_payload):
            #simple model - add more weight and recalculate power
            pass

        def humidity(self,rel_hum):
            #we will want to create a best-fit line from Table 2 based on temperature and relative humidity
            
            #super simple version for now (assume 25 deg)
            if rel_hum < 25:
                del_power = 0.34 #percent
            elif rel_hum < 50:
                del_power = 0.59
            elif rel_hum < 70:
                del_power = 0.84
            else:
                del_power = 1.1
            
            new_power = old_power * (1 - del_power/100.0)

        def ice(self):
            #we may not even try to model this one
            print("Nope, you don't want to do that...")

print("Successfully imported `PowerCorrection.py`")