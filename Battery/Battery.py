class Battery:
    'Class used to track battery characteristics and performance during simulation'

    # class variables go here:
    capacity        = None          # Ampere-hours
    soc             = None          # state of charge (in percent nominal capacity)
    soh             = None          # state of health (actual capacity divided by ideal capacity)
    batterytype     = None          # possible values include LiPo, Li-ion, NiCd, NiMH, SLA
    voltage         = None          # Volts; this is the instantaneous voltage
    voltage_mean    = None          # Volts; this is the average voltage used for time-invariant simulations
    voltage_charged = None          # Volts
    voltage_dead    = None          # Volts
    current         = None          # Amperes; this is the instantaneous current

    # methods go here:
    def __init__(self, drone):                      # class constructor
        # estimate list lengths for prior memory allocation
        self.batterytype = drone.battery['type']
        self.voltage     = drone.battery['voltage']

    # time-invariant methods go here (indicated by suffix *_ti):
    def getEndurance_ti(self, power):
        
        # check that objects are time-invariant
        if power.timevariant == True:
            print("")
            sys.exit("ERR: attempted to run a time-invariant simulation with a time-variant `Power` object")
        else:

            # calculate endurance
            return self.capacity * self.voltage_mean / power.power

    # time-variant methods go here:
    def updateLoad(self, power):
        self.current    = power.power / self.voltage

    def discharge(self, power, timestep):
        self.soc        = self.soc - power.power * timestep
        self.current    = power.power / self.voltage



    

                

print("Successfully imported `Battery.py`")

