import Drone

mav2 = Drone.Drone("dji-Mavic2")
print("Drone name is: ", mav2.name)
mav2.getparams()
print("Battery type is: ", mav2.params['batt_type'])
