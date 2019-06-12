import sys
sys.path.append('../')
sys.path.append('./')
sys.path.append('../Drone/')
from Drone import Drone
from Power import Power
from Weather import Weather 

alt = 100

mav2 = Drone.Drone("dji-Mavic2")
mav2.Drone.getParams()
mav2_w = Weather(alt)
test = Power("dji-Mavic2",mav2,mav2_w)
