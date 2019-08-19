import numpy as np
import matplotlib.pyplot as plt

I0 = 1e-12 # W/m2
firstdistance = 10 # meters from first measurement

altitude = [
	1,
	2,
	3,
	4,
	5,
	7.6,
	10,
	15,
	20,
	25,
	30,
	35,
	40,
	# 45,
	# 50,
	# 60,
	# 70
]

altitudepredicted = altitude + np.linspace(80,100).tolist()

# for element in altitude:
# 	altitude[altitude.index(element)] -= 1
# altitude[0] = 0

decibels = [
	66.3,
	64.6,
	65.7,
	64.9,
	64.4,
	60.5,
	57.9,
	52.5,
	50.7,
	48.7,
	46.1,
	45.5,
	44.1,
	# 45.4,
	# 49.6,
	# 45.5,
	# 44.5
]

ambient = [
    37.9,
    37.9,
    37.9,
    37.9,
    37.9,
    37.9,
    37.9,
    37.9,
    37.9,
    37.9,
    37.9,
    37.9,
    37.9,
    # 37.9,
    # 37.9,
    # 37.9,
    # 37.9
]

# calculate drone power
ambientintensity    = 10**(ambient[0]/10) * I0
dronepower          = (10**(decibels[1]/10) * I0 - ambientintensity) * 4*np.pi*firstdistance

# calculate predicted intensity as a function of distance
droneintensity      = []
decibelspredicted   = []
for i in range(len(altitude+altitudepredicted)):
    droneintensity.append(dronepower / (4*np.pi*(altitude + altitudepredicted)[i]**2))
    decibelspredicted.append(10 * np.log10(droneintensity[i]/I0 + ambientintensity/I0))

fig = plt.figure(1)
fig.patch.set_facecolor('w')
plt.xlabel('Altitude [m]')
plt.ylabel('Noise Level [dB]')
plt.plot(altitude, decibels)
plt.plot(altitude, ambient)
# plt.plot(altitude+altitudepredicted, decibelspredicted)
plt.legend(['Drone Noise','Ambient Noise'])

plt.show()
