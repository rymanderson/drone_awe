# testing wind math
import numpy as np

for dt in [0,30,60,90,120,150,180]:
    velocityvector = [5.0,0.0,0.0]      # starts blowing north at 5m/s
    theta          = dt*np.pi/180.0               # it will turn 30deg in the positive direction (towards the east)
    rotationmatrix = np.array([[np.cos(theta),-np.sin(theta),0.0],[np.sin(theta),np.cos(theta),0.0],[0.0,0.0,0.0]])
    velocitynew    = np.dot(rotationmatrix,velocityvector)

    print("with a ",str(dt)," deg ccw angle adjustment, the velocity vector is: ",velocitynew)