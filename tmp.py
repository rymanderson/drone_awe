import drone_awe
a = drone_awe.model({'debug':True,'validation':False,'rain':True,'dropsize':0.004,'zlabel':'payload','zvals':[0,.1,.2],'liquidwatercontent':0.2},plot=False)
a.simulate()
