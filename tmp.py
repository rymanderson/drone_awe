import drone_awe
a = drone_awe.model({'validation':False,'rain':True,'dropsize':0.004,'liquidwatercontent':0.2},plot=True)
a.simulate()
