import drone_awe
a = drone_awe.model({'validation':True,'dropsize':0.000,'zlabel':'payload','zvals':[0,.1,.2],'liquidwatercontent':0.0},plot=True,debug={'drone':False,'battery':False,'power':True,'weather':False,'mission':False,'model':False})
a.simulate()
