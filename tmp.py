import drone_awe
a = drone_awe.model({'validation':True,'validationcase':'Stolaroff2018','dropsize':0.000,'zlabel':'payload','zvals':[0,.1,.2],'liquidwatercontent':0.0},plot=False,debug={'drone':False,'battery':False,'power':True,'weather':False,'mission':False,'model':False})
a.simulate()

#a = drone_awe.model({'validation':True,'validationcase':'Ostler2009','dropsize':0.000,'zlabel':'payload','zvals':[0,.1,.2],'liquidwatercontent':0.0},plot=False,debug={'drone':False,'battery':False,'power':True,'weather':False,'mission':False,'model':False})
#a.simulate()

a = drone_awe.model({'validation':False,'validationcase':'Ostler2009','dropsize':0.004,'zlabel':'payload','zvals':[0,.1,.2],'liquidwatercontent':0.0},plot=False,debug={'drone':False,'battery':False,'power':True,'weather':False,'mission':False,'model':False})
a.simulate()

