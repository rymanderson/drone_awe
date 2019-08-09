import drone_awe
a = drone_awe.model({'validation':False,'dronename':'aerovironment-Puma3AE'},plot=True,debug={'drone':False,'battery':False,'power':True,'weather':False,'mission':False,'model':False})
a.simulate()
