import numpy as np

validationdatabase = [
    {
        'id': 'Stolaroff2018',
        'xvalid': [1.358974359,1.7179487179,1.6923076923,1.7692307692,1.7692307692,1.8205128205,1.8717948718,1.8974358974,1.9230769231,1.9743589744,2.0512820513,2.1025641026,2.1025641026,2.1794871795,2.2820512821,2.3076923077,2.3333333333,2.4358974359,2.5384615385,2.5384615385,2.6666666667,2.7948717949,2.9487179487,3.1025641026,3.2564102564,3.358974359,3.4615384615,3.641025641,3.7435897436,3.8717948718,4.0512820513,4.2564102564,4.358974359,4.641025641,4.8974358974,5.8461538462,6,6.2820512821,6.5384615385,7.1538461538,7.4871794872,7.8205128205,8.4871794872,9.1794871795,9.7692307692,10.2564102564,10.6923076923,11.1282051282,11.5897435897],
        'yvalid': [133.7510729614,149.2375075128,142.6477385276,144.5897859156,140.9450186657,141.7583868756,142.3874173587,140.4340604922,139.3656195241,139.5335348046,140.2321151941,139.9589946754,136.1689649626,138.3618186589,141.3711557508,139.4001574523,136.7492021569,136.5131929807,137.712167001,133.136399421,135.7755034665,136.58240428,141.8694500174,143.7278953027,145.8751386173,144.4503136349,143.5492800366,144.4966689523,140.8944307591,139.7885398414,140.0842115956,140.6905892611,140.7593265104,140.4412051028,148.3297356325,161.0751453894,160.7738527567,158.8079335653,157.0439596719,161.2248774665,165.1381601781,164.7032531681,177.1408013138,175.6982333173,190.095233258,197.4485952036,209.0736216573,203.7684942987,215.7551870381],
        'drone': {
            'validationcase':'Stolaroff2018',
            'wingtype': 'rotary',
            'rotorquantity': 4,
            'takeoffweight': 1.3,
            'batterytype': 'LiPo',
            'batteryvoltage': 11.1,
            'batterymass': 0.262,
            'props': '10x4.7',
            'endurancemaxhover': 16,
            'payloadmax': 0.4,
            'batterycapacity': 5500,
            'payload': 0.0,
            'rotordiameter': 0.254,
            'batterycells': 3,
            'length': 0.280,
            'width': 0.140,
            'height': 0.100,
            'waterproof': 'no'
            },
        'settings': {
            'dronename': 'drone',
            'stateofhealth': 100.0,
            'startstateofcharge': 100.0,
            'altitude': 100.0,
            'temperaturesealevel': 15.0,
            'rain': False,
            'dropsize': 0.0,
            'liquidwatercontent': 1.0,
            'temperature': 15.0,
            'wind': False,
            'windspeed': 0.0,
            'winddirection': 0.0,
            'relativehumidity': 85.0,
            'icing': False,
            "mission": {
                "missionspeed": 10.0
            },
            'timestep': 1,
            'plot': True,
            'xlabel': 'missionspeed',
            'ylabel': 'power',
            'title': 'Stolaroff',
            'simulationtype': 'simple',
            'model': 'abdilla',
            'xvals': [0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0],
            # 'xbegin': 0.0,
            # 'xend': 15.0,
            # 'xnumber': 10,
            'validation': True,
            'validationcase': 'Stolaroff2018',
            'batterytechnology': 'current'
            }
    },
    {
        'id': 'Ostler2009',
        'xvalid': [10.6301249769294,14.0337688190471,15.8927017691882,17.1717378121127,18.2963324280855,19.4996308698289,20.8198974345453,21.6747633611939,22.5301236586073,23.3466264139004,23.930478287236,24.553286576845],
        'yvalid': [8.92583120204603,10.153452685422,12.6854219948849,15.2173913043478,17.9028132992327,21.6624040920716,26.1892583120205,29.5652173913043,33.3248081841432,36.9309462915601,40,43.2992327365729],
        'drone': {
            'validationcase':'Ostler2009',
            'wingtype': 'fixed',
            'wingarea': 0.321,
            'wingspan': 1.06,
            'rotorquantity': 1,
            'batteryvoltage': 11.1,
            'batterycapacity': 1500,
            'batterytype': 'LiPo',
            'numbatteries': 2,
            'numbatteriesconnection': 'parallel',
            'takeoffweight': 0.95,
            'props': '7x4',
            'rotordiameter': 0.1778,
            'payload': 0.0,
            'motorio': 0.65
            },
        'settings': {
            'validation': True,
            'validationcase': 'Ostler2009',
            'drone': True,
            'dronename': 'drone',
            'stateofhealth': 100.0,
            'startstateofcharge': 100.0,
            'altitude': 100.0,
            'temperaturesealevel': 15.0,
            'rain': False,
            'dropsize': 0.0,
            'liquidwatercontent': 1.0,
            'temperature': 15.0,
            'wind': False,
            'windspeed': 0.0,
            'winddirection': 0.0,
            'relativehumidity': 89.0,
            'icing': False,        
            "mission": {
                "missionspeed": 10.0
            },
            'timestep': 1,
            'plot': True,
            'xlabel': 'missionspeed',
            'ylabel': 'power',
            'title': 'Ostler2009',
            'simulationtype': 'simple',
            'xvals': np.linspace(0,30.0,20).tolist(),
            # 'xbegin': 0.0,
            # 'xend': 30.0,
            # 'xnumber': 20,
            'batterytechnology': 'current'
            }
    },
    {
        'id': 'Gatti2015',
        'batteryratio': [0.298319327731092,0.636134453781512,1.00420168067227,1.3672268907563,1.72521008403361,2.10336134453782,2.40588235294118,2.85966386554622,3.20252100840336,3.53025210084034,3.80252100840336],
        'endurance': [13.1838565022422,22.2085201793722,27.2533632286996,30.0560538116592,31.4013452914798,32.2421524663677,31.6255605381166,32.1860986547085,31.5695067264574,30.6165919282511,30.5044843049327],
        'drone': None,
        'settings': None
    },
    {
        'id': 'FreeflyAlta8',
        'xvalid': [0,0.4535929094,0.9071858189,1.3607787283,1.8143716377,2.2679645472,2.7215574566,3.175150366,3.6287432755,4.0823361849,4.5359290944,4.9895220038,5.4431149132,5.8967078227,6.3503007321,6.8038936415,7.257486551,7.7110794604,8.1646723698,8.6182652793,9.0718581887],
        'yvalid': [1421.052631579,1310.5263157895,1215.7894736842,1144.7368421053,1073.6842105263,1010.5263157895,947.3684210526,900,852.6315789474,805.2631578947,765.7894736842,726.3157894737,702.6315789474,663.1578947368,639.4736842105,607.8947368421,576.3157894737,568.4210526316,536.8421052632,513.1578947368,497.3684210526],
        'drone': {
            'validationcase':'FreeflyAlta8',
            'wingtype': 'rotary',
            'rotorquantity': 8,
            'diagonal': 1.325,
            'batterycells': 6,
            'batteryvoltage': 22.2,
            'batterycapacity': 10000,
            'batterytype': 'LiPo',
            'numbatteries': 2,
            'numbatteriesconnection': 'parallel',
            'max_takeoffweight': 18.1,
            'takeoffweight': 6.2,
            'max_payload': 9.1,
            'specific_power': 145,
            'props': '18x6_Folding',
            'motor_max_power_continuous': 350,
            'motor_max_power_peak': 950,
            'temperaturemin': -20,
            'temperaturemax': 45,
            'thrust_ratio_at_max_takeoffweight': '1.85:1',
            'rotordiameter': 0.4572,
            'payload': 0.0
            },
        'settings': {
            'drone': True,
            'dronename': 'drone',
            'stateofhealth': 100.0,
            'startstateofcharge': 100.0,
            'altitude': 100.0,
            'temperaturesealevel': 15.0,
            'rain': False,
            'dropsize': 0.0,
            'liquidwatercontent': 1.0,
            'temperature': 15.0,
            'wind': False,
            'windspeed': 0.0,
            'winddirection': 0.0,
            'relativehumidity': 89.0,
            'icing': False,
            "mission": {
                "missionspeed": 10.0
            },
            'timestep': 1,
            'plot': True,
            'xlabel': 'payload',
            'ylabel': 'endurance',
            'title': 'FreeFLY_Alta8',
            'simulationtype': 'simple',
            'xvals': np.linspace(0.0,10.0,100).tolist(),
            # 'xbegin': 0.0,
            # 'xend': 10.0,
            # 'xnumber': 100,
            'validation': True,
            'validationcase': 'FreeFLYAlta8',
            'batterytechnology': 'current'
            }
    },
    {
        'id': 'FireFLY6Pro',
        'xvalid': [0.232643087239583,0.27304114453125,0.3448008515625,0.400082403645833,0.53137608984375,0.612172204427083,0.67170618359375],
        'yvalid': [3533.27586206897,3417.00431034483,3296.27693965517,3233.64762931034,3174.38577586207,3057.70474137931,2995.03232758621],
        'drone': {
            'validationcase':'FireFLY6Pro',
            'wingtype': 'fixed',
            'span': 1.524,
            'takeoffweight': 4.5,
            'endurancemax': 59,
            'max_payload': 0.7,
            'cruise_speed': 18,
            'VTOL': 'yes',
            'batterytype': 'LiHV'
            },
        'settings': {
            'dronename': 'drone',
            'stateofhealth': 100.0,
            'startstateofcharge': 100.0,
            'altitude': 100.0,
            'temperaturesealevel': 15.0,
            'rain': False,
            'dropsize': 0.0,
            'liquidwatercontent': 1.0,
            'temperature': 15.0,
            'wind': False,
            'windspeed': 0.0,
            'winddirection': 0.0,
            'relativehumidity': 0.0,
            'icing': False,
            "mission": {
                "missionspeed": 10.0
            },
            'timestep': 1,
            'plot': True,
            'xlabel': 'payload',
            'ylabel': 'endurance',
            'title': 'FreeFLY_Alta8',
            'simulationtype': 'simple',
            'xvals': np.linspace(0.0,1.0,20).tolist(),
            # 'xbegin': 0.0,
            # 'xend': 1.0,
            # 'xnumber': 20,
            'validation': True,
            'validationcase': 'FireFLY6Pro',
            'batterytechnology': 'current'
            }
    },
    {
        'id': 'Dorling2017_4S',
        'xvalid': [0.551111111111111,0.551111111111111,0.673333333333333,0.673333333333333,0.793333333333333,0.8,0.924444444444444,0.913333333333333,1.03333333333333,1.04888888888889,1.16,1.15333333333333,1.14666666666667,1.28222222222222,1.4,1.39555555555556,1.51333333333333,1.51333333333333],
        'yvalid': [0.313180515759312,0.309169054441261,0.322922636103152,0.329799426934097,0.34297994269341,0.359598853868195,0.381375358166189,0.385959885386819,0.407736389684814,0.412893982808023,0.436389684813754,0.440974212034384,0.451862464183381,0.462750716332378,0.482234957020057,0.491404011461318,0.497134670487106,0.520630372492837],
        'drone': {
            'validationcase':'Dorling2017_4S',
            'wingtype': 'rotary',
            'rotorquantity': 6,
            'takeoffweight': 1.5,
            'batterytype': 'LiPo',
            'batteryvoltage': 14.8,
            'batterymass': 0.5,
            'props': '10x4.7'
            },
        'settings': {
            'dronename': 'drone',
            'stateofhealth': 90.0,
            'startstateofcharge': 100.0,
            'altitude': 100.0,
            'temperaturesealevel': 15.0,
            'rain': False,
            'dropsize': 0.0,
            'liquidwatercontent': 1.0,
            'temperature': 15.0,
            'wind': False,
            'windspeed': 0.0,
            'winddirection': 0.0,
            'relativehumidity': 0.0,
            'icing': False,
            "mission": {
                "missionspeed": 10.0
            },
            'timestep': 1,
            'plot': True,
            'xlabel': 'payload',
            'ylabel': 'power',
            'title': 'First_test',
            'simulationtype': 'simple',
            'model': 'abdilla',
            'xvals': np.linspace(0.5,1.5,10).tolist(),
            # 'xbegin': 0.5,
            # 'xend': 1.5,
            # 'xnumber': 10,
            'validation': True,
            'validationcase': 'Dorling2017_4S',
            'batterytechnology': 'current'
            }
    },
    {
        'id': 'Dorling2017_3S',
        'xvalid': [0.451111111111111,0.566666666666667,0.571111111111111,0.688888888888889,0.695555555555556,0.813333333333333,0.82,0.935555555555556,0.944444444444444,1.05777777777778,1.06,1.18222222222222,1.18444444444444,1.29555555555556,1.3,1.41111111111111,1.41333333333333],
        'yvalid': [0.268481375358166,0.283381088825215,0.286246418338109,0.307449856733524,0.301719197707736,0.326361031518625,0.32865329512894,0.344699140401146,0.337822349570201,0.36189111747851,0.360171919770774,0.385386819484241,0.386532951289398,0.403151862464183,0.412893982808023,0.428939828080229,0.443839541547278],
        'drone': {
            'validationcase':'Dorling2017_3S',
            'wingtype': 'rotary',
            'rotorquantity': 6,
            'takeoffweight': 1.5,
            'batterytype': 'LiPo',
            'batteryvoltage': 14.8,
            'batterymass': 0.5,
            'props': '10x4.7',
            'endurancemax': 15
            },
        'settings': {
            'dronename': 'drone',
            'stateofhealth': 90.0,
            'startstateofcharge': 100.0,
            'altitude': 100.0,
            'temperaturesealevel': 15.0,
            'rain': False,
            'dropsize': 0.0,
            'liquidwatercontent': 1.0,
            'temperature': 15.0,
            'wind': False,
            'windspeed': 0.0,
            'winddirection': 0.0,
            'relativehumidity': 0.0,
            'icing': False,
            "mission": {
                "missionspeed": 10.0
            },
            'timestep': 1,
            'plot': True,
            'xlabel': 'payload',
            'ylabel': 'power',
            'title': 'First_test',
            'simulationtype': 'simple',
            'model': 'abdilla',
            'xvals': np.linspace(0.5,1.5,10).tolist(),
            # 'xbegin': 0.5,
            # 'xend': 1.5,
            # 'xnumber': 10,
            'validation': True,
            'validationcase': 'Dorling2017_3S',
            'batterytechnology': 'current'
            }
    },
    {
        'id': 'DiFranco2016',
        'xvalid': [0.264285714285714,1.74285714285714,2.87857142857143,4.87142857142857,5.9,6.88571428571429,7.85,8.85714285714286,9.65,10.8071428571429,12.0071428571429,13.2928571428571,13.5071428571429,14.6214285714286,15.0285714285714],
        'yvalid': [220.714707868963,220.507429922323,219.820584262074,209.694360013509,204.750084430935,204.060283687943,209.989446133063,211.427727119216,216.880699088146,213.593802769335,239.385764944276,256.905184059439,279.131627828436,290.264690982776,332.116683552854],
        'drone': {
            'validationcase': 'DiFranco2016',
            'wingtype': 'rotary',
            'rotorquantity': 4,
            'takeoffweight': 1.3,
            'batterytype': 'LiPo',
            'batteryvoltage': 11.1,
            'batterymass': 0.262,
            'props': '10x4.7',
            'endurancemaxhover': 16,
            'payloadmax': 0.4,
            'batterycapacity': 5500,
            'payload': 0.0,
            'rotordiameter': 0.254,
            'batterycells': 3,
            'length': 0.280,
            'width': 0.140,
            'height': 0.100,
            'waterproof': 'no'
            },
        'settings': {
            'dronename': 'drone',
            'stateofhealth': 100.0,
            'startstateofcharge': 100.0,
            'altitude': 100.0,
            'temperaturesealevel': 15.0,
            'rain': False,
            'dropsize': 0.0,
            'liquidwatercontent': 1.0,
            'temperature': 15.0,
            'wind': False,
            'windspeed': 0.0,
            'winddirection': 0.0,
            'relativehumidity': 0.0,
            'icing': False,
            "mission": {
                "missionspeed": 10.0
            },
            'timestep': 1,
            'plot': True,
            'xlabel': 'missionspeed',
            'ylabel': 'power',
            'title': 'DiFranco',
            'simulationtype': 'simple',
            'model': 'abdilla',
            'xvals': np.linspace(0.0,16.0,20).tolist(),
            # 'xbegin': 0.0,
            # 'xend': 16.0,
            # 'xnumber': 20,
            'validation': True,
            'validationcase': 'DiFranco2016',
            'batterytechnology': 'current'
            }
    },
    {
        'id': 'Chang2016',
        'xvalid': [0.183691373832219,0.797979797979798,2.81338028169014,5.75096030729834,8.75045051453502,11.068347797221,12.9491392801252,15.9902072366861],
        'yvalid': [18.1818181818182,17.4242424242424,17.6515151515152,16.0606060606061,17.9545454545455,22.0454545454545,25.1515151515152,29.0909090909091],
        'drone': {
            'validationcase':'Chang2016',
            'wingtype': 'rotary',
            'rotorquantity': 4,
            'takeoffweight': 1.3,
            'batterytype': 'LiPo',
            'batteryvoltage': 11.1,
            'batterymass': 0.262,
            'props': '10x4.7',
            'endurancemax': 15,
            'payloadmax': 0.4,
            'batterycapacity': 5500,
            'payload': 0.095,
            'rotordiameter': 0.254
            },
        'settings': {
            'dronename': 'drone',
            'stateofhealth': 100.0,
            'startstateofcharge': 100.0,
            'altitude': 100.0,
            'temperaturesealevel': 15.0,
            'rain': False,
            'dropsize': 0.0,
            'liquidwatercontent': 1.0,
            'temperature': 15.0,
            'wind': False,
            'windspeed': 0.0,
            'winddirection': 0.0,
            'relativehumidity': 0.0,
            'icing': False,
            "mission": {
                "missionspeed": 10.0
            },
            'timestep': 1,
            'plot': True,
            'xlabel': 'missionspeed',
            'ylabel': 'power',
            'title': 'Chang2016',
            'simulationtype': 'simple',
            'model': 'abdilla',
            'xvals': np.linspace(0.0,16.0,20).tolist(),
            # 'xbegin': 0.0,
            # 'xend': 16.0,
            # 'xnumber': 20,
            'validation': True,
            'validationcase': 'Chang2016',
            'batterytechnology': 'current'
            }
    },
    {
        'id': 'Abdilla2015endurance',
        'xvalid': [0.409590526228004,0.415107163706113,0.437878408806605,0.440750006949657,0.442819614711033,0.450689405943346,0.475877741639563,0.504592333138743,0.505405442971117,0.549256386734495,0.55135379312262],
        'yvalid': [4.94708994708995,4.94708994708995,8.12169312169312,8.80952380952381,10.026455026455,10.4497354497354,11.6137566137566,18.5449735449735,17.5925925925926,18.3333333333333,18.4920634920635],
        'drone': {
            'validationcase':'Abdilla2015endurance',
            'wingtype': 'rotary',
            'rotorquantity': 4.0,
            'endurancemax': 12.0,
            'rangemax': 50.0,
            'takeoffweight': 1.81,
            'batterycapacity': 2200.0,
            'batteryvoltage': 11.1,
            'batterytype': 'LiPo',
            'batterymass': 0.191,
            'payload': 0.0,
            'rotordiameter': 0.2
            },
        'settings': {
            'dronename': 'drone',
            'stateofhealth': 90.0,
            'startstateofcharge': 100.0,
            'altitude': 100.0,
            'temperaturesealevel': 15.0,
            'rain': False,
            'dropsize': 0.0,
            'liquidwatercontent': 1.0,
            'temperature': 15.0,
            'wind': False,
            'windspeed': 0.0,
            'winddirection': 0.0,
            'relativehumidity': 0.0,
            'icing': False,
            "mission": {
                "missionspeed": 10.0
            },
            'timestep': 1,
            'plot': True,
            'xlabel': 'payload',
            'ylabel': 'endurance',
            'title': 'Abdilla 2015 Endurance vs Payload Validation Test',
            'simulationtype': 'simple',
            'model': 'abdilla',
            'xvals': np.linspace(0.4,0.55,20).tolist(),
            # 'xbegin': 0.4,
            # 'xend': 0.55,
            # 'xnumber': 20,
            'validation': True,
            'validationcase': 'Abdilla2015endurance',
            'batterytechnology': 'current'
            }
    }
]