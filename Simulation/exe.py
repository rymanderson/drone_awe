import importer

# Dronename = 'dji-Mavic2'
# SoH = 90

with open('settings.txt','r') as fp:
    #instantiate drone
    drone = readline()
    if drone == true:
        dronename = fp.readline()
        drone = Drone(dronename)
    else:
        raise Exception('Must specify drone name')

    # instantiate battery
    battery = Battery(drone)

    #get class initialization info
        rainlength = 2 #for now
        rain = fp.readline()
        if rain == true:
            Dropsize = readline
            WVC = readline
            # …
            rain = Rain(dropsize,WVC,...)
        else:
            for i = 1:rainlength:
                fp.readline()

        #simulation variables
        Timestep = fp.readline()
        
        #plot settings
        ##time plots
        Totalpower	true
        Rainintensity false

        while eof == false:
        nextline = readline()
        if nextline == plot:
            plotx.push(fp.readline())
            ploty.push(fp.readline())
        else:
            raise Exception('Invalid plot syntax')

simulation = Simulation(drone,weather,power,...timestep,...)
simulation.run()

plotter = Plotter(inputs….)
plotter.plot()
