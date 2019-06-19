import importer

# Dronename = 'dji-Mavic2'
# SoH = 90

sim_params = functions.getParams('Simulation','settings_list.txt','settings.txt'," ")
weather_list = []
    
#instantiate drone
if sim_params['drone'] == True:
    dronename = sim_params['dronename']
    drone_params = functions.getPArams('Drone','paramlist.param',dronename + '.param',' ')
    drone = classes.Drone(dronename,drone_params)
else:
    raise Exception('Must specify drone name')

# instantiate battery
soh = sim_params['soh']
start_soc = sim_params['start_soc']
battery = classes.Battery(drone,soh,start_soc)

#get class initialization info
# rainlength = 2 #for now
rain_test = sim_params['rain']
if rain_test == True:
    weather_list.append('rain')
    dropsize = sim_params['rain_dropsize']
    LWC = sim_params['rain_LWC']
    # …
    rain = classes.Rain(dropsize,LWC)
# else:
#     for i = 1:rainlength:
#         fp.readline()

temp_test = sim_params['temperature']
if temp_test == True:
    weather_list.append('temperature')
    new_temp = sim_params['new_temp']
    temperature = classes.Temperature(new_temp)

wind_test = sim_params['wind']
if wind_test == True:
    weather_list.append('wind')
    speed = sim_params['wind_speed']
    direction = sim_params['wind_direction']
    direction = classes.Temperature(direction)

humid_test = sim_params['humidity']
if humid_test == True:
    weather_list.append('humidity')
    rel_hum = sim_params['rel_hum']
    direction = classes.Humidity(rel_hum)

icing_test = sim_params['icing']
if icing_test == True:
    weather_list.append('icing')
    icing = classes.Icing()

print("Weather parameters are: ")
print(str(weather_list)[1:-1]) 

weather = classes.Weather(weather_list)

#simulation variables
timestep = sim_params['timestep']

power = classes.Power(drone,battery,weather)

simulation = classes.Simulation(drone,battery,weather,power,timestep)
simulation.run()

if sim_params['plot'] == True:
    xlabel = sim_params['xlabel']
    ylabel = sim_params['ylabel']
    axistitle = sim_params['title']
    plotter = classes.Plotter(simulation.x,xlabel,simulation.y,ylabel,axistitle)
    plotter.plot_line()

# with open('textfile.txt','r') as fp:
#     #instantiate drone
#     drone = readline()
#     if drone == true:
#         dronename = fp.readline()
#         drone = Drone(dronename)
#     else:
#         raise Exception('Must specify drone name')

#     # instantiate battery
#     battery = Battery(drone)

#     #get class initialization info
#         rainlength = 2 #for now
#         rain = fp.readline()
#         if rain == true:
#             Dropsize = fp.readline()
#             WVC = fp.readline()
#             # …
#             rain = Rain(dropsize,WVC,...)
#         else:
#             for i = 1:rainlength:
#                 fp.readline()

#         #simulation variables
#         Timestep = fp.readline()
#         # …
#         #plot settings
#         ##time plots
#         Totalpower	true
#         Rainintensity false
#         …
#         # VVV In text file itself VVV
#         # ##time invariant plots: select two parameters from list:
#         # ### range payload endurance …
#         # #### example:
#         # #### range	payload
#         # plot range payload
#         # ^^^ 		         ^^^
#         while eof == false:
#         nextline = readline()
#         if nextline == plot:
#             plotx.push(fp.readline())
#             ploty.push(fp.readline())
#         else:
#             raise Exception('Invalid plot syntax')

# simulation = Simulation(drone,weather,power,...timestep,...)
# simulation.run()

# plotter = Plotter(inputs….)
# plotter.plot()
