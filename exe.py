import importer

# Dronename = 'dji-Mavic2'
# SoH = 90

params = functions.getParams('Simulation','settings_list.txt','settings'," ")
weather_list = []
    
#instantiate drone
drone = params['drone']
if params[drone] == True:
    dronename = params['dronename']
    drone = Classes.Drone(dronename)
else:
    raise Exception('Must specify drone name')

# instantiate battery
battery = Classes.Battery(drone)

#get class initialization info
# rainlength = 2 #for now
rain_test = params['rain']
if rain_test == True:
    weather_list.append('rain')
    dropsize = params['rain_dropsize']
    LWC = params['rain_LWC']
    # …
    rain = Classes.Rain(dropsize,LWC)
# else:
#     for i = 1:rainlength:
#         fp.readline()

temp_test = params['temperature']
if temp_test == True:
    weather_list.append('temperature')
    new_temp = params['new_temp']
    temperature = Classes.Temperature(new_temp)

wind_test = params['wind']
if wind_test == True:
    weather_list.append('wind')
    speed = params['wind_speed']
    direction = params['wind_direction']
    direction = Classes.Temperature(direction)

humid_test = params['humidity']
if humid_test == True:
    weather_list.append('humidity')
    rel_hum = params['rel_hum']
    direction = Classes.Humidity(rel_hum)

icing_test = params['icing']
if icing_test == True:
    weather_list.append('icing')
    icing = Classes.Icing()

print("Weather parameters are: ")
print(str(weather_list)[1:-1]) 

weather = Classes.Weather(weather_list)

#simulation variables
timestep = params['timestep']

simulation = Classes.Simulation(drone,battery,weather,power,timestep)
simulation.run()

if params['plot'] == True:
    xlabel = params['xlabel']
    ylabel = params['ylabel']
    axistitle = params['title']
    plotter = Classes.Plotter(simulation.x,xlabel,simulation.y,ylabel,axistitle)
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
