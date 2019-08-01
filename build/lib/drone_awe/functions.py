# import required modules
import os


def getParams(params,case):
    for mydict in params:
        if mydict['id'] == case:
            return mydict
    raise(Exception("~~~~~ ERROR: CASE NOT FOUND ~~~~~"))
    
    # # initialize `params` dictionary
    # params = {}
    # # stringlist = ['wingtype','batterytype','waterproof','VTOL','drone','dronename','rain','temperature','wind','humidity','icing','plot','xlabel','ylabel','title','simulationtype']
    # # find name.param in the params/ directory

    # if listpath == "":  # if a different listpath is specified, use it instead. Otherwise use the default path
    #     listpath = "params/" + classname

    # pwd = os.path.dirname(os.path.abspath(__file__))

    # with open(os.path.join(pwd, listpath, listname), 'r') as paramlist:
    #     with open(os.path.join(pwd, 'params', classname, filename), 'r') as paramfile:
    #         paramlist = paramlist.readlines()
    #         paramfile = paramfile.readlines()  # open file with read privilege

    #         specs = []  # initialize lists
    #         values = []

    #         for line in paramfile:  # separates columns from .param file
    #             parts = line.strip().split(delimiter)
    #             specs.append(parts[0])
    #             values.append(parts[1])
    #         for line in paramlist:
    #             # gets rid of quotation marks when comparing strings
    #             line = line.strip().split(delimiter)
    #             for spec in specs:
    #                 if line[0] == spec:
    #                     # if values[specs.index(spec)].isnumeric():
    #                     # 	params[spec] = float(values[specs.index(spec)])
    #                     # else:
    #                     # 	if values[specs.index(spec)] == 'True' or values[specs.index(spec)] == 'true':
    #                     # 		params[spec] = True
    #                     # 	elif values[specs.index(spec)] == 'False' or values[specs.index(spec)] == 'false':
    #                     # 		params[spec] = False
    #                     # 	else:
    #                     # 		params[spec] = values[specs.index(spec)]

    #                     try:
    #                         params[spec] = float(values[specs.index(spec)])
    #                     except:
    #                         if values[specs.index(spec)] == 'True' or values[specs.index(spec)] == 'true':
    #                             params[spec] = True
    #                         elif values[specs.index(spec)] == 'False' or values[specs.index(spec)] == 'false':
    #                             params[spec] = False
    #                         else:
    #                             params[spec] = values[specs.index(spec)]

    #         if len(paramlist) > len(params):
    #             print("~~~~~ WARNING: ", str(len(paramlist)-len(params)),
    #                   " parameter(s) missing from '", filename, "' parameters file ~~~~~")


    return params

# reads in two columns into x and y lists, ignoring the top row (labels)
def getXandY(validationcase, delimiter):
    pwd = os.path.dirname(os.path.abspath(__file__))
    path = "params/Validation/" + validationcase
    filename = "data.txt"
    with open(os.path.join(pwd, path, filename), "r") as paramfile:
        paramfile = paramfile.readlines()  # open file with read privilege
        x = []  # initialize lists
        y = []
        for line in paramfile:  # separates columns from .param file
            parts = line.strip().split(delimiter)
            # print(parts)
            x.append(parts[0])
            y.append(parts[1])
        x = x[1:]
        y = y[1:]
        counter = 0
        for val in x:  # convert to floats
            x[counter] = float(val)
            counter += 1
        counter = 0
        for val in y:
            y[counter] = float(val)
            counter += 1
        return x, y

def interpolate(x1,x2,y1,y2,x):
    if x1 == x2:
        y = y1
    else:
        y = (x-x1)*(y2-y1) / (x2-x1) + y1   
    return y