# name = "Mavic2"

# paramfilepath = "params/" + name

# from Mavic2 import params
# print(params)

params = {
        'wingtype':None,
        'TOW':None,
        'max_speed':None,
        'max_alt':None,
        'max_t':None,
        'max_t_hover':None,
        'max_tilt':None,
        'min_temp':None,
        'max_temp':None,
        'power_rating':None,	
        'batt_type':None,
        'batt_capacity':None,
        'batt_voltage':None,
        'batt_cells':None,
        'batt_energy':None,
        'batt_mass':None
        }

paramlist = open("paramlist.param","r").readlines()
paramfile = open("Mavic2" + ".param", "r").readlines() #open file with read privilege

specs = []
values = []

for line in paramfile:
    parts = line.split()
    specs.append(parts[0])
    values.append(parts[1])

for line in paramlist:
    line = line.strip()
    for spec in specs:
        if line == spec:
            if spec == "wingtype" or spec == "batt_type":
                params[spec] = values[specs.index(spec)]
            else:
                params[spec] = float(values[specs.index(spec)])

print(params)

# from __future__ import print_function
# import csv

# infile = "Mavic2.param"
# tsvfile = csv.reader(infile, delimiter=" ")

# lines = []
# for line in tsvfile:
#     lines.append(line)
# print("Col1", [line[0] for line in lines])

# f = open("Mavic2.param","r")
# csvReader = csv.reader(f)
# specs = []
# values = []

# for row in csvReader:
#     specs.append(row[0])
#     # values.append(row[1])

# print(specs)
# print(values)

# inFile = "Mavic2.param"
# cols = {}
# indexToName = {}
# delim = " "
# for lineNum, line in enumerate(inFile):
#     if lineNum == 0:
#         headings = line.split(delim)
#         i = 0
#         for heading in headings:
#             heading = heading.strip()
#             cols[i] = [heading]
#             indexToName[i] = i
#             i += 1
#     else:
#         cells = line.split(delim)
#         i = 0
#         for cell in cells:
#             cell = cell.strip()
#             cols[indexToName[i]] += [cell]
#             i += 1
                
# print(cols)
# print(indexToName)

# f = open("Mavic2.param","r")

# # with f as inf:
# reader = csv.reader(f,delimiter=" ")
# specs = list(zip(*reader))[0]

# # with f as inf:
# #     reader = csv.reader(inf,delimiter=" ")
# values = list(zip(*reader))[1]

# # print(specs)
# # print(values)

# with open('Mavic2.param') as f:
#     values = zip(*[line.split() for line in f])[1]

# cols, indexToName = getColumns("Mavic2.param", delim=" ", header=False)