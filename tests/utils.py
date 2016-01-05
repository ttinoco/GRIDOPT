#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2016, Tomas Tinoco De Rubira.         #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

import csv
from os import listdir
from os.path import join

# Test cases
RES_DIR = './tests/resources'
files = [join(RES_DIR,f) for f in listdir(RES_DIR)]
test_cases = [f for f in files if f.split('.')[-1] not in ['sol1','sol2']]

def read_solution_data(sol_file):

    bus_data = {}
    sol_data = {'v_mag_tol': 0.,
                'v_ang_tol': 0.,
                'bus_data': bus_data}

    try:

        f = open(sol_file)
        reader = csv.reader(f,delimiter=',')
        
        v_mag_tol,v_ang_tol = map(float,reader.next())

        sol_data['v_mag_tol'] = v_mag_tol
        sol_data['v_ang_tol'] = v_ang_tol

        reader.next() # header

        for row in reader:
            bus_number,code,v_mag,v_ang = int(row[0]),int(row[1]),float(row[2]),float(row[3])
            bus_data[bus_number] = {'v_mag': v_mag, # p.u.
                                    'v_ang': v_ang} # degrees
            
    except IOError:
        pass

    return sol_data

