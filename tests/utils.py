#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

import csv
from os import listdir
from os.path import join

DIR = './tests/resources'
test_cases = [join(DIR+'/cases',f) for f in listdir(DIR+'/cases')]
test_cases.sort()

def get_pf_solution_file(case,sol):

    return join(DIR+'/pf_solutions',case.split('/')[-1]+'.'+sol)

def read_pf_solution_file(sol_file):

    try:

        bus_data = {}
        sol_data = {'v_mag_tol': 0.,
                    'v_ang_tol': 0.,
                    'bus_data': bus_data}

        f = open(sol_file)
        reader = csv.reader(f,delimiter=',')
        
        v_mag_tol,v_ang_tol = list(map(float,next(reader)))

        sol_data['v_mag_tol'] = v_mag_tol
        sol_data['v_ang_tol'] = v_ang_tol

        next(reader) # header

        for row in reader:
            bus_number,code,v_mag,v_ang = int(row[0]),int(row[1]),float(row[2]),float(row[3])
            bus_data[bus_number] = {'v_mag': v_mag, # p.u.
                                    'v_ang': v_ang} # degrees
           
        return sol_data

    except IOError:
        return None
