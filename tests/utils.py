#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015-2017, Tomas Tinoco De Rubira.    #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

import os
import csv
from os import listdir
from os.path import join

DIR_CASES = join('.', 'tests', 'resources', 'cases')
DIR_PFSOL = join('.', 'tests', 'resources', 'pf_solutions')

test_cases = [join(DIR_CASES, f) for f in listdir(DIR_CASES)]
test_cases.sort()

def get_pf_solution_file(case_file, dir, sol):
    """
    Gets power flow solution file path.

    Parameters
    ----------
    case_file : path to case file (string)
    dir : path to solution file directory (string)
    sol : solution code extension (string)

    Returns
    -------
    sol_file : path to solution file (string)
    """

    return join(dir, case_file.split(os.sep)[-1]+'.'+sol)

def read_pf_solution_file(sol_file):
    """
    Reads contents of power flow solution file.
    
    Parameters
    ----------
    sol_file : path to solution file.

    Returns
    -------
    sol_data : solution data (dictionary)
    """

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
