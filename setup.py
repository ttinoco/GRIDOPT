#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015, Tomas Tinoco De Rubira.         #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

from setuptools import setup

setup(name='GRIDOPT',
      zip_safe=False,
      version='1.3.3',
      description='Power Grid Optimization Library',
      author='Tomas Tinoco De Rubira',
      author_email='ttinoco5687@gmail.com',
      license='BSD 2-Clause License',
      packages=['gridopt',
                'gridopt.power_flow'],
      scripts=['./script/gridopt'],
      install_requires=['numpy>=1.11.2',
                        'scipy>=0.18.1',
                        'pfnet',
                        'optalg'])
