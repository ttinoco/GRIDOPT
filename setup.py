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
      version='1.3.5rc1',
      description='Power Grid Optimization Library',
      url='https://github.com/ttinoco/GRIDOPT',
      author='Tomas Tinoco De Rubira',
      author_email='ttinoco5687@gmail.com',
      license='BSD 2-Clause License',
      packages=['gridopt',
                'gridopt.power_flow'],
      scripts=['./script/gridopt'],
      classifiers=['Development Status :: 5 - Production/Stable',
                   'License :: OSI Approved :: BSD License',
                   'Programming Language :: Python :: 2.7',
                   'Programming Language :: Python :: 3.5'],
      install_requires=['cython>=0.20.1',
                        'numpy>=1.11.2',
                        'scipy>=0.18.1',
                        'pfnet==1.3.3rc1',
                        'optalg==1.1.6rc1',
                        'nose'])
