#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 10:32:17 2019

@author: mhumbert
"""

from numpy import f2py
import os

os.chdir('src')

f = open('calcdistances.f90','r').read()
f2py.compile(f,modulename='calcdistances',source_fn='calcdistances.f90',verbose=False)
f = open('calcCOM.f90','r').read()
f2py.compile(f,modulename='calccomf',source_fn='calcCOM.f90',verbose=False)
f = open('ipcorr.f90','r').read()
f2py.compile(f,modulename='ipcorr',source_fn='ipcorr.f90',verbose=False)

os.chdir('..')