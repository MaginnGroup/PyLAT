#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 10:32:17 2019
Modified by Rohit Goswami (HaoZeke) <rog32@hi.is>

@author: mhumbert
"""

from numpy import f2py
import os, glob

os.chdir('src')

f = open('calcdistances.f90','r').read()
f2py.compile(f,modulename='calcdistances',source_fn='calcdistances.f90',verbose=False)
f = open('calcCOM.f90','r').read()
f2py.compile(f,modulename='calccomf',source_fn='calcCOM.f90',verbose=False)
f = open('ipcorr.f90','r').read()
f2py.compile(f,modulename='ipcorr',source_fn='ipcorr.f90',verbose=False)

os.symlink(glob.glob('ipcorr.*.so')[0],'ipcorr.so')
os.symlink(glob.glob('calcdistances.*.so')[0],'calcdistances.so')
os.symlink(glob.glob('calccomf.*.so')[0],'calccomf.so')

os.chdir('..')
