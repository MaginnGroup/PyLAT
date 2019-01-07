# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 14:45:35 2015

@author: mhumbert
"""

import os
import random


samplefile = open('test_1/in.visc','r').readlines()
filelist= range(2,51)
for num in filelist:
    os.system('mkdir test_{0}'.format(num))
    output = open('test_{0}/in.visc'.format(num),'w')
    for line in range(0,len(samplefile)):
        if line == 31:
            output.write('velocity        all create  ${{mytemp}} {0} units box \n'.format(random.randint(1,999999999)))
        else:
            output.write(samplefile[line])
    output.close()
    os.system('cp test_1/mol.data test_{0}/'.format(num))
