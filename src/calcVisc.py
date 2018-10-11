# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 09:16:20 2015

@author: mhumbert
"""

from viscio import LammpsLog
import numpy as np
import sys

class calcVisc:
    
    def calcvisc(self,numtrj,numskip,dirbase,logname,output,ver):
        '''
        
        Calculates average and standard deviation of the integral of the 
        pressure tensor autocorrelation function over numtrj lammps trajectories
        
        '''
        output['Viscosity']={}
        output['Viscosity']['Units']='cP'
        if dirbase==None:
            dirbase='./'
        filename=dirbase+'1/'+logname
        Log = LammpsLog.from_file(filename)
        (Time,visco)=Log.viscosity(numskip)
        trjlen = len(Time)
        viscosity = np.zeros((numtrj,trjlen))
        for i in range(0,len(visco)):
            viscosity[0][i] += visco[i]
        if ver>=1:
            sys.stdout.write('Viscosity Trajectory 1 of {} complete'.format(numtrj))
        
        
        for i in range(2,numtrj+1):
            filename=dirbase+str(i)+'/'+logname
            Log = LammpsLog.from_file(filename)
            (Time,visco)=Log.viscosity(numskip)
            if len(visco) < trjlen:
                trjlen = len(visco)
            for j in range(0,trjlen):
                viscosity[i-1][j] += visco[j]
            if ver>=1:
                sys.stdout.write('\rViscosity Trajectory {} of {} complete'.format(i,numtrj))
        if ver>=1:
            sys.stdout.write('\n')
        
        
        average = np.zeros(trjlen)
        stddev = np.zeros(trjlen)
        for i in range(0,trjlen):
            average[i] = np.average(viscosity.transpose()[i])
            stddev[i] = np.std(viscosity.transpose()[i])
        
        output['Viscosity']['Average Integral']=average.tolist()
        output['Viscosity']['Standard Deviation']=stddev.tolist()
        output['Viscosity']['Time']=Time[:trjlen].tolist()
        
        return(output)
