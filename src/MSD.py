# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 14:08:27 2015

@author: mhumbert

PyLAT: Python LAMMPS Analysis Tools
Copyright (C) 2018  Michael Humbert, Yong Zhang and Ed Maginn

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""

import numpy as np
#import matplotlib.pyplot as plt
import os
import time
import copy
import warnings
import sys

class MSD:
    
    def runMSD(self, comx, comy, comz, Lx, Ly, Lz, Lx2, Ly2, Lz2, moltype, moltypel, dt, tsjump, output, ver, skip, num_init):
        
        """
        This function calculates the mean square displacement for all molecule
        types in the system from center of mass positions
        """
        
        (comx, comy, comz) = self.unwrap(comx,comy,comz,Lx,Ly,Lz,Lx2,Ly2,Lz2)
        if ver > 0:        
            print('unwrap complete')
        num_timesteps = len(comx)
        (num_init, len_MSD, MSD, diffusivity) = self.gettimesteps(num_timesteps, moltypel,skip, num_init)
        (molcheck,nummol) = self.setmolarray(moltype,moltypel)
        for i in range(skip,num_init+skip):
            for j in range(i,i+len_MSD):
                r2 = self.calcr2(comx, comy, comz, i, j)
                MSD = self.MSDadd(r2, MSD, molcheck, i, j)
            if ver:
                sys.stdout.write('\rMSD calculation {:.2f}% complete'.format((i+1-skip)*100.0/num_init))
        if ver:
            sys.stdout.write('\n')
        MSD = self.MSDnorm(MSD, num_init, nummol)
        Time = self.createtime(dt, tsjump, len_MSD)
        self.append_dict(MSD, moltypel, output, Time)
        return output
    
    def unwrap(self, comx, comy, comz, Lx, Ly, Lz, Lx2, Ly2, Lz2):
        #unwraps the coordintes of the molecules 
        #assumes if a molecule is more than half a box length away from its 
        #previous coordinte that it passed through a periodic boundary
        for i in range(1,len(comx)):
            for j in range(0,len(comx[i])):
                if (comx[i][j]-comx[i-1][j])>Lx2:
                    while (comx[i][j]-comx[i-1][j])>Lx2:
                        comx[i][j] -= Lx
                elif (comx[i][j]-comx[i-1][j]) < (-Lx2):
                    while (comx[i][j]-comx[i-1][j]) < (-Lx2):
                        comx[i][j] += Lx
                    
                if (comy[i][j]-comy[i-1][j])>Ly2:
                    while (comy[i][j]-comy[i-1][j])>Ly2:
                        comy[i][j] -= Ly
                elif (comy[i][j]-comy[i-1][j]) < (-Ly2):
                    while (comy[i][j]-comy[i-1][j]) < (-Ly2):
                        comy[i][j] += Ly
                    
                if (comz[i][j]-comz[i-1][j])>Lz2:
                    while (comz[i][j]-comz[i-1][j])>Lz2:
                        comz[i][j] -= Lz
                elif (comz[i][j]-comz[i-1][j])< (-Lz2):
                    while (comz[i][j]-comz[i-1][j])< (-Lz2):
                        comz[i][j] += Lz         
        return (comx, comy, comz)
        
    def gettimesteps(self, num_timesteps,moltypel,skip,num_init):
        #Calculates the length of the trajectory
        #Uses length to determine length of MSD and number of initial timesteps
        if num_init==None:
            num_init = int(np.floor((num_timesteps-skip)/2))
        else:
            num_init=int(num_init)
        len_MSD = num_timesteps-skip-num_init
        MSD = np.zeros((len(moltypel), len_MSD))
        diffusivity = []
        return (num_init, len_MSD, MSD, diffusivity)
        
    def setmolarray(self, moltype, moltypel):
        #Generates arrays for dot product calculation
        #Array is MxN where M is number of molecule types and N is number of molecules
        #value is 1 if molecule N is of type M else is 0
        molcheck = np.zeros((len(moltypel), len(moltype)))
        for i in range(0,len(moltype)):
            molcheck[moltype[i]][i]=1
        nummol = np.zeros(len(moltypel))
        for i in range(0,len(nummol)):
            nummol[i] = np.sum(molcheck[i])        
        return (molcheck,nummol)
        
    def calcr2(self, comx, comy, comz, i, j):
        #Calculates distance molecule has traveled between steps i and j
        r2 = (comx[j]-comx[i])**2 + (comy[j]-comy[i])**2 + (comz[j]-comz[i])**2
        
        return r2
        
    def MSDadd(self, r2, MSD, molcheck, i, j):
        #Uses dot product to calculate average MSD for a molecule type
        for k in range(0,len(molcheck)):        
            sr2 = np.dot(r2, molcheck[k])
            MSD[k][j-i] += sr2
        return MSD
        
    def MSDnorm(self, MSD,MSDt,nummol):
        #Normalize the MSD by number of molecules and number of initial timesteps
        for i in range(0,len(nummol)):        
            MSD[i] /= MSDt*nummol[i]
        
        return MSD
        
    def createtime(self, dt, tsjump, MSDt):
        #Creates an array of time values
        Time = np.arange(0,MSDt,dtype=float)
        Time *= dt*tsjump
        return Time
    
            
    def append_dict(self, MSD, moltypel, output, Time):
        #Write MSD to output dictionary
        output['MSD'] = {}
        output['MSD']['units']='Angstroms^2, fs'
        for i in range(0,len(moltypel)):
            output['MSD'][moltypel[i]] = copy.deepcopy(MSD[i].tolist())
            
        output['MSD']['time'] = copy.deepcopy(Time.tolist())
