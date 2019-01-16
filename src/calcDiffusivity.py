# -*- coding: utf-8 -*-
"""
Created on Wed May 27 16:05:33 2015

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
from scipy import stats
import warnings

class calcdiffusivity:
        
    
    def calcdiffusivity(self, output, moltypel, dt, tol):
        
        """
        This function fits the mean square displacement to calculate the
        diffusivity for all molecule types in the system
        """
        
        output['Diffusivity'] = {}
        output['Diffusivity']['units'] = 'm^2/s'
        for i in range(0,len(moltypel)):
            MSD = output['MSD'][moltypel[i]]
            lnMSD = np.log(MSD[1:])
            time = output['MSD']['time']
            lntime = np.log(time[1:])
            firststep = self.findlinearregion(lnMSD, lntime, dt, tol)
            #self.writeLogLog(lnMSD,lntime,moltypel[i])
            diffusivity = self.getdiffusivity(time, MSD, firststep)
            output['Diffusivity'][moltypel[i]] = diffusivity
            
    
    
    def findlinearregion(self, lnMSD, lntime, dt, tol):
        #Uses the slope of the log-log plot to find linear regoin of MSD
        timestepskip=np.ceil(500/dt)
        linearregion=True
        maxtime = len(lnMSD)
        numskip=1
        while linearregion == True:
            if numskip*timestepskip+1 > maxtime:
                firststep = maxtime-1-(numskip-1)*timestepskip
                return firststep
                linearregion= False
            else:
                t1=int(maxtime-1-(numskip-1)*timestepskip)
                t2=int(maxtime-1-numskip*timestepskip)
                slope = (lnMSD[t1]-lnMSD[t2])/(lntime[t1]-lntime[t2])
                if abs(slope-1.) < tol:
                    numskip += 1
                else:
                    firststep=t1
                    return firststep
                    linearregion = False
                
                
    def getdiffusivity(self, Time, MSD, firststep):
        #Fits the linear region of the MSD to obtain the diffusivity
        calctime = []        
        calcMSD = []
        for i in range(int(firststep), len(Time)):
            calctime.append(Time[i])
            calcMSD.append(MSD[i])
        if len(calctime)==1:
            diffusivity = 'runtime not long enough'
        else:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                line = stats.linregress(calctime,calcMSD)
            slope = line[0]
            diffusivity=slope/600000
        return diffusivity
    
    def writeLogLog(self,lnMSD,lntime,moltype):
        outfile = open('LogLog{}.dat'.format(moltype),'w')
        for i in range(0,len(lnMSD)):
            outfile.write('{}\t{}\n'.format(lntime[i],lnMSD[i]))
        outfile.close()
            
