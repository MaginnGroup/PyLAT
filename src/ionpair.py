# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 10:49:36 2015

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
import src.calcdistances as calcdistances
import numpy as np
from scipy.integrate import cumtrapz
import src.ipcorr as ipcorr
from scipy.optimize import curve_fit
import copy
import sys
import warnings

class ionpair:

    
    def runionpair(self,comx,comy,comz,Lx,Ly,Lz,moltypel,moltype,tsjump,dt,output,ver,skipframes):
        
        """
        Calculates the ion pair lifetime for all combinations of molecule types
        in the system. An ion pair for types A and B is defined as the closest 
        molecule of type B around a molecule of type A
        
        The integral is fit using multiexponentials from one to five 
        exponentials to obtain a good fit without overfitting
        """
        
        output['Ion_Pair_Lifetime'] = {}
        output['Ion_Pair_Lifetime']['Units'] = 'picoseconds'
        output['Ion_Pair_Lifetime']['Explanation'] = 'The Ion Pair Lifetime correlation function is fit to a single exponential, a double exponential up to 5 exponentials. The results shown are the result of these successive fittings'
        (closest,begin,end,C)=self.init(len(comx[0]),moltypel,len(comx),moltype)
        for step in range(0,len(comx)):
            r = self.calcdistance(comx[step],comy[step],comz[step],Lx,Ly,Lz)
            self.findclosest(r,closest,begin,end,step)
            if ver:
                sys.stdout.write('\rIPL distance calculation {:.2f}% complete'.format((step+1)*100.0/len(comx)))
        
        if ver:
            sys.stdout.write('\n')
        
        correlation = self.correlation(closest,moltype,moltypel,ver,skipframes)
        if ver:
            print('correlation complete')
        time = []
        for i in range(0,len(correlation)):
            time.append(float(i*tsjump*dt/1000))
        begin = int(1000/dt/tsjump)
        end = len(time)
        for i in range(0,len(moltypel)):
            for j in range(0,len(moltypel)):
                if i != j:
                    y = []
                    end = 0
                    for k in range(0,len(correlation)):
                        y.append(float(correlation[k][i][j]))
                        if correlation[k][i][j] <= 0.04 and end ==0:
                            end = k
                    if end==0:
                        end = len(y)
                    (IPL,r2) = self.curvefit(y,time,begin,end)
                    output['Ion_Pair_Lifetime']['{0} around {1}'.format(moltypel[j],moltypel[i])]=IPL
                    output['Ion_Pair_Lifetime']['{0} around {1} r2'.format(moltypel[j],moltypel[i])]=r2
                    output['Ion_Pair_Lifetime']['{0} around {1} correlation'.format(moltypel[j],moltypel[i])]=copy.deepcopy(y)                   
        output['Ion_Pair_Lifetime']['Correlation_Time']=copy.deepcopy(time)
    
    def init(self,nummol, moltypel, numtimesteps,moltype):
        #initializes arrays for the calculations
        closest = np.zeros((numtimesteps,nummol,len(moltypel)))
        C = np.zeros((len(moltypel),len(moltypel)))
        begin = [0]
        end = []
        for i in range(1,len(moltype)):
            if moltype[i]!= moltype[i-1]:
                end.append(i)
                begin.append(i)
                
        end.append(len(moltype))
        return (closest,begin,end,C)
        
    def calcdistance(self,comx,comy,comz,Lx,Ly,Lz):
        #Runs a fortran script calculating the distance between all molecules
        r = calcdistances.calcdistances(len(comx),comx,comy,comz,Lx,Ly,Lz)
        return r
        
    def findclosest(self,r,closest,begin,end,timestep):
        #Search molecules to find the closest molecules at each timestep
        for i in range(0,len(r)):
            for j in range(0,len(begin)):
                distance = 10000
                for k in range(begin[j],end[j]):
                    if r[i][k]<distance:
                        distance = r[i][k]
                        closest[timestep][i][j] = k
    
    def correlation(self,closest,moltype,moltypel,ver,skipframes):
        #Runs a fortran script perfroming the correlation function
        correlation = ipcorr.ipcorr(closest,skipframes,len(closest),len(closest[0]),len(closest[0][0]),(len(closest)-skipframes)/2,moltype)
        return correlation
    
    def exponential1(self,x,A1,B1):
        return A1*np.exp(-x/B1)
        
    def exponential2(self,x,A1,A2,B1,B2):
        return A1*np.exp(-x/B1) + A2*np.exp(-x/B2)
        
    def exponential3(self,x,A1,A2,A3,B1,B2,B3):
        return A1*np.exp(-x/B1) + A2*np.exp(-x/B2) + A3*np.exp(-x/B3)
        
    def exponential4(self,x,A1,A2,A3,A4,B1,B2,B3,B4):
        return A1*np.exp(-x/B1) + A2*np.exp(-x/B2) + A3*np.exp(-x/B3) + A4*np.exp(-x/B4)
        
    def exponential5(self,x,A1,A2,A3,A4,A5,B1,B2,B3,B4,B5):
        return A1*np.exp(-x/B1) + A2*np.exp(-x/B2) + A3*np.exp(-x/B3) + A4*np.exp(-x/B4) + A5*np.exp(-x/B5)
        
    def exponential6(self,x,A1,A2,A3,A4,A5,A6,B1,B2,B3,B4,B5,B6):
        return A1*np.exp(-x/B1) + A2*np.exp(-x/B2) + A3*np.exp(-x/B3) + A4*np.exp(-x/B4) + A5*np.exp(-x/B5) + A6*np.exp(-x/B6)
        
    def curvefit(self,correlation,time, begin,end):
        #Fit the exponential functions to the correlation function to estimate the ion pair lifetime 
        funlist = [self.exponential1,self.exponential2,self.exponential3,self.exponential4,self.exponential5]
        IPL = []
        r2 = []
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for fun in funlist:
                popt, pcov = curve_fit(fun, time[begin:end], correlation[begin:end],maxfev=100000000)
                fit = []
                for i in time:
                    fit.append(fun(i,*popt))
                yave = np.average(correlation[begin:end])
                SStot = 0
                SSres = 0
                for l in range(begin,end):
                    SStot += (correlation[l]-yave)**2
                    SSres += (correlation[l]-fit[l])**2
                r2.append(1-SSres/SStot)
                IPL.append(0)
                for i in range(0,int(len(popt)/2)):
                    IPL[-1] += popt[i]*popt[i+int(len(popt)/2)]
        return (IPL,r2)
