# -*- coding: utf-8 -*-
"""
Created on Wed May 18 16:05:34 2016

@author: mhumbert
"""

import numpy as np
from scipy import optimize

class fitVisc:

    def doubexp(self,x,A,alpha,tau1, tau2):
        return A*alpha*tau1*(1-np.exp(-x/tau1))+A*(1-alpha)*tau2*(1-np.exp(-x/tau2))
    
    def fitvisc(self,output):
        time = output['Viscosity']['Time']
        visc = output['Viscosity']['Average Integral']
        stddev = output['Viscosity']['Standard Deviation']
        popt2=[1e-3,1.5e-1,3e4,1e3]
        
        foundcutoff = False
        foundstart = False
        start = 1
        while not foundstart and start<len(visc):
            if time[start] > 3000:
                foundstart = True
            else:
                start+=1
        cut = 1
        while not foundcutoff and cut<len(visc):
            if stddev[cut] > 0.4*visc[cut]:
                foundcutoff = True
            else:
                cut += 1
        #cut = len(visc)
        popt2,pcov2 = optimize.curve_fit(self.doubexp, time[start:cut], visc[start:cut],maxfev=1000000,p0=popt2, sigma=stddev[start:cut])
        
        fit = []
        for t in time:
            fit.append(self.doubexp(t,*popt2))
            
        output['Viscosity']['Fit']=fit
        output['Viscosity']['Value']=popt2[0]*popt2[1]*popt2[2]+popt2[0]*(1-popt2[1])*popt2[3]
        
        return(output)
