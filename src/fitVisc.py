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
    
    def doubexp1(self,x,A,alpha,tau1, tau2):
        return A*alpha*tau1*(1-np.exp(-x/tau1))
    
    def doubexp2(self,x,A,alpha,tau1, tau2):
        return A*(1-alpha)*tau2*(1-np.exp(-x/tau2))
    
    def singleexp(self, x, A, tau):
        return A*(1-np.exp(-x/tau))
    
    def fitvisc(self,time,visc,stddev,plot):
        popt2=[1e-3,1.5e-1,3e4,1e3]
        #popt2=[2e-3,5e-2,2e3,2e2]
        #popt2=[1e-4,1e2]
        foundcutoff = False
        foundstart = False
        start = 1
        while not foundstart and start<len(visc):
            if time[start] > 2000:
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
        #popt2,pcov2 = optimize.curve_fit(self.doubexp, time[start:cut], visc[start:cut],maxfev=1000000,p0=popt2, sigma=stddev[start:cut])
        popt2,pcov2 = optimize.curve_fit(self.doubexp, time[start:cut], visc[start:cut],maxfev=1000000,p0=popt2, sigma=stddev[start:cut],bounds=(0,[np.inf,1,np.inf,np.inf]))
        
        fit = []
        fit1 = []
        fit2 = []
        for t in time:
            fit.append(self.doubexp(t,*popt2))
            fit1.append(self.doubexp1(t,*popt2))
            fit2.append(self.doubexp2(t,*popt2))
        Value = popt2[0]*popt2[1]*popt2[2]+popt2[0]*(1-popt2[1])*popt2[3]
        #Value = popt2[0]
        
        if plot:
            from matplotlib import pyplot as plt
            from matplotlib import rcParams
            rcParams.update({'font.size':14})
            print('Viscosity estimate is {}'.format(Value))
            print('A={}, alpha={}, tau1={}, tau2={}'.format(popt2[0],popt2[1],popt2[2],popt2[3]))
            print('Time cutoff is {}'.format(time[cut]))
            plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
            plt.plot(time[:len(visc)],visc,label='Viscosity')
            plt.plot(time[:len(fit)],fit,label='Double Exponential fit')
            plt.plot(time[:len(fit1)],fit1,label='Contribution of tau1')
            plt.plot(time[:len(fit2)],fit2,label='Contribution of tau2')
            plt.axvline(time[cut])
            plt.ylabel('Viscosity(cP)')
            plt.xlabel('Time (fs)')
            plt.legend()
            plt.show()
        
        return(Value)
