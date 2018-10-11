# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 10:41:46 2015

@author: mhumbert
"""
import sys

class gettimedata:
    
    def getdt(self, logfilename):
        #determines the timestep of the simulation
        numlines=int(sum(1 for line in open(logfilename)))
        logfile = open(logfilename)
        foundtimestep=False
        count=0
        while foundtimestep==False:
            if count > numlines:
                sys.exit('could not find timestep size in log file')
            inline = logfile.readline()
            inline = inline.split()
            if len(inline) > 0:
                if inline[0] == 'timestep':
                    dt = float(inline[1])
                    foundtimestep=True
            count+=1
        logfile.close()
        
        return dt
        
    def getjump(self, trjfilename):
        #determines the frequency of output for the trajectory file
        trjfile = open(trjfilename)
        trjfile.readline()
        t1 = trjfile.readline()
        t1 = int(t1)
        trjfile.readline()
        n = int(trjfile.readline())
        for i in range(0,n+6):
            trjfile.readline()
        t2 = int(trjfile.readline())
        tsjump = t2-t1
        return tsjump