# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 10:41:46 2015

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
