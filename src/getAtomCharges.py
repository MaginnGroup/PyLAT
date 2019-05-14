# -*- coding: utf-8 -*-
"""
Created on Thu May  7 14:10:12 2015

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
class getatomcharges:
    
    def findnumatoms(self, datfilename):
        #returns number of atoms in the simulation
        datfile = open(datfilename)
        foundnumatoms=False
        datfile.readline()
        while foundnumatoms==False:
            line = datfile.readline()
            line = line.split()
            if len(line) >= 2:
                if line[1] == 'atoms':
                    n=int(line[0])
                    foundnumatoms = True
        datfile.close()
        return n
    
    def getmolcharges(self, datfilename, n):
        #Returns arrays with the charge of all atoms and molecules in the system
        datfile = open(datfilename)
        for j in range(0,4):
            datfile.readline()
        atomcharges = np.zeros(n)
        mol = np.zeros(n)
        foundatoms= False
        readingcharges = True
        
        while foundatoms == False:
            line = datfile.readline()
            line = line.split()
            
            if len(line) > 0:
                if line[0] == 'Atoms':
                    foundatoms = True
                    datfile.readline()
                
        while readingcharges == True:
            line = datfile.readline()
            line = line.split()
            if len(line) == 10 or len(line) == 7:
                atomcharges[int(line[0])-1] = float(line[3])
                mol[int(line[0])-1] = int(line[1])
                
            else:
                readingcharges = False
                
        nummol = int(max(mol))
        molcharges = np.zeros(nummol)
        for atom in range(0,n):
            molcharges[int(mol[int(atom)])-1] += atomcharges[int(atom)]
            
        datfile.close()
        return (molcharges,atomcharges,n)
        
    def molchargedict(self, molcharges, moltypel, moltype):
        #creates a dictionary assigning charge to a molecule type
        molcharge = {}
        for molecules in range(0,len(moltypel)):
            molcharge[moltypel[molecules]] = molcharges[moltype.index(molecules)]
        return molcharge
