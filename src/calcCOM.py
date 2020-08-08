# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 10:07:34 2015

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
import os
import src.calccomf as calccomf
import sys


class calcCOM:
    
    def calcCOM(self, trjfilename, datfilename, ver):
        """
        This function will read in the x,y,z positions of all atoms for a 
        timestep and then calculate the center of mass for all molecules in the
        system. Returns the x, y and z coordinates of all molecules for all
        timesteps. 
        
        Center of mass coordinates are often used in place of atomic positions
        to minimize the amount of memory required for calculations
        """
        
        
        (num_lines, n, num_timesteps, count, line)=self.getnum(trjfilename)
        (Lx, Lx2, Ly, Ly2, Lz, Lz2) = self.getdimensions(trjfilename[0])  
        (x,y,z,mol,atype) = self.createarrays(n)
        (xcol, ycol, zcol, molcol, typecol) = self.getcolumns(trjfilename[0])
        atommass = self.getmass(datfilename)
        for i in range(0,len(trjfilename)):
            trjfile = open(trjfilename[i])
            while line[i] < num_lines[i]:
                (x,y,z,mol,atype,line) = self.readdata(trjfile, n, line, x, y, z, mol, atype, xcol, ycol, zcol, molcol, typecol,i)
                if count == 0:
                    (nummol, comx, comy, comz, molmass) = self.comprep(mol, n, atype, atommass, num_timesteps)
                (comx, comy, comz, count) = self.calccom(comx, comy, comz, x, y, z, mol, atype, atommass, molmass, Lx, Ly, Lz, Lx2, Ly2, Lz2, n, count, nummol)
                if ver:                
                    sys.stdout.write('\rCOM calculation {:.2f}% complete'.format(count*100.0/num_timesteps))
            trjfile.close()
        if ver:
            sys.stdout.write('\n')
        return (comx, comy, comz, Lx, Ly, Lz, Lx2, Ly2, Lz2)
        
    def getnum(self,trjfilename):
        # uses the trjectory file and returns the number of lines and the number of atoms
        trjfile = open(trjfilename[0])
        for i in range(0,3):
            trjfile.readline()
        n = int(trjfile.readline())
        trjfile.close()
        num_timesteps=1
        num_lines=[]
        for i in range(0,len(trjfilename)):
            num_lines.append(int(sum(1 for line in open(trjfilename[i]))))
            num_timesteps += int(num_lines[i] / (n+9))-1
        line = [10 for x in trjfilename]
        for j in range(1,len(trjfilename)):
            line[j] += n+9
        count = 0
        return (num_lines, n, num_timesteps, count, line)
        
    def getdimensions(self,trjfilename):
        # uses trjectory file to get the length of box sides
        trjfile = open(trjfilename)
        for i in range(0,5):
            trjfile.readline()
        xbounds = trjfile.readline()
        xbounds = xbounds.split()
        ybounds = trjfile.readline()
        ybounds = ybounds.split()
        zbounds = trjfile.readline()
        zbounds = zbounds.split()
        Lx = float(xbounds[1])-float(xbounds[0])
        Lx2 = Lx/2
        Ly = float(ybounds[1])-float(ybounds[0])
        Ly2 = Ly/2
        Lz = float(zbounds[1])-float(zbounds[0])
        Lz2 = Lz/2 
        trjfile.close()
        return (Lx, Lx2, Ly, Ly2, Lz, Lz2)
        
    def createarrays(self,n):
        #creates numpy arrays for data reading
        x = np.zeros(n)
        y = np.zeros(n)
        z = np.zeros(n)
        mol = np.zeros(n)
        atype = np.zeros(n)
        return (x,y,z,mol,atype)

    def getcolumns(self,trjfilename):
        # defines the columns each data type is in in the trjectory file
        trjfile = open(trjfilename)
        for j in range(0,8):
            trjfile.readline()
        inline = trjfile.readline()
        inline = inline.split()
        inline.remove('ITEM:')
        inline.remove('ATOMS')
        xcol = inline.index('x')
        ycol = inline.index('y')
        zcol = inline.index('z')
        molcol = inline.index('mol')
        typecol = inline.index('type')
        trjfile.close()
        return (xcol, ycol, zcol, molcol, typecol)
        
    def getmass(self, datfilename):
        # returns a dictionary of the mass of each atom type
        atommass = {}
        foundmass= False
        readingmasses = True
        atomnum = 1
        datfile = open(datfilename)
        for i in range(0,4):
            datfile.readline()
        
        while foundmass == False:
            line = datfile.readline()
            line = line.split()
            
            if len(line) > 0:
                if line[0] == 'Masses':
                    foundmass = True
                    datfile.readline()
                
        while readingmasses == True:
            line = datfile.readline()
            line = line.split()
            if len(line) > 0:
                if int(line[0]) == atomnum:
                    atommass[int(line[0])] = float(line[1])
                    atomnum += 1
                    
                else:
                    readingmasses = False
                
            else:
                readingmasses = False
        datfile.close()
        return atommass
        
    def readdata(self, trjfile, n, line, x, y, z, mol, atype, xcol, ycol, zcol, molcol, typecol, i):
        # reads data from trjectory file into precreated arrays
        for j in range(0,9):
            trjfile.readline()
        for a in range(0,n):
            inline = trjfile.readline()
            inline = inline.split()
            x[a]= inline[xcol]
            y[a]= inline[ycol]
            z[a]= inline[zcol]
            mol[a]= inline[molcol]
            atype[a]= inline[typecol]
            
        line[i] += n+9
        return (x,y,z,mol,atype,line)
        
    def comprep(self, mol, n, atype, atommass, num_timesteps):
        #creates arrays to prepare for center of mass calculations
        nummol = int(max(mol))
        comx = [[0 for x in range(nummol)]for x in range(num_timesteps)]
        comy = [[0 for x in range(nummol)]for x in range(num_timesteps)]
        comz = [[0 for x in range(nummol)]for x in range(num_timesteps)]

        molmass = np.zeros(nummol)
        for atom in range(0,n):
            molmass[int(mol[atom]-1)] += atommass[atype[atom]]
            
        return (nummol, comx, comy, comz, molmass)
        
    def calccom(self, comx, comy, comz, x, y, z, mol, atype, atommass, molmass, Lx, Ly, Lz, Lx2, Ly2, Lz2, n, count, nummol):
        #calculates the center of mass for each molecule
        amass = np.zeros(n)
        for i in range(0,n):
            amass[i] = atommass[atype[i]]
            
            
        #Calls a fortran code to increase the efficiency of the calculations    
        (comxt, comyt, comzt)= calccomf.calccom(n, nummol, x, y, z, mol, amass, molmass, Lx, Ly, Lz, Lx2, Ly2, Lz2)
        comx[count] += comxt
        comy[count] += comyt
        comz[count] += comzt
        count += 1
            
        return (comx, comy, comz, count)
