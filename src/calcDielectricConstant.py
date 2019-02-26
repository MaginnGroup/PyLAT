# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 14:49:27 2016

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
import sys

class calcDielectricConstant:
    def calcDEC(self, atomcharges, trjfilename, T, output, V, ver,start):
        
        """
        This function calculates the dielectric constant of the system using 
        the fluctuations in the dipole moment
        """
        
        (num_lines, n, num_timesteps, count, line) = self.getnum(trjfilename)
        (Lx, Lx2, Ly, Ly2, Lz, Lz2) = self.getdimensions(trjfilename)
        (x,y,z,mol,Mx,My,Mz) = self.createarrays(n,num_timesteps)
        (xcol, ycol, zcol, molcol, idcol) = self.getcolumns(trjfilename)
        if ver >= 1:
            print('beginning dipole moment calculation')
        for i in range(0,len(trjfilename)):
            trjfile = open(trjfilename[i])
            while line[i] < num_lines[i]:
                (x,y,z,mol,line) = self.readdata(trjfile, n, line, x, y, z, mol, xcol, ycol, zcol, molcol, idcol, i)
                self.unwrap(x, y, z, mol, Lx, Lx2, Ly, Ly2, Lz, Lz2)
                (Mx[count],My[count],Mz[count]) = self.calcdipolemoment(x, y, z, atomcharges)
                count += 1
                if ver > 1:
                    sys.stdout.write('\rdipole moment calculation {:.2f}% complete'.format(count*100.0/num_timesteps))
        if ver > 1:
            sys.stdout.write('\n')
        if ver >= 1:
            print('calculating dielectric constant')
        (AveM2, AveM) = self.calcaverage(Mx, My, Mz,start)
        dielectric = self.calcdielectric(AveM2, AveM, V, T)
        output['Dielectric Constant'] = {}
        output['Dielectric Constant']['units'] = 'unitless'
        output['Dielectric Constant']['explanation'] = 'returns a list of the cumulative dielectric constant for the number of timesteps used in the averages to show convergence.'
        output['Dielectric Constant']['list'] = dielectric
        return output
                
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
        trjfile = open(trjfilename[0])
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
        
    def createarrays(self,n, num_timesteps):
        #creates numpy arrays for data reading
        x = np.zeros(n)
        y = np.zeros(n)
        z = np.zeros(n)
        mol = np.zeros(n)
        Mx = np.zeros(num_timesteps)
        My = np.zeros(num_timesteps)
        Mz = np.zeros(num_timesteps)
        return (x,y,z,mol,Mx,My,Mz)
        
    def getcolumns(self,trjfilename):
        # defines the columns each data type is in in the trjectory file
        trjfile = open(trjfilename[0])
        for j in range(0,8):
            trjfile.readline()
        inline = trjfile.readline()
        inline = inline.split()
        inline.remove('ITEM:')
        inline.remove('ATOMS')
        try:
            xcol = inline.index('x')
        except:
            xcol = inline.index('xu')
        try:
            ycol = inline.index('y')
        except:
            ycol = inline.index('yu')
        try:
            zcol = inline.index('z')
        except:
            zcol = inline.index('zu')
        molcol = inline.index('mol')
        idcol = inline.index('id')
        trjfile.close()
        return (xcol, ycol, zcol, molcol, idcol)
        
    def readdata(self, trjfile, n, line, x, y, z, mol, xcol, ycol, zcol, molcol, idcol, i):
        # reads data from trjectory file into precreated arrays
        for j in range(0,9):
            trjfile.readline()
        for a in range(0,n):
            inline = trjfile.readline()
            inline = inline.split()
            aid = inline[idcol]
            x[int(aid)-1]= inline[xcol]
            y[int(aid)-1]= inline[ycol]
            z[int(aid)-1]= inline[zcol]
            mol[int(aid)-1]= inline[molcol]
            
        line[i] += n+9
        return (x,y,z,mol,line)
        
    def unwrap(self, x, y, z, mol, Lx, Lx2, Ly, Ly2, Lz, Lz2):
        #Ensures all atoms in a molecule are in the same image
        for atom in range(1,len(x)):
            if mol[atom] == mol[atom-1]:
                if x[atom] - x[atom-1] > Lx2:
                    x[atom] -= Lx
                elif x[atom] - x[atom-1] < -1*Lx2:
                    x[atom] += Lx
                if y[atom] - y[atom-1] > Ly2:
                    y[atom] -= Ly
                elif y[atom] - y[atom-1] < -1*Ly2:
                    y[atom] += Ly
                if z[atom] - z[atom-1] > Lz2:
                    z[atom] -= Lz
                elif z[atom] - z[atom-1] < -1*Lz2:
                    z[atom] += Lz
                    
    def calcdipolemoment(self, x, y, z, atomcharges):
        #calculates the dipole moment for a given timestep
        mx = np.dot(x,atomcharges)
        my = np.dot(y,atomcharges)
        mz = np.dot(z,atomcharges)
        return (mx,my,mz)
        
    def calcaverage(self,Mx,My,Mz,start):
        #Calculates the cumulative average dipole moment and dipole moment squared for fluctuation calculations
        normarray = np.arange(1,len(Mx)-start+1)
        AveM2 = np.cumsum(Mx[start:]**2+My[start:]**2+Mz[start:]**2)/normarray
        AveM = (np.cumsum(Mx[start:])**2+np.cumsum(My[start:])**2+np.cumsum(Mz[start:])**2)/normarray**2
        return (AveM2, AveM)
        
    def calcdielectric(self, AveM2, AveM, V, T):
        #Calculates the running average of the dielectric constant
        kb = 1.38065e-23
        epsilon0 = 8.85418e-12
        AveM2 *= 1e-20*1.6022e-19**2 #correct units
        AveM *= 1e-20*1.6022e-19**2
        V *= 1e-30
        dielectric = (1+1/3./epsilon0/V/kb/T*(AveM2-AveM)).tolist()
        return dielectric
    
    def getvolume(self, trjfilename):
        #Calculates volume of the system
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
        Ly = float(ybounds[1])-float(ybounds[0])
        Lz = float(zbounds[1])-float(zbounds[0])
        trjfile.close()
        V = Lx*Ly*Lz
        return (V, Lx, Ly, Lz)

