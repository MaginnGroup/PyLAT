# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 14:09:31 2017

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
import __future__
import src.calccomf as calccomf

class distSearch:
    def distSearch(self, mol1, mol2, dist, deltaDist, firstFrame, numSamples, trjfilename, datfilename, moltype, moltypel,output):
            
        
        """
        This function will search through a trajectory to find examples where two 
        molecules of given types are at a certain distance apart. Useful for 
        generation of images
        """
        
        output['Distance_Search'] = {}      
        (num_lines, n, num_timesteps, count, line, numFound, frame)=self.getnum(trjfilename)
        (Lx, Lx2, Ly, Ly2, Lz, Lz2) = self.getdimensions(trjfilename[0])  
        (x,y,z,mol,atype,aid) = self.createarrays(n)
        (xcol, ycol, zcol, molcol, typecol, idcol) = self.getcolumns(trjfilename[0])
        atommass = self.getmass(datfilename)
        for i in range(0,len(trjfilename)):
            trjfile = open(trjfilename[i])
            for j in range(0,firstFrame-1):
                for k in range(0,n+9):
                    trjfile.readline()
                line[i] += n+9
                #print 'frame {}'.format(frame)
                frame+=1
            while line[i] < num_lines[i] and numFound < numSamples:
                #print 'frame {}'.format(frame)
                (x,y,z,mol,atype,line,aid) = self.readdata(trjfile, n, line, x, y, z, mol, atype, aid, xcol, ycol, zcol, molcol, typecol,i,idcol)
                if count == 0:
                    (nummol, comx, comy, comz, molmass) = self.comprep(mol, n, atype, atommass, num_timesteps)
                    molid = self.molID(mol,aid,moltype)
                (comx, comy, comz, count) = self.calccom(comx, comy, comz, x, y, z, mol, atype, atommass, molmass, Lx, Ly, Lz, Lx2, Ly2, Lz2, n, count, nummol)
                (numFound,output) = self.distCalc(comx, comy, comz, mol1, mol2, dist, deltaDist, moltype, molid, numFound, numSamples, Lx, Ly, Lz,frame, moltypel, output)
                frame += 1
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
        numFound = 0
        frame = 1
        return (num_lines, n, num_timesteps, count, line, numFound, frame)
        
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
        aid = np.zeros(n)
        return (x,y,z,mol,atype,aid)

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
        idcol = inline.index('id')
        trjfile.close()
        return (xcol, ycol, zcol, molcol, typecol, idcol)
        
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
        
    def readdata(self, trjfile, n, line, x, y, z, mol, atype, aid, xcol, ycol, zcol, molcol, typecol, i, idcol):
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
            aid[a] = inline[idcol]
            
        line[i] += n+9
        return (x,y,z,mol,atype,line,aid)
        
    def comprep(self, mol, n, atype, atommass, num_timesteps):
        #creates arrays to prepare for center of mass calculations
        nummol = int(max(mol))
        comx = [0 for x in range(nummol)]
        comy = [0 for x in range(nummol)]
        comz = [0 for x in range(nummol)]

        molmass = np.zeros(nummol)
        for atom in range(0,n):
            molmass[int(mol[atom]-1)] += atommass[atype[atom]]
            
        return (nummol, comx, comy, comz, molmass)
        
    def molID(self, mol, aid, moltype):
        #generates arrays for ids of example molecules
        molid = [[] for i in range(0,len(moltype))]
        for i in range(0,len(mol)):
            molid[int(mol[i])-1].append(int(aid[i]))
        return molid
        
    def calccom(self, comx, comy, comz, x, y, z, mol, atype, atommass, molmass, Lx, Ly, Lz, Lx2, Ly2, Lz2, n, count, nummol):
        #calculates the center of mass for each molecule
        amass = np.zeros(n)
        for i in range(0,n):
            amass[i] = atommass[atype[i]]
            
        (comxt, comyt, comzt)= calccomf.calccom(n, nummol, x, y, z, mol, amass, molmass, Lx, Ly, Lz, Lx2, Ly2, Lz2)
        comx = np.array(comxt)
        comy = np.array(comyt)
        comz = np.array(comzt)
    
        count += 1
            
        return (comx, comy, comz, count)
    def distCalc(self, comx, comy, comz, mol1, mol2, dist, deltaDist, moltype, molid, numFound, numSamples, Lx, Ly, Lz, frame, moltypel, output):
        #Finds molecules at the given distance apart
        ind = []       
        for i in range(0,len(moltype)):
            ind.append(moltype[i]==moltypel.index(mol2))
        indid = [l for l, x in enumerate(ind) if x]
        for i in range(0,len(moltype)):
            if moltype[i]==moltypel.index(mol1):
                dx = comx[ind]-np.array([comx[i] for j in range(0,len(indid))])
                dy = comy[ind]-np.array([comy[i] for j in range(0,len(indid))])
                dz = comz[ind]-np.array([comz[i] for j in range(0,len(indid))])
                
                dx -= Lx * np.around(dx / Lx)
                dy -= Ly * np.around(dy / Ly)
                dz -= Lz * np.around(dz / Lz)
                
                r2 = dx ** 2 + dy ** 2 + dz ** 2
                r = np.sqrt(r2)
                
                for j in range(0,len(r)):
                    if (np.abs(r[j]-dist) < deltaDist) and (numFound < numSamples):
                        if ((mol1 == mol2) and (i < indid[j])) or (not mol1==mol2):
                            numFound += 1
                            molid[i].sort()
                            molid[indid[j]].sort()
                            print('Sample {} Found'.format(numFound))
                            output['Distance_Search']['Sample_{}'.format(numFound)] = {}
                            output['Distance_Search']['Sample_{}'.format(numFound)]['Distance'] = float(r[j])
                            output['Distance_Search']['Sample_{}'.format(numFound)]['Frame'] = int(frame)
                            output['Distance_Search']['Sample_{}'.format(numFound)]['Molecule_1_IDs'] = molid[i]
                            output['Distance_Search']['Sample_{}'.format(numFound)]['molecule_2_IDs'] = molid[indid[j]]
                            
        return (numFound,output)
        
        
        
