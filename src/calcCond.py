# -*- coding: utf-8 -*-
"""
Created on Fri May  8 10:17:55 2015

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
from scipy.optimize import curve_fit
from scipy.integrate import cumtrapz
from src.getTimeData import gettimedata
import src.calccomf as calccomf
import copy
import sys

class calcCond:
    g = gettimedata()
    def calcConductivity(self, molcharges, trjfilename,logfilename, datfilename, T, output, moltype, moltypel,ver,firstpoint,tol,Jout):
        
        """
        This function calculates the ionic conductivity of the system using the 
        Green-Kubo formalism
        
        returns both the total conductivity as well as the contribution of each
        species
        
        """
        
        dt = self.g.getdt(logfilename)
        tsjump = self.g.getjump(trjfilename[0])
        (num_lines, n, num_timesteps, count, line) = self.getnum(trjfilename)
        atommass = self.getmass(datfilename)
        V = self.getdimensions(trjfilename[0])
        (vx, vy, vz, mol, atype) = self.createarrays(n, num_timesteps, moltype)
        (vxcol, vycol, vzcol, idcol, molcol, typecol) = self.getcolumns(trjfilename[0])
        dotlist = self.getchargearrays(molcharges,moltype)
        if ver >= 1:
            print('beginning COM velocity calculation')
        for i in range(0,len(trjfilename)):
            trjfile = open(trjfilename[i])
            while count < num_timesteps:
                (vx,vy,vz,line,mol, atype) = self.readdata(trjfile, n, line, vx, vy, vz, vxcol, vycol, vzcol, idcol, i, mol, molcol, atype, typecol)
                if count == 0:
                    (comvx, comvy, comvz, nummol, molmass, jx, jy, jz) = self.COMprep(mol, atype, atommass, n, num_timesteps, moltypel)
                (comvx, comvy, comvz) = self.calcCOMv(comvx, comvy, comvz, vx, vy, vz, mol, atype, atommass, molmass, n, nummol)
                (jx, jy, jz, count) = self.calcj(dotlist, comvx, comvy, comvz, jx, jy, jz, count)
                if ver == 2:
                    sys.stdout.write('\rCOM velocity calculation {:.2f}% complete'.format(count*100.0/num_timesteps))
        if ver == 2:
            sys.stdout.write('\n')
        if ver >= 2:
            print('begining GK conductivity correlation')
        J = self.clacJ(jx, jy, jz, dt, tsjump, firstpoint,ver)
        if Jout:
            self.writeJ(J,tsjump,dt)
        integral = self.integrateJ(J, tsjump*dt)
        (begcon, endcon) = self.findconvergance(J[0],tol)
        time = []
        for i in range(0,len(integral[0])):
            time.append(i*tsjump*dt)
        print(time[begcon])
        print(time[endcon])
        cond = np.zeros(len(J))
        for i in range(0,len(J)):
            ave = self.fitcurve(time, integral[i], begcon, endcon)    
            cond[i] = self.greenkubo(ave, T, V)
        GKintegral = self.greenkubo(integral,T,V)
        fit = []
        for i in range(0,len(time)):
            fit.append(ave)
        output['Conductivity']['Green_Kubo'] = cond[0]
        a = GKintegral[0]
        a = a.tolist()
        output['Conductivity']['GK_Integral'] = a
        output['Conductivity']['Time'] = time
        for i in range(1,len(cond)):
            output['Conductivity']['Green_Kubo_{0}'.format(moltypel[i-1])] = cond[i]
            a = GKintegral[i]
            a = a.tolist()
            output['Conductivity']['GK_Integral_{0}'.format(moltypel[i-1])] = a
        return output

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
        Ly = float(ybounds[1])-float(ybounds[0])
        Lz = float(zbounds[1])-float(zbounds[0])
        V = Lx * Ly * Lz /10**30
        trjfile.close()
        return V
    
    def getnum(self,trjfilename):
        # uses the trjectory file and returns the number of lines and the number of atoms
        trjfile = open(trjfilename[0])
        for j in range(0,3):
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
        
    def createarrays(self,n, num_timesteps, moltype):
        #creates numpy arrays for data reading
        vx = np.zeros(n)
        vy = np.zeros(n)
        vz = np.zeros(n)
        mol = np.zeros(n, dtype=int)
        atype = np.zeros(n)
        return (vx, vy, vz, mol, atype)
        
    def getcolumns(self,trjfilename):
        # defines the columns each data type is in in the trjectory file
        trjfile = open(trjfilename)
        for j in range(0,8):
            trjfile.readline()
        inline = trjfile.readline()
        inline = inline.split()
        inline.remove('ITEM:')
        inline.remove('ATOMS')
        vxcol = inline.index('vx')
        vycol = inline.index('vy')
        vzcol = inline.index('vz')
        idcol = inline.index('id')
        molcol = inline.index('mol')
        typecol = inline.index('type')
        trjfile.close()
        return (vxcol, vycol, vzcol, idcol, molcol, typecol)
        
    def readdata(self, trjfile, n, line, vx, vy, vz, vxcol, vycol, vzcol, idcol, i, mol, molcol, atype, typecol):
        # reads data from trjectory file into precreated arrays
        for j in range(0,9):
            trjfile.readline()
        for a in range(0,n):
            inline = trjfile.readline()
            inline = inline.split()
            vx[int(inline[idcol])-1]= float(inline[vxcol])
            vy[int(inline[idcol])-1]= float(inline[vycol])
            vz[int(inline[idcol])-1]= float(inline[vzcol])
            mol[int(inline[idcol])-1] = int(inline[molcol])
            atype[int(inline[idcol])-1] = int(inline[typecol])     
            
        line[i] += n+9
        return (vx,vy,vz,line, mol, atype)
        
    def calcj(self, dotlist, comvx, comvy, comvz, jx, jy, jz, count):
        #calculates the charge flux for a timestep
        #seperated into the contributions by different molecule types
        for i in range(0,len(dotlist)):        
            jx[i][count] = np.dot(dotlist[i], comvx)
            jy[i][count] = np.dot(dotlist[i], comvy)
            jz[i][count] = np.dot(dotlist[i], comvz)
        count += 1
        return (jx, jy, jz, count)
        
    def clacJ(self, jx, jy, jz, dt, tsjump, firstpoint,ver):
        #Calculates the charge flux correlation function for all timesteps
        #J[0] is the total charge correlation function
        #J[m] is the charge correlation function for the mth species
        counter = 0
        J = np.zeros((len(jx)+1,len(jx[0])))
        for l in range(0,len(jx)):
            for m in range(0, len(jx)):
                Jtest = self.correlate(jx[l],jx[m])
                J[0] += Jtest
                J[m+1] += Jtest
                Jtest = self.correlate(jy[l],jy[m])
                J[0] += Jtest
                J[m+1] += Jtest
                Jtest = self.correlate(jz[l],jz[m])
                J[0] += Jtest
                J[m+1] += Jtest
                counter += 1
                if ver == 2:
                    sys.stdout.write('\rGK conductivity correlation {0}% complete'.format(100.*float(counter)/len(jx)**2))
        if ver == 2:
            sys.stdout.write('\n')
        return J

    def integrateJ(self, J, dt):
        #Integrates the charge flux correlation function to calculate conductivity
        integral = np.zeros((len(J),len(J[0])))
        for i in range(0,len(J)):
            integral[i][1:] = cumtrapz(J[i], dx=dt)
        return integral

    def fitcurve(self, time, integral, begin, end):
        #calculates average value of the integral over the range begin:end
        ave = (np.average(integral[begin:end]))
        return ave

    def greenkubo(self, ave, T, V):
        #normalizes the average value of the integral to obtain the conductivity
        k = 1.38e-23
        el = 1.60217e-19
        cond = ave/3/k/T/V*el**2/10**5
        return cond


    def COMprep(self, mol, atype, atommass, n, num_timesteps, moltypel):
        #creates arrays to prepare for center of mass velocity calculations
        nummol = int(max(mol))
        comvx = np.zeros(nummol)
        comvy = np.zeros(nummol)
        comvz = np.zeros(nummol)
        
        molmass = np.zeros(nummol)
        for atom in range(0,n):
            molmass[mol[atom]-1] += atommass[atype[atom]]
        jx = np.zeros((len(moltypel),num_timesteps))
        jy = np.zeros((len(moltypel),num_timesteps))
        jz = np.zeros((len(moltypel),num_timesteps))
        return (comvx, comvy, comvz, nummol, molmass, jx, jy, jz)
        
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
        
    def calcCOMv(self, comvx, comvy, comvz, vx, vy, vz, mol, atype, atommass, molmass, n, nummol):
        #calculates the center of mass velocity of all molecules for the timestep
        amass = np.zeros(n)
        for i in range(0,n):
            amass[i] = atommass[atype[i]]
        (comvxt, comvyt, comvzt)= calccomf.calccom(n, nummol, vx, vy, vz, mol, amass, molmass, 0, 0, 0, 100000, 100000, 100000)
        comvx = copy.deepcopy(comvxt)
        comvy = copy.deepcopy(comvyt)
        comvz = copy.deepcopy(comvzt)
        
        return (comvx, comvy, comvz)
        
    def findconvergance(self, J,tol):
        #Determines range of indeces for fitting the integral
        #The range is for when the sum of the square of five consecutive values is less than the tolerance
        converged = False
        i = 3
        while not converged:
            if i >= len(J):
                print('J did not converge')
                begcon = i-1
                converged = True
            else:
                test = J[i-3]**2 + J[i-2]**2 + J[i-1]**2 + J[i]**2
                if test < tol * J[0]**2:
                    begcon = i
                    i += 1
                    converged = True
                else:
                    i += 1
        
        while converged:
            if i >= len(J):
                endcon = i-1
                converged = False
            else:
                test = J[i-3]**2 + J[i-2]**2 + J[i-1]**2 + J[i]**2
                if test > tol * J[0]**2:
                    endcon = i
                    converged = False
                else:
                    i += 1
        return (begcon, endcon)
        
    def getchargearrays(self, molcharges, moltype):
        #generates an array with the charge on each molecule
        dotlist = np.zeros((int(max(moltype)+1),int(len(molcharges))))
        for i in range(0,len(dotlist)):
            for k in range(0,len(moltype)):
                if moltype[k] == i:
                    dotlist[i][k] = molcharges[k]
                    
        return dotlist
        
    def writeJ(self,J, tsjump, dt):
        #writes the values of the charge flux correlation function to J.dat
        #First column is time, second is the total correlation function
        #Remaining columns are the contribution of different molecule types
        outfile = open('J.dat','w')
        for i in range(0,len(J[0])):
            outfile.write('{}\t{}'.format(tsjump*dt*i,J[0][i]))
            for k in range(1,len(J)):
                outfile.write('\t{}'.format(J[k][i]))
            outfile.write('\n')
    
    def correlate(self,a,b):
        #Use fast Fourier transforms to calculate the correlation function between a and b
        #Zeros are added to the end of the arrays to avoid wrapping the calculation
        al = np.concatenate((a,np.zeros(len(a))),axis=0)
        bl = np.concatenate((b,np.zeros(len(b))),axis=0)
        c= np.fft.ifft(np.fft.fft(al)*np.conjugate(np.fft.fft(bl))).real
        d = c[:int(len(c)/2)]
        d/=(np.arange(len(d))+1)[::-1]
        return d
