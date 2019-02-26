# -*- coding: utf-8 -*-
"""
Created on Wed May 27 12:38:04 2015

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
import numpy
from scipy.integrate import cumtrapz
import copy

class getcoordinationnumber:
    def calccoordinationnumber(self, output, nummoltype, moltypel, V):
        
        """
        This function integrates the radial distribution function to obtain the 
        coordination number for all molecule types
        
        This code determines the first three minima in the radial distribution 
        function and returns the coordination number for these three points. 
        The cumulative integral is also returned
        """
        
        output['Coordination_Number'] = {}
        output['Coordination_Number']['units'] = 'Minima in angstroms, Coordination numbers in Angstroms'
        output['Coordination_Number']['explanation'] = 'This program finds the first three local minima and finds the coordination number integrating until there. Na-H20 represents the coordination number for water around sodium.'
        pairlist = list(output['RDF'].keys())
        pairlist.remove('units')
        pairlist.remove('distance')
        r = output['RDF']['distance']
        for i in range(0,len(pairlist)):
            g = output['RDF'][pairlist[i]]
            split = pairlist[i].split('-')
            mol1 = split[0]
            mol2 = split[1]
            (minima,index) = self.findfirst3minima(g,r)
            output['Coordination_Number']['{0} around {1}'.format(mol1, mol2)] = {}
            integral = self.integrate(g,r,nummoltype, moltypel, V, mol1)
            output['Coordination_Number']['{0} around {1}'.format(mol1, mol2)]['Cumulative_Integral'] = copy.deepcopy(integral)
            output['Coordination_Number']['{0} around {1}'.format(mol1, mol2)]['Minima'] = minima
            coord = []
            for j in range(0,len(minima)):
                coord.append(integral[index[j]])
            output['Coordination_Number']['{0} around {1}'.format(mol1, mol2)]['Coordination_Numbers'] = coord
            if mol2 != mol1:
                output['Coordination_Number']['{0} around {1}'.format(mol2, mol1)] = {}
                integral = self.integrate(g,r,nummoltype, moltypel, V, mol2)
                output['Coordination_Number']['{0} around {1}'.format(mol2, mol1)]['Cumulative_Integral'] = copy.deepcopy(integral)
                output['Coordination_Number']['{0} around {1}'.format(mol2, mol1)]['Minima'] = minima
                coord = []
                for j in range(0,len(minima)):
                    coord.append(integral[index[j]])
                output['Coordination_Number']['{0} around {1}'.format(mol2, mol1)]['Coordination_Numbers'] = coord
        return output
            
    def findfirst3minima(self, g, r):
        #determines the first three minima in the radial distribution function
        foundpositive = False
        minima = []
        index = []
        i=0
        while not foundpositive:
            if g[i] > 1:
                foundpositive = True
                i += 1
            else:
                i +=1
        
        while len(minima) < 3 and i < len(g)-2:
            if g[i-1] > g[i] and g[i+1] > g[i]:
                minima.append(r[i])
                index.append(i)
            i += 1
            
        return (minima,index)
        
        
    def integrate(self, g,r, nummoltype, moltypel, V, mol):
        #integrates the radial distribution functions
        integrallist = []
        for i in range(0,len(g)):
            integrallist.append(g[i]*nummoltype[moltypel.index(mol)]/V*4*numpy.pi*r[i]**2)
        integral = cumtrapz(integrallist, x = r)
        integral= integral.tolist()
        return integral

    def getvolume(self, trjfilename):
        #calculates the volume of the system (assumes NVT or NVE ensemble)
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
