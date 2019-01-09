# -*- coding: utf-8 -*-
"""
Created on Thu May  7 15:06:47 2015

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

class calcNEconductivity:
    
    def calcNEconductivity(self, output, molcharge, Lx, Ly, Lz, nummoltype, moltypel, T):
        '''
        This function uses the Nernst-Einstien equation to estimate the ionic 
        conductivity of the system from the diffusivities
        
        '''
        V = Lx*Ly*Lz*10**-30
        e = 1.60217657e-19
        k = 1.3806488e-23
        NEcond = 0
        for i in range(0,len(moltypel)):
            q = float(molcharge[moltypel[i]])
            if q != 0:
                try:
                    D = float(output['Diffusivity'][moltypel[i]])
                except ValueError:
                    output['Nernst Einstien Conductivity in S/m'] = 'runtime not long enough'
                    return output
                N = int(nummoltype[i])
                NEcond += N*q**2*D
                output['Conductivity']['Nernst_Einstein_{0}'.format(moltypel[i])]=N*q**2*D*e**2/k/T/V
        NEcond *= e**2/k/T/V
        
        output['Conductivity']['Nernst_Einstein'] = NEcond
        
        
        return output
