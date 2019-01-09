# coding: utf-8

"""
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

from __future__ import division, print_function, unicode_literals, \
    absolute_import

import copy

import numpy as np
from six.moves import range

__author__ = "mhumbert"


class COMradialdistribution:
    def runradial(self, datfilename, comx, comy, comz, Lx, Ly, Lz, Lx2, Ly2,
                  Lz2, output, nummoltype, moltypel, moltype, timesteps, ver,
                  maxr, binsize):
        
        """
        
        This function calculates the radial distribution function between the 
        center of mass for all species in the system
        
        """
        (maxr, numbins, count, g,firststep) = self.setgparam(Lx2, Ly2, Lz2, timesteps,
                                                            moltypel, maxr,
                                                            binsize,len(comx))
        (count) = self.radialdistribution(g, len(comx[1]), moltype, comx,
                                              comy, comz, Lx, Ly, Lz, binsize,
                                              numbins, maxr, count)
        (radiuslist) = self.radialnormalization(numbins, binsize, Lx, Ly, Lz,
                                                nummoltype, count, g,
                                                firststep)
        self.append_dict(radiuslist, g, output, moltypel)
        return output

    def setgparam(self, Lx2, Ly2, Lz2, timesteps, moltypel, maxr, binsize,numsteps):
        # uses side lengths to set the maximum radius for box and number of bins
        # also sets the first line using data on firststep and number of atoms
        firststep = numsteps-timesteps
        if maxr == None:
            maxr = min(Lx2, Ly2, Lz2)
        else:
            maxr = float(maxr)
        numbins = int(np.ceil(maxr / binsize))
        count = firststep
        g = np.zeros((len(moltypel), len(moltypel), numbins))
        return maxr, numbins, count, g, firststep

    def radialdistribution(self, g, nummol, moltype, comx, comy, comz, Lx, Ly,
                           Lz, binsize, numbins, maxr, count):
        # calculates the number of molecules within each shell
        comxt = np.array(comx[count:]).transpose()
        comyt = np.array(comy[count:]).transpose()
        comzt = np.array(comz[count:]).transpose()
        indexlist = []
        
        # change indeces order to com*[molecule][timestep]
        
        for i in range(0,len(g)):
            indexlist.append(np.array(moltype) == i)
            #creates two dimensional array
            #first dimension is molecule type
            #second dimension is over molecules
            #contains boolean for if that molecule is of the molecule type
        
        for molecule in range(0, nummol - 1):
            dx = comxt[molecule+1:] - np.tile(comxt[molecule],
                                              (len(comxt)-molecule-1,1))
            dy = comyt[molecule+1:] - np.tile(comyt[molecule],
                                              (len(comyt)-molecule-1,1))
            dz = comzt[molecule+1:] - np.tile(comzt[molecule],
                                              (len(comzt)-molecule-1,1))

            dx -= Lx * np.around(dx / Lx)
            dy -= Ly * np.around(dy / Ly)
            dz -= Lz * np.around(dz / Lz)
            #minimum image convention

            r2 = dx ** 2 + dy ** 2 + dz ** 2
            r = np.sqrt(r2)
            for i in range(0,len(indexlist)):
                gt,dist = np.histogram(r[indexlist[i][molecule+1:]].ravel(),
                                       bins=numbins,
                                       range=(0.5*binsize,binsize*(numbins+0.5)))
                g[moltype[molecule]][i]+= gt
                g[i][moltype[molecule]]+= gt
                
        count = len(comx)
        return count

    def radialnormalization(self, numbins, binsize, Lx, Ly, Lz, nummol, count,
                            g, firststep):
        # normalizes g to box density
        radiuslist = (np.arange(numbins) + 1) * binsize
        radiuslist = np.around(radiuslist,decimals = 1)
        for i in range(0, len(g)):
            for j in range(0, len(g)):
                g[i][j] *= Lx * Ly * Lz / nummol[i] / nummol[j] / 4 / np.pi / (
                               radiuslist) ** 2 / binsize / (
                               count - firststep)
        return radiuslist

    def append_dict(self, radiuslist, g, output, moltypel):
        for i in range(0, len(moltypel)):
            for j in range(i, len(moltypel)):
                if not all([v == 0 for v in g[i][j]]):
                    output['RDF']['{0}-{1}'.format(moltypel[i],
                                                   moltypel[
                                                       j])] = copy.deepcopy(
                        g[i][j].tolist())
        if 'distance' not in list(output['RDF'].keys()):
            output['RDF']['distance'] = copy.deepcopy(radiuslist.tolist())
