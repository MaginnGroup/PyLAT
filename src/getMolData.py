# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 10:36:25 2015

@author: mhumbert
"""


class getmoldata:
    
    def getmoltype(self, datfilename):
        # determines molecule types and number of each molecule type
        #also creates a list of molecule type of each molecule
        datfile = open(datfilename)
        datfile.readline()
        datfile.readline()
        nummoltype = []
        moltypel = []
        moltype = []
        readingmolecules = True
        while readingmolecules == True:
            line = datfile.readline()
            line = line.split()
            if len(line) == 4:
                nummoltype.append(int(line[1]))
                moltypel.append(line[2])
                
            else:
                readingmolecules = False
                
        for i in range(0,len(moltypel)):
            for j in range(0,nummoltype[i]):
                moltype.append(int(i))
        
        datfile.close()
        return (nummoltype, moltypel, moltype)
        
