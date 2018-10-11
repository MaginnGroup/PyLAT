# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 16:06:40 2017

@author: mhumbert
"""

"""
This script is used to parse/visualize the output of the PyLAT codes

To use, call python plots.py {json file} and follow prompts
"""

from matplotlib import pyplot as plt
import sys
import json
import numpy as np

def aver(x,y):
    xstart = float(raw_input("Select Starting x value (-1 for last point): "))
    if xstart==-1:
        print y[-1]
    else:
        xend = float(raw_input("Select ending x value: "))
        xstartind = None
        xendind = None
        for i in range(0,len(x)):
            if x[i] > xstart and xstartind==None:
                xstartind=i
            if x[i] > xend and xendind==None:
                xendind=i
        print np.average(y[xstartind:xendind])
    
    
def writefile(x,y,filename):
    out=open(filename,'w')
    for i in range(0,len(x)):
        out.write('{}\t{}\n'.format(x[i],y[i]))
    out.close()

def dic(data):
    keylist = data.keys()
    for i in range(0,len(keylist)):
        print('{}\t{}'.format(i,keylist[i]))
    select = int(raw_input("Select number of next level or y variable: "))
    if type(data[keylist[select]])==list:
        select2 = int(raw_input("select x variable or if no x variable select -1 for plot vs count: "))
        if select2==-1:
            plt.plot(range(0,len(data[keylist[select]])),data[keylist[select]])
            x = range(0,len(data[keylist[select]]))
            y = data[keylist[select]]
        else:
            plt.plot(data[keylist[select2]],data[keylist[select]])
            x = data[keylist[select2]]
            y = data[keylist[select]]
        plt.show()
        select3 = raw_input("Would you like to calculate an average or get last point? ")
        if select3.lower()=="yes" or select3.lower()=="y":
            aver(x,y)
        select4 = raw_input("Would you like to save a data file? ")
        if select4.lower()=="yes" or select4.lower()=="y":
            filename = raw_input("Please give file name: ")
            writefile(x, y, filename)
    elif type(data[keylist[select]])==int:
        print("You have selected an int. the value is {}".format(data[keylist[select]]))
    elif type(data[keylist[select]])==float:
        print("You have selected a float. the value is {}".format(data[keylist[select]]))
    elif type(data[keylist[select]])==dict:
        dic(data[keylist[select]])
    elif keylist[select]=="units":
        print("The units are {}.".format(data[keylist[select]]))
    elif type(data[keylist[select]])==str:
        print(data[keylist[select]])
    elif type(data[keylist[select]])==unicode:
        print(data[keylist[select]])
    else:
        print type(data[keylist[select]])
        print("Type unknown. The value is {}".format(data[keylist[select]]))
        

        
if __name__=="__main__":
    data = json.load(open(sys.argv[1]))
    dic(data)