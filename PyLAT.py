# -*- coding: utf-8 -*-
"""
Created on Wed May 20 13:37:57 2015

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

This is the driver file for the PyLAT program. This file parsers the user 
input, and then calls the proper class files
"""
if __name__ == '__main__':
    
    import argparse as args
    import sys
    import json
    
    try:
        import src.calccomf
    except ImportError:
        print("Please use either 'sh compile.sh' or 'python compile.py' to compile the fortran modules")
        sys.exit(1)
    
    #sys.path.append('~/Projects/Parser/publish/src')
    
    from src.MSD import MSD
    from src.calcDiffusivity import calcdiffusivity
    from src.calcCOM import calcCOM
    from src.getTimeData import gettimedata
    from src.getMolData import getmoldata
    from src.COMradial import COMradialdistribution
    from src.getAtomCharges import getatomcharges
    from src.calcNEconductivity import calcNEconductivity
    from src.calcCond import calcCond    
    from src.getCoordinationNumber import getcoordinationnumber
    from src.ionpair import ionpair
    from src.calcDielectricConstant import calcDielectricConstant
    from src.distSearch import distSearch
    from src.calcVisc import calcVisc
    from src.fitVisc import fitVisc
    
    parser = args.ArgumentParser(formatter_class=args.RawTextHelpFormatter, description='''
    Analyze output from LAMMPS simulation
    
    Required files are the log file, the data file and at least one trajectory
    file in that order
    
    Also required is the name and number of each molecule type in the system. 
    To do this, use the --mol and --nummol flags. Use the same order as they
    were entered into the data file. 
    
    If you want to calculate an ionic conductivity or dielectric constant, a 
    temperature must be included using a --temp flag. 
    
    The output of this program is a json file with all of the requested 
    information. Default file name is output.json. To use data in this file,
    use the following python code:
        import json
        openfile = open("filename")
        dict = json.load(openfile)
        
    This yields a dictionary with all of the information that you can query.
    
    plot.py is an interactive program written to parse and display the 
    information in the json output files. Proper usage is
    plots.py {output file}
    ''', epilog = '''
    
    Example: For diffusivity and RDF of a system with 500 water and 500 
        methanol in the directory methanol/inwater
    
    python PyLAT.py -d -g --mol H2O --mol CH3OH --nummol 500 
        --nummol 500 -p methanol/inwater/ -f methanolwater.json -v 2 mol.log 
        restart.dat mol1.lammpstrj mol2.lammpstrj
    ''')
    File_Options = parser.add_argument_group("File arguments")
    File_Options.add_argument('LOG', nargs=1, help = 'LAMMPS log file.')
    File_Options.add_argument('DAT', nargs=1, help = 'LAMMPS data dile.')
    File_Options.add_argument('TRJ', action='append', nargs='+',
                        help = 'LAMMPS trajectory file. Can include multiple trajectory files.')
    File_Options.add_argument('-p', '--path', 
                        help = 'Path to directory containing files. Will also save output in same directory. Default is current directory. Example: path/to/directory/',
                        default = './')
    File_Options.add_argument('-i', '--input', help = 'json file if properties have already been calculated for this system')
    File_Options.add_argument('-f', '--output', 
                        help = 'File name for output file. Defult file name is output.json', default = 'output.json')
    properties = parser.add_argument_group("Possible Properties")
    properties.add_argument('-g', '--RDF', help = 'Radial Distribution Function.', action = 'store_true')
    properties.add_argument('-c', '--coord', help = 'Coordination Number. Will determine the first three minima of the RDF and calculates the coordination numbers integrating to those points.', action = 'store_true')
    properties.add_argument('-m', '--MSD', help = 'Mean Square Displacement.', action = 'store_true')
    properties.add_argument('-d', '--D', 
                        help = 'Diffusivity using MSD. Automatically includes MSD calculations.', 
                        action = 'store_true')
    properties.add_argument('--NEC', 
                        help = 'Nernst Einstein ionic conductivity. Automatically includes MSD and diffusivity calculations.', 
                        action = 'store_true')
    properties.add_argument('--GKC', help = 'Green-Kubo ionic conductivity.', action = 'store_true')
    properties.add_argument('--IPL', help = 'Ion Pair Lifetime.', action = 'store_true')
    properties.add_argument('--DEC', help = 'Dielectric Constant', action = 'store_true')    
    properties.add_argument('--DS', help = 'Search for molecules of given types at a given distance', action = 'store_true')
    properties.add_argument('--Visc',help='Calculate the viscosity of the system. Requires multiple log files in different directories. The name of the log files will be taken from the LOG option', action='store_true')
    
    system = parser.add_argument_group("System Properties")    
    system.add_argument('--mol', help = 'Add molecule name to list of molecules. Required', action = 'append', required=True)
    system.add_argument('--nummol', help = 'Add number of molecule type to list. Must have same number as --mol. Required', 
                        action = 'append', required=True, type=int)
    system.add_argument('-T', '--temp', help = 'Temperature of system. Required for conductivity calculations', type=float)    
    parser.add_argument('-v', '--verbose', help = '''Verbosity of program: \n0: No print statements \n1: Print after each major step \n2: Print after each major step and progress of the major step at each minor step. Not recommended for cluster calculations''', type=int)
    
    MSDoptions = parser.add_argument_group('MSD Options')
    MSDoptions.add_argument('--MSD_skip', help = 'Number of timesteps to skip at the beginning of the trajectory file before calculating MSD. Default is 0', default = 0, type=int)
    MSDoptions.add_argument('--MSD_num_init', help = 'Number of initial timesteps to consider in MSD calculation. Default is half of frames being used in MSD calculation', default = None)    
    
    Difoptions = parser.add_argument_group('Diffusivity Options')
    Difoptions.add_argument('--D_Tolerance', help = 'Tolerance for determining linear region in Log-Log plot. Default is 0.075', type=float, default=0.075)
    
    RDFoptions = parser.add_argument_group('RDF Options')
    RDFoptions.add_argument('--RDF_Timesteps', help = 'Number of timesteps to use in the RDF calculation. \nWill use the last n timesteps. \nDefault is to use all timesteps', type = int)    
    RDFoptions.add_argument('--RDF_maxr', help = 'Maximum r for RDF calculation. Default is half the shortest box length.', default = None)
    RDFoptions.add_argument('--RDF_binsize', help = 'Bin size for RDF calculation. Default is 0.1', default = 0.1, type = float)
    
    
    GKCoptions = parser.add_argument_group('GK Conductivity options')
    GKCoptions.add_argument('--GKC_skip', help = 'Number of timesteps to skip at the beginning of the trajectory file before calculating GK Conductivity. Default is 0', default = 0, type=int)    
    GKCoptions.add_argument('--GKC_Tolerance', help= 'Tolerance for finding the converged region for GK conductivity calculation. Default is 0.001', default = 0.001, type=float)
    GKCoptions.add_argument('--GKC_J_Output', help='Option to output the charge flux correlation function. If included, a file J.dat will be created with the values of the correlation function.', action='store_true')
    
    IPLoptions = parser.add_argument_group('Ion Pair Lifetime options')
    IPLoptions.add_argument('--IPL_skip', help = 'Number of timesteps to skip at the beginning of the trajectory file before calculating Ion Pair Lifetime. Default is 0', default = 0, type=int)    
    
    DECoptions = parser.add_argument_group('Dielectric Constant options')
    DECoptions.add_argument('--DEC_start', help = 'Starting frame for dielectric constant calculations. Default is 1', default = 1, type = int)    
    
    DSoptions = parser.add_argument_group('Distance Search Options')
    DSoptions.add_argument('--DS_MOL1', help = 'Molecule name for the first molecule type. Required for Distance Search')
    DSoptions.add_argument('--DS_MOL2', help = 'Molecule name for the second molecule type. Required for Distance Search')    
    DSoptions.add_argument('--DS_Dist', help = 'Distance between molecules. Required for Distance Search', type = float)
    DSoptions.add_argument('--DS_Dist_Tol', help = 'Tolerance for distance. Required for Distance Search', type = float)
    DSoptions.add_argument('--DS_First_Frame', help = 'First frame for distance search calculation. Default is 1', default = 1, type = int)    
    DSoptions.add_argument('--DS_Num_Samples', help = 'Number of samples to find. Default is 1', default = 1, type = int)

    Viscoptions = parser.add_argument_group('Viscosity options')
    Viscoptions.add_argument('--Visc_Dir', help = 'Basename for the viscosity directories. For example, if basename is "visc" the directories should be visc1, visc2, ... \n If no base name is given, directories should be 1, 2, 3, ...',
                             default = None)
    Viscoptions.add_argument('--Visc_Num', help = 'Number of Trajectories to use in the viscosity calculation. Each trajectory should have a seperate trajectory',type=int)
    Viscoptions.add_argument('--Visc_Skip', help = 'Number of timesteps to skip in the viscosity calculation. Default of 0 should only be used if the simulations were equilibrated before the start of the production run', default=0, type=int)
    Viscoptions.add_argument('--Visc_Num_Boots', help = 'Number of times to sample trajectories for bootstrapping. Default is 10', default = 10, type = int)
    Viscoptions.add_argument('--Visc_Samples', help = 'Number of samples for each bootstrapping iteration. Default is 30', default = 30, type = int)
    Viscoptions.add_argument('--Visc_Plot', help = 'Include plotting for viscosity fitting. Matplotlib required for fitting. Will output to x11', action = 'store_true')
    Viscoptions.add_argument('--Visc_Guess', help = 'Initial guess for viscosity fitting parameters. Must add all four parameters with one call. Order of parameters is A, alpha, tau1, tau2. Default is the equivalent of calling "--Visc_Guess 1e-3 1.5e-1 1e2 1e3".', type = float, nargs = 4)
    
    arg = parser.parse_args()
    
    c = calcCOM()
    m = MSD()
    cd = calcdiffusivity()
    gt = gettimedata()
    gm = getmoldata()
    crd = COMradialdistribution()
    gc = getatomcharges()
    ne = calcNEconductivity()
    cc = calcCond()
    gcn = getcoordinationnumber()
    ip = ionpair()
    dec = calcDielectricConstant()
    ds = distSearch()
    cv = calcVisc()
    fv = fitVisc()
    
    if arg.path[-1] != '/':
        arg.path = arg.path + '/'
     
    for i in range(0,len(arg.TRJ[0])):
        arg.TRJ[0][i] = arg.path + arg.TRJ[0][i]
        
    trjfilename = arg.TRJ[0]
    datfilename = arg.path + arg.DAT[0]
    logfilename = arg.path + arg.LOG[0]
    
    if arg.verbose == 2:
        ver = True
    else:
        ver = False

        
    if arg.input == None:    
        output = {}
        
    else:
        inputfile = arg.path + arg.input
        openfile = open(inputfile, 'r')
        output = json.load(openfile)
        openfile.close()
    
    nummoltype = arg.nummol
    moltypel = arg.mol
    moltype = []
    for i in range(0,len(moltypel)):
        for j in range(0,nummoltype[i]):
            moltype.append(int(i))
       
    if arg.NEC or arg.GKC:
        if arg.temp == None:
            sys.exit('Need temperature for conductivity calculations')
        if 'Conductivity' not in output.keys():
            output['Conductivity'] = {}
            output['Conductivity']['units'] = 'S/m'
            
    if arg.DEC:
        if arg.temp == None:
            sys.exit('Need temperature for dielectric constant calculations')
        
                
    if len(arg.mol) != len(arg.nummol):
        sys.exit('Must enter same number of arguments for nummol and mol')                
                
    if arg.NEC and 'Diffusivity' not in output.keys():
        arg.D = True
    
    if arg.D and 'MSD' not in output.keys():
        arg.MSD = True
        
    if arg.coord and 'RDF' not in output.keys():
        arg.RDF = True
    
    if arg.MSD or arg.D or arg.IPL:
        tsjump = gt.getjump(trjfilename[0])
        dt = gt.getdt(logfilename)
    if arg.NEC or arg.GKC or arg.GKC or arg.DEC:
        n = gc.findnumatoms(datfilename)
        (molcharges, atomcharges,n) = gc.getmolcharges(datfilename,n)
        molcharge = gc.molchargedict(molcharges, moltypel, moltype)
    if arg.RDF or arg.MSD or arg.IPL or arg.coord or arg.DEC:
        (V, Lx, Ly, Lz) = gcn.getvolume(trjfilename[0])

    if arg.RDF or arg.MSD or arg.IPL:
        if arg.verbose >= 1:    
            print('beginning COM calculation')        
        (comx, comy, comz, Lx, Ly, Lz, Lx2, Ly2, Lz2) = c.calcCOM(trjfilename,datfilename, ver)
        if arg.verbose == 1:
            print('COM calculation complete')
        
    if arg.RDF:
        if arg.RDF_Timesteps == None:
            arg.RDF_Timesteps = len(comx)
        if arg.RDF_Timesteps > len(comx):
            print('Number of RDF timesteps requested longer than simulation. Using all timesteps.')
            arg.RDF_Timesteps = len(comx)
        if arg.verbose >= 1:    
            print('beginning RDF calculation')
        output['RDF'] = {}
        output['RDF']['units'] = 'unitless, angstroms'
        output = crd.runradial(datfilename, comx, comy, comz, Lx, Ly, Lz, Lx2, Ly2, Lz2, output, nummoltype, moltypel, moltype, arg.RDF_Timesteps, ver,arg.RDF_maxr,arg.RDF_binsize)
        if arg.verbose >= 1:
            print('RDF calculation complete')
    
    if arg.coord:
        if arg.verbose >= 1:    
            print('beginning coordination number calculation')
        output = gcn.calccoordinationnumber(output, nummoltype, moltypel, V)
        if arg.verbose >= 1:
            print('coordination number calculation complete')
             
    if arg.DEC:
        output = dec.calcDEC(atomcharges, trjfilename, arg.temp, output, V, arg.verbose, arg.DEC_start)
             
    if arg.GKC:
        if arg.verbose >= 1:    
            print('beginning GK conductivity calculation')
        output = cc.calcConductivity(molcharges, trjfilename, logfilename, datfilename, arg.temp, output, moltype, moltypel,arg.verbose,arg.GKC_skip,arg.GKC_Tolerance,arg.GKC_J_Output)
        if arg.verbose >= 1:
            print('GK conductivity calculation complete')
            
    if arg.IPL:
        if arg.verbose >= 1:    
            print('beginning Ion Pair Lifetime calculation')
        ip.runionpair(comx,comy,comz,Lx,Ly,Lz,moltypel,moltype,tsjump,dt,output,ver,arg.IPL_skip)
        if arg.verbose >= 1:
            print('Ion Pair Lifetime calculation complete')
            
            
    if arg.DS:
        if arg.verbose >= 1:    
            print('beginning Distance Search')
        output = ds.distSearch(arg.DS_MOL1, arg.DS_MOL2, arg.DS_Dist, arg.DS_Dist_Tol, arg.DS_First_Frame, arg.DS_Num_Samples, trjfilename, datfilename, moltype, moltypel, output)
        
        
    if arg.MSD:
        if arg.verbose >= 1:    
            print('beginning MSD calculation')
        output = m.runMSD(comx, comy, comz, Lx, Ly, Lz, Lx2, Ly2, Lz2, moltype, moltypel, dt, tsjump, output, ver, arg.MSD_skip, arg.MSD_num_init)
        if arg.verbose == 1:
            print('MSD calculation complete')
            
    if arg.D:
        if arg.verbose >= 1:    
            print('beginning diffusivity calculation')
        cd.calcdiffusivity(output, moltypel, dt,arg.D_Tolerance)
        if arg.verbose:
            print('diffusivity calculation complete')
        
    if arg.NEC:
        if arg.verbose >= 1:    
            print('beginning NE conductivity calculation')
        output = ne.calcNEconductivity(output, molcharge, Lx, Ly, Lz, nummoltype, moltypel, arg.temp)
        if arg.verbose >= 1:
            print('NE conductivity calculation complete')
            
    if arg.Visc:
        if arg.verbose >= 1:    
            print('beginning viscosity calculation')
        if arg.Visc_Guess == None:
            arg.Visc_Guess = [1e-3,1.5e-1,1e2,1e3]
        if len(arg.Visc_Guess) != 4:
            print('need four inputs for viscosity initial parameters.')
        output = cv.calcvisc(arg.Visc_Num,arg.Visc_Skip,arg.Visc_Dir,arg.LOG[0],output, arg.verbose, arg.Visc_Samples, arg.Visc_Num_Boots, arg.Visc_Plot, arg.Visc_Guess)
        #output = fv.fitvisc(output)
        if arg.verbose >= 1:
            print('viscosity calculation complete')

    if output == {}:
        sys.exit('Please select at least one property to calculate')
    
    if arg.verbose >= 1:    
        print('beginning file generation')
    outputfile = open(arg.path + arg.output, 'w')
    json.dump(output, outputfile, indent=4, sort_keys=True)
    outputfile.close()
    if arg.verbose >= 1:    
        print('file generated')
