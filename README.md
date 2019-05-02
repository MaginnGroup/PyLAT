# PyLAT

This is the github repository for the Python LAMMPS Analysis Tools. 

Some of the properties utilize Fortran for speedup. Before using this code, run either the shell script compile.sh or the python script compile.py

## Requirements
numpy>=1.14.1

scipy>=1.0.0

The following information can be accessed by running the command "python PyLAT.py -h"

usage: 

                PyLAT.py [-h] [-p PATH] [-i INPUT] [-f OUTPUT] [-g] [-c] [-m] [-d]
                [--NEC] [--GKC] [--IPL] [--DEC] [--DS] [--Visc] --mol MOL
                --nummol NUMMOL [-T TEMP] [-v VERBOSE] [--MSD_skip MSD_SKIP]
                [--MSD_num_init MSD_NUM_INIT] [--D_Tolerance D_TOLERANCE]
                [--RDF_Timesteps RDF_TIMESTEPS] [--RDF_maxr RDF_MAXR]
                [--RDF_binsize RDF_BINSIZE] [--GKC_skip GKC_SKIP]
                [--GKC_Tolerance GKC_TOLERANCE] [--GKC_J_Output]
                [--IPL_skip IPL_SKIP] [--DEC_start DEC_START]
                [--DS_MOL1 DS_MOL1] [--DS_MOL2 DS_MOL2] [--DS_Dist DS_DIST]
                [--DS_Dist_Tol DS_DIST_TOL] [--DS_First_Frame DS_FIRST_FRAME]
                [--DS_Num_Samples DS_NUM_SAMPLES] [--Visc_Dir VISC_DIR]
                [--Visc_Num VISC_NUM] [--Visc_Skip VISC_SKIP]
                [--Visc_Num_Boots VISC_NUM_BOOTS]
                [--Visc_Samples VISC_SAMPLES] [--Visc_Plot]
                [--Visc_Guess VISC_GUESS VISC_GUESS VISC_GUESS VISC_GUESS]
                LOG DAT TRJ [TRJ ...]

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
    

optional arguments:

    -h, --help            show this help message and exit
    -v VERBOSE, --verbose VERBOSE
                        Verbosity of program: 
                        0: No print statements 
                        1: Print after each major step 
                        2: Print after each major step and progress of the major step at each minor step. 
                           Not recommended for cluster calculations

File arguments:

    LOG                   LAMMPS log file.
    DAT                   LAMMPS data dile.
    TRJ                   LAMMPS trajectory file. Can include multiple trajectory files.
    -p PATH, --path PATH  Path to directory containing files. Will also save output in same directory. 
                          Default is current directory. Example: path/to/directory/
    -i INPUT, --input INPUT
                        json file if properties have already been calculated for this system
    -f OUTPUT, --output OUTPUT
                        File name for output file. Defult file name is output.json

Possible Properties:

    -g, --RDF             Radial Distribution Function.
    -c, --coord           Coordination Number. Will determine the first three minima of the RDF 
                          and calculates the coordination numbers integrating to those points.
    -m, --MSD             Mean Square Displacement.
    -d, --D               Diffusivity using MSD. Automatically includes MSD calculations.
    --NEC                 Nernst Einstein ionic conductivity. 
                          Automatically includes MSD and diffusivity calculations.
    --GKC                 Green-Kubo ionic conductivity.
    --IPL                 Ion Pair Lifetime.
    --DEC                 Dielectric Constant
    --DS                  Search for molecules of given types at a given distance
    --Visc                Calculate the viscosity of the system. 
                          Requires multiple log files in different directories. 
                          The name of the log files will be taken from the LOG option

System Properties:

    --mol MOL             Add molecule name to list of molecules. Required
    --nummol NUMMOL       Add number of molecule type to list. Must have same number as --mol. Required
    -T TEMP, --temp TEMP  Temperature of system. Required for conductivity calculations

MSD Options:

    --MSD_skip MSD_SKIP   Number of timesteps to skip at the beginning of the 
                          trajectory file before calculating MSD. 
                          Default is 0
    --MSD_num_init MSD_NUM_INIT
                        Number of initial timesteps to consider in MSD calculation. 
                        Default is half of frames being used in MSD calculation

Diffusivity Options:

    --D_Tolerance D_TOLERANCE
                        Tolerance for determining linear region in Log-Log plot. Default is 0.075

RDF Options:

    --RDF_Timesteps RDF_TIMESTEPS
                        Number of timesteps to use in the RDF calculation. 
                        Will use the last n timesteps. 
                        Default is to use all timesteps
    --RDF_maxr RDF_MAXR   Maximum r for RDF calculation. Default is half the shortest box length.
    --RDF_binsize RDF_BINSIZE
                        Bin size for RDF calculation. Default is 0.1

GK Conductivity options:
  
    --GKC_skip GKC_SKIP   Number of timesteps to skip at the beginning of the trajectory file 
                          before calculating GK Conductivity. Default is 0
    --GKC_Tolerance GKC_TOLERANCE
                        Tolerance for finding the converged region for GK conductivity calculation. 
                        Default is 0.001
    --GKC_J_Output        Option to output the charge flux correlation function. 
                          If included, a file J.dat will be created with the values 
                          of the correlation function.

Ion Pair Lifetime options:

    --IPL_skip IPL_SKIP   Number of timesteps to skip at the beginning of the trajectory file 
                          before calculating Ion Pair Lifetime. Default is 0

Dielectric Constant options:

    --DEC_start DEC_START
                        Starting frame for dielectric constant calculations. Default is 1

Distance Search Options:

    --DS_MOL1 DS_MOL1     Molecule name for the first molecule type. Required for Distance Search
    --DS_MOL2 DS_MOL2     Molecule name for the second molecule type. Required for Distance Search
    --DS_Dist DS_DIST     Distance between molecules. Required for Distance Search
    --DS_Dist_Tol DS_DIST_TOL
                        Tolerance for distance. Required for Distance Search
    --DS_First_Frame DS_FIRST_FRAME
                        First frame for distance search calculation. Default is 1
    --DS_Num_Samples DS_NUM_SAMPLES
                        Number of samples to find. Default is 1

Viscosity options:
  
    --Visc_Dir VISC_DIR   Basename for the viscosity directories. 
                          For example, if basename is "visc" the directories should be visc1, visc2, ... 
                          If no base name is given, directories should be 1, 2, 3, ...
    --Visc_Num VISC_NUM   Number of Trajectories to use in the viscosity calculation. 
                          Each trajectory should have a seperate trajectory
    --Visc_Skip VISC_SKIP
                          Number of timesteps to skip in the viscosity calculation. 
                          Default of 0 should only be used if the simulations were equilibrated 
                          before the start of the production run
    --Visc_Num_Boots VISC_NUM_BOOTS
                        Number of times to sample trajectories for bootstrapping. Default is 10
    --Visc_Samples VISC_SAMPLES
                        Number of samples for each bootstrapping iteration. Default is 30
    --Visc_Plot           Include plotting for viscosity fitting. 
                          Matplotlib required for fitting. Will output to x11
    --Visc_Guess VISC_GUESS VISC_GUESS VISC_GUESS VISC_GUESS
                        Initial guess for viscosity fitting parameters. 
                        Must add all four parameters with one call. 
                        Order of parameters is A, alpha, tau1, tau2. 
                        Default is the equivalent of calling "--Visc_Guess 1e-3 1.5e-1 1e2 1e3".

    
    Example: For diffusivity and RDF of a system with 500 water and 500 
        methanol in the directory methanol/inwater
    
    python PyLAT.py -d -g --mol H2O --mol CH3OH --nummol 500 
        --nummol 500 -p methanol/inwater/ -f methanolwater.json -v 2 mol.log 
        restart.dat mol1.lammpstrj mol2.lammpstrj
    

    
