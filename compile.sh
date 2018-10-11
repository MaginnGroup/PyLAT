cd src
f2py -c -m calccomf calcCOM.f90
f2py -c -m calcdistances calcdistances.f90
f2py -c -m comradialall comradialall.f90
f2py -c -m ipcorr ipcorr.f90
#f2py -c -m src/elradial src/elradial.f90
#f2py -c -m src/siteradial src/siteradial.f90
cd ..
