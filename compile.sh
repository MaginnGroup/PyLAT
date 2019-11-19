cd src
python -m numpy.f2py -c -m calccomf calcCOM.f90
python -m numpy.f2py -c -m calcdistances calcdistances.f90
python -m numpy.f2py -c -m ipcorr ipcorr.f90
#f2py -c -m src/elradial src/elradial.f90
#f2py -c -m src/siteradial src/siteradial.f90
cd ..
