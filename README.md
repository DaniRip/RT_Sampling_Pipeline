# Automated Radiation Plan Evaluation

## Installation guide for Windows

1. Install CPLEX

* Make sure to keep track of the path to the following folders

* `C:\…\CPLEX_Studio201\cplex`

* `C:\…\CPLEX_Studio201\concert`

2. Install the Visual Studio C/C++ Compiler

* This application relies on the use of MATLAB mex files, which will not compile properly on other compilers like gcc or MINGW

* In MATLAB, make sure the default C++ compiler is set to the Visual Studio C++ compiler

3. Run mex -setup cpp


## Compiling and running the application

1. To compile the mex file, the following line is run via MATLAB

* `mex('-IC:\Program Files\IBM\ILOG\CPLEX_Studio201\cplex\include','-IC:\Program Files\IBM\ILOG\CPLEX_Studio201\concert\include','-LC:\Program Files\IBM\ILOG\CPLEX_Studio201\concert\lib\x64_windows_msvc14\stat_mda','-LC:\Program Files\IBM\ILOG\CPLEX_Studio201\cplex\lib\x64_windows_msvc14\stat_mda','-lcplex2010.lib','-lilocplex.lib','-lconcert.lib','run_FMO.cpp')`

* If the mexfile was successfully compiled, you can verify by looking for a file with a .mexw64 extension

3. To run the application, run the GUI, enter your custom parameters, and select .mat files
