lte_mpi
=======

Parallel LTE synthesis

This code carries out the synthesis of the spectrum in any given
spectral window using as many computing nodes as available. 

Compilation
-----------

In order to compile the code, you need to have an MPI distribution installed
in your system and a suitable Fortran compiler. I have only tested it with GFortran
and MPICH. To compile the code, just enter in the Source directory and type:

 make
 
Running the code
----------------

The code is run with the following command

 ./run.py conf.ini n
 
where conf.ini is the configuration file and n is the number of cpus you want to use
for the synthesis.

The code needs a list of lines that should include all atomic and molecular lines of interest.
The necessary data to generate these files is in the DATA subdirectory. Initially, only molecular
lines are available (kurucz_expanded.list). There is a small IDL script (merge.pro) that allows
the user to generate a larger file with CO, OH and TiO lines added. Just run the routine "merge"
and it will generate the appropriate files.

Configuration
-------------

The configuration of each run is modified through a simple configuration file (conf.ini in the
previous example). The configuration should be self-explanatory.