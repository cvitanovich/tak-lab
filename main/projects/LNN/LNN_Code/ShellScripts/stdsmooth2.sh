#!/bin/sh
#
# job submission script for matlab
# Request 1 node
#PBS -l nodes=1:ppn=8,walltime=7:00:00:00,other=matlab
# Join the STDOUT and STDERR output streams into STDOUT
#PBS -j oe
# Write STDOUT to stdsmooth2.txt
#PBS -o stdsmooth2.txt
# Send mail when the job aborts or terminate
#
#PBS -m ae
#
# Mail address
#
#PBS -M cvitanov@uoregon.edu

# Add matlab module

#. use_modules
#module add matlab

# Goto the directory from which you submitted the job

cd $PBS_O_WORKDIR

# Specify the executable ...

matlab -nodisplay -nosplash < stdsmooth2.m > stdsmooth2.dat
