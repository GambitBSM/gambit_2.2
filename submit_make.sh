#!/bin/bash
#SBATCH -p batch   
##SBATCH --qos=xlong
#SBATCH --mem=8GB         
#SBATCH -n 32
##SBATCH -N 1
#SBATCH --time=2:00:00
CB
cd /fast/users/a1607156/gambitgit/build
make -j32 gambit 