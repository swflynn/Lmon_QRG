#!/bin/bash
#$ -N intensity
#$ -q free64
#$ -ckpt blcr

module load intel-parallel-studio-xe/2018.3.051
\time -o timeout ~/intensity_Lmon/src/intensity < input
