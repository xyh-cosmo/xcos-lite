#!/bin/bash
qsub -q x120.q -pe mpich 20  job_h8
qsub -q x120.q -pe mpich 20  job_h10
qsub -q x120.q -pe mpich 20  job_u_h8
qsub -q x120.q -pe mpich 20  job_u_h10
qsub -q x120.q -pe mpich 20  job_j_h8
qsub -q x120.q -pe mpich 20  job_j_h10
#qsub -q b130.q -pe mpich 50  job
#qsub -q x120.q@compute-0-9.local -pe mpich 10  job
