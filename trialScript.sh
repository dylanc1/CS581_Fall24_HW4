#!/bin/bash
source /apps/profiles/modules_asax.sh.dyn
module load openmpi/4.1.4-gcc11

mpicc -g -Wall -o hw4-v1 hw4-v1.c
mpicc -g -Wall -o hw4-v2 hw4-v2.c

mpiexec -n 1 ./hw4-v1 5000 5000 /scratch/$USER
mpiexec -n 2 ./hw4-v1 5000 5000 /scratch/$USER
mpiexec -n 4 ./hw4-v1 5000 5000 /scratch/$USER
mpiexec -n 8 ./hw4-v1 5000 5000 /scratch/$USER
mpiexec -n 10 ./hw4-v1 5000 5000 /scratch/$USER
mpiexec -n 16 ./hw4-v1 5000 5000 /scratch/$USER
mpiexec -n 20 ./hw4-v1 5000 5000 /scratch/$USER

mpiexec -n 1 ./hw4-v2 5000 5000 /scratch/$USER
mpiexec -n 2 ./hw4-v2 5000 5000 /scratch/$USER
mpiexec -n 4 ./hw4-v2 5000 5000 /scratch/$USER
mpiexec -n 8 ./hw4-v2 5000 5000 /scratch/$USER
mpiexec -n 10 ./hw4-v2 5000 5000 /scratch/$USER
mpiexec -n 16 ./hw4-v2 5000 5000 /scratch/$USER
mpiexec -n 20 ./hw4-v2 5000 5000 /scratch/$USER


diff /scratch/$USER/output.5000.5000.1 /scratch/$USER/output.5000.5000.2
diff /scratch/$USER/output.5000.5000.1 /scratch/$USER/output.5000.5000.4
diff /scratch/$USER/output.5000.5000.1 /scratch/$USER/output.5000.5000.8
diff /scratch/$USER/output.5000.5000.1 /scratch/$USER/output.5000.5000.10
diff /scratch/$USER/output.5000.5000.1 /scratch/$USER/output.5000.5000.16
diff /scratch/$USER/output.5000.5000.1 /scratch/$USER/output.5000.5000.20
diff /scratch/$USER/output.5000.5000.1 /scratch/$USER/output.5000.5000.1-2
diff /scratch/$USER/output.5000.5000.1 /scratch/$USER/output.5000.5000.2-2
diff /scratch/$USER/output.5000.5000.1 /scratch/$USER/output.5000.5000.4-2
diff /scratch/$USER/output.5000.5000.1 /scratch/$USER/output.5000.5000.8-2
diff /scratch/$USER/output.5000.5000.1 /scratch/$USER/output.5000.5000.10-2
diff /scratch/$USER/output.5000.5000.1 /scratch/$USER/output.5000.5000.16-2
diff /scratch/$USER/output.5000.5000.1 /scratch/$USER/output.5000.5000.20-2
