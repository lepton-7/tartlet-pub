#!/bin/bash

#SBATCH --time=00:10:00
#SBATCH --nodes=2 --cpus-per-task=2

#SBATCH --account=PDS0325
#SBATCH --mail-type=BEGIN,END,FAIL

# Change the output file to reflect input directory used
## sbatch --output=jobs/2.record_downstream.out.$(date +"%Y-%m-%d_%H-%M-%S").%j 2.record_downstream.sh

set echo on

RT="$HOME/packages/tart"
PROD="$RT/validation/orfs/calls"
LEDGER="$RT/validation/tables/inf_results.csv"

mpiexec tart-utils record-downstream --ledger $LEDGER -i $PROD &&
    wait

echo "Done"
