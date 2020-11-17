#!/bin/sh
## job name
#PBS -N GA_MSY_ow_group_low
## maximum runtime
#PBS -l walltime=06:00:00
## select number of nodes, cpus (per node) and memory (per node)
###PBS -l select=7:ncpus=24:mpiprocs=15:mem=100gb
### for 21-29
###PBS -l select=13:ncpus=24:mpiprocs=8:mem=100gb
### for 13-20
#PBS -l select=13:ncpus=24:mpiprocs=4:mem=100gb
### for 1-12
## standard output standard error
#PBS -o reports
#PBS -e reports
## request disk space for temporary directory
###PBS -l tmpspace=10gb
## array job
####PBS -J 1-4
####PBS_ARRAY_INDEX=5
## start job after another job finished?
#PBS -W depend=afterany:2542173.pbs


### print details about job
echo ""
echo "This is job $PBS_JOBID index $PBS_ARRAY_INDEX"
echo "The following ressources have been allocated"
cat $PBS_NODEFILE
echo ""

### set working directory
cd $HOME/git/wklife9_GA_tmp

## load modules
## anaconda includes R and OpenMPI
module purge
module load mpi/intel-2018.1.163 anaconda3/personal
### activate MPI environment
source activate R_2020

echo "starting the simulations..."
### run job
### when running with MPI, worker count is one less than selected cores
mpiexec R CMD BATCH --vanilla --quiet "--args use_MPI=TRUE n_workers=0 n_blocks=1 popSize=100 maxiter=100 run=10 stock_id=1:12 n_iter=500 n_yrs=50 fhist='one-way' catch_rule='catch_rule' ga_search=TRUE lag_idx=TRUE range_idx_1=TRUE range_idx_2=TRUE range_catch=FALSE exp_r=TRUE exp_f=TRUE exp_b=TRUE interval=TRUE multiplier=TRUE upper_constraint=FALSE lower_constraint=FALSE obj_SSB=TRUE obj_F=FALSE obj_C=TRUE obj_risk=TRUE obj_ICV=TRUE obj_ICES_PA=FALSE obj_ICES_PA2=FALSE obj_ICES_MSYPA=FALSE collate=TRUE scenario='MSY' stat_yrs='all' add_suggestions=FALSE" $HOME/git/wklife9_GA_tmp/run_ms.R $HOME/reports/$PBS_JOBID.$PBS_ARRAY_INDEX.Rout
## $PBS_ARRAY_INDEX

echo ""
echo "R job finished!"
echo ""

## print details about job
echo "job details from PBS:"
echo "==============================================================================="
qstat -f
echo "==============================================================================="
qstat -t
echo "==============================================================================="

