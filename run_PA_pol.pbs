#!/bin/sh
## job name
#PBS -N GA_PA_pol_cap
## maximum runtime
#PBS -l walltime=24:00:00
## select number of nodes, cpus (per node) and memory (per node)
#PBS -l select=5:ncpus=24:mpiprocs=21:mem=100gb
### for 6-29
####PBS -l select=10:ncpus=24:mpiprocs=11:mem=100gb
### for 1-5
## standard output standard error
#PBS -o reports
#PBS -e reports
## request disk space for temporary directory
###PBS -l tmpspace=10gb
## array job
###PBS -J 6-29
PBS_ARRAY_INDEX=12
## start job after another job finished?
###PBS -W depend=afterany:2521439.pbs


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
### multiplier
#mpiexec R CMD BATCH --vanilla --quiet "--args use_MPI=TRUE n_workers=0 n_blocks=1 popSize=100 maxiter=100 run=10 stock_id=$PBS_ARRAY_INDEX n_iter=500 n_yrs=50 fhist='one-way' catch_rule='catch_rule' ga_search=TRUE lag_idx=FALSE range_idx_1=FALSE range_idx_2=FALSE range_catch=FALSE exp_r=FALSE exp_f=FALSE exp_b=FALSE  interval=FALSE multiplier=TRUE upper_constraint=FALSE lower_constraint=FALSE  obj_SSB=FALSE obj_F=FALSE obj_C=FALSE obj_risk=FALSE obj_ICV=FALSE obj_ICES_PA=FALSE obj_ICES_PA2=FALSE obj_ICES_MSYPA=TRUE collate=TRUE scenario='PA' stat_yrs='more' add_suggestions=TRUE" $HOME/git/wklife9_GA_tmp/run_ms.R $HOME/reports/$PBS_JOBID.$PBS_ARRAY_INDEX.Rout
### cap
mpiexec R CMD BATCH --vanilla --quiet "--args use_MPI=TRUE n_workers=0 n_blocks=1 popSize=100 maxiter=100 run=10 stock_id=$PBS_ARRAY_INDEX n_iter=500 n_yrs=50 fhist='one-way' catch_rule='catch_rule' ga_search=TRUE lag_idx=FALSE range_idx_1=FALSE range_idx_2=FALSE range_catch=FALSE exp_r=FALSE exp_f=FALSE exp_b=FALSE  interval=FALSE multiplier=FALSE upper_constraint=TRUE lower_constraint=TRUE  obj_SSB=FALSE obj_F=FALSE obj_C=FALSE obj_risk=FALSE obj_ICV=FALSE obj_ICES_PA=FALSE obj_ICES_PA2=FALSE obj_ICES_MSYPA=TRUE collate=TRUE scenario='PA' stat_yrs='more' add_suggestions=TRUE" $HOME/git/wklife9_GA_tmp/run_ms.R $HOME/reports/$PBS_JOBID.$PBS_ARRAY_INDEX.Rout
### multiplier & cap
#mpiexec R CMD BATCH --vanilla --quiet "--args use_MPI=TRUE n_workers=0 n_blocks=1 popSize=100 maxiter=100 run=10 stock_id=$PBS_ARRAY_INDEX n_iter=500 n_yrs=50 fhist='one-way' catch_rule='catch_rule' ga_search=TRUE lag_idx=FALSE range_idx_1=FALSE range_idx_2=FALSE range_catch=FALSE exp_r=FALSE exp_f=FALSE exp_b=FALSE  interval=FALSE multiplier=TRUE upper_constraint=TRUE lower_constraint=TRUE  obj_SSB=FALSE obj_F=FALSE obj_C=FALSE obj_risk=FALSE obj_ICV=FALSE obj_ICES_PA=FALSE obj_ICES_PA2=FALSE obj_ICES_MSYPA=TRUE collate=TRUE scenario='PA' stat_yrs='more' add_suggestions=TRUE" $HOME/git/wklife9_GA_tmp/run_ms.R $HOME/reports/$PBS_JOBID.$PBS_ARRAY_INDEX.Rout
### all without cap
#mpiexec R CMD BATCH --vanilla --quiet "--args use_MPI=TRUE n_workers=0 n_blocks=1 popSize=100 maxiter=100 run=10 stock_id=$PBS_ARRAY_INDEX n_iter=500 n_yrs=50 fhist='one-way' catch_rule='catch_rule' ga_search=TRUE lag_idx=TRUE range_idx_1=TRUE range_idx_2=TRUE range_catch=FALSE exp_r=TRUE exp_f=TRUE exp_b=TRUE  interval=TRUE multiplier=TRUE upper_constraint=FALSE lower_constraint=FALSE  obj_SSB=FALSE obj_F=FALSE obj_C=FALSE obj_risk=FALSE obj_ICV=FALSE obj_ICES_PA=FALSE obj_ICES_PA2=FALSE obj_ICES_MSYPA=TRUE collate=TRUE scenario='PA' stat_yrs='more' add_suggestions=TRUE" $HOME/git/wklife9_GA_tmp/run_ms.R $HOME/reports/$PBS_JOBID.$PBS_ARRAY_INDEX.Rout
### all with cap
#mpiexec R CMD BATCH --vanilla --quiet "--args use_MPI=TRUE n_workers=0 n_blocks=1 popSize=100 maxiter=100 run=10 stock_id=$PBS_ARRAY_INDEX n_iter=500 n_yrs=50 fhist='one-way' catch_rule='catch_rule' ga_search=TRUE lag_idx=TRUE range_idx_1=TRUE range_idx_2=TRUE range_catch=FALSE exp_r=TRUE exp_f=TRUE exp_b=TRUE  interval=TRUE multiplier=TRUE upper_constraint=TRUE lower_constraint=TRUE  obj_SSB=FALSE obj_F=FALSE obj_C=FALSE obj_risk=FALSE obj_ICV=FALSE obj_ICES_PA=FALSE obj_ICES_PA2=FALSE obj_ICES_MSYPA=TRUE collate=TRUE scenario='PA' stat_yrs='more' add_suggestions=TRUE" $HOME/git/wklife9_GA_tmp/run_ms.R $HOME/reports/$PBS_JOBID.$PBS_ARRAY_INDEX.Rout


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

