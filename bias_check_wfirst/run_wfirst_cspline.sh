#!/bin/bash

################################################################################
#   Youhua Xu
#   @2017-11-30
################################################################################

# total number of mock samples
NumMockTot=1000

# counter (start from 1)
n=101

# current work directory
PWD=`pwd`

# remove '/export'
PWD_OLD=$PWD
PWD=${PWD#*/export}
echo "reset current working dir from $PWD_OLD to $PWD"

# CLASS_INI directory
DIR_INPUTS="$PWD/w_fixed_num_inputs"

if [ -d $DIR_INPUTS ]
then
    echo "folder $DIR_INPUTS already exist"
else
    echo "folder $DIR_INPUTS does not exist, so I will create it"
    mkdir $DIR_INPUTS
fi

# JOB script directory
DIR_JOBS="$PWD/jobs_jw"

if [ -d $DIR_JOBS ]
then
    echo "folder $DIR_JOBS already exist"
else
    echo "folder $DIR_JOBS does not exist, so I will create it"
    mkdir $DIR_JOBS
fi

# LOGFILE directory
DIR_LOGS="$PWD/logs_jw"

if [ -d $DIR_LOGS ]
then
    echo "folder $DIR_LOGS already exist"
else
    echo "folder $DIR_LOGS does not exist, so I will create it"
    mkdir $DIR_LOGS
fi

# chain directory
CHAIN_DIR="WFIRST_fixed_Num"
if [ -d $CHAIN_DIR ]
then
    echo "folder $CHAIN_DIR already exist"
else
    echo "folder $CHAIN_DIR does not exist, so I will create it"
    mkdir $CHAIN_DIR
fi

# some template *ini files
TEMP_INPUT="DH_WFIRST_mock.ini"

# some key parameters in the input *ini file
CHAIN_ROOT=""
KEY_CHAIN_ROOT=""

# some key parameters in the job script
KEY_JOB_ID=""
KEY_LOGFILE=""
KEY_INI=""


while [ $n -le $NumMockTot ]
do
# create chain folder for n-th reconstruction
CHAIN_DIR_TMP="$PWD/$CHAIN_DIR/chain_$n"
if [ -d $CHAIN_DIR_TMP ]
then
    echo "folder '$CHAIN_DIR_TMP' already exist"
else
    echo "folder '$CHAIN_DIR_TMP' does not exist, so I will create it"
    mkdir $CHAIN_DIR_TMP
fi

# create input parameter files (*.ini) and place them in chain_n
# 1) update the output chain root
# 2) update SN dataset file

CURRENT_INPUT="$DIR_INPUTS/bias_test_cspline_$n.ini"
CHAIN_ROOT="$CHAIN_DIR_TMP/Hz_from_WFIRST_cspline"
KEY_CHAIN_ROOT="chain_root = $CHAIN_ROOT"
KEY_WFIRST_MOCK="mock_wfirst_datafile = $PWD/SNe4Hz.fixed_Num/mock_WFIRST_1000/WFIRST_SN_$n.txt"
KEY_JLA_MOCK="mock_snls_jla_datafile = $PWD/SNe4Hz.fixed_Num/mock_JLA_1000/MOCK_JLA_$n.txt"
KEY_JLA_COV="mock_snls_jla_covmat = $PWD/SNe4Hz.fixed_Num/data/JLA_cov.txt"
KEY_E_INTERP="E_table_interp_method = cspline"
m4  -DKEY_CHAIN_ROOT="$KEY_CHAIN_ROOT" \
    -DKEY_WFIRST_MOCK="$KEY_WFIRST_MOCK" \
	-DKEY_E_INTERP="$KEY_E_INTERP" \
    $TEMP_INPUT > $CURRENT_INPUT

# create job script:
# 1) update KEY_JOB_ID
# 2) update KEY_LOGFILE
# 3) update KEY_INI

CURRENT_JOB="$DIR_JOBS/job_cspline.$n"

python gen_job.py "Hz_$n" "$DIR_LOGS/Hz_log_cspline_$n" "$CURRENT_INPUT" temp_job
mv temp_job $CURRENT_JOB

# Now submit the job and wait for the results
qsub -q x120.q -pe mpich 4 $CURRENT_JOB

# update counter
#echo "> n = $n"
n=`expr 1 + $n`

sleep 1
done
