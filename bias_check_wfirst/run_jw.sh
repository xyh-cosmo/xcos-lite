#!/bin/bash

################################################################################
#   Youhua Xu
#   @2017-11-30
################################################################################

# total number of mock samples
NumMockTot=1000

# counter (start from 1)
n=1

# current work directory
PWD=`pwd`

# CLASS_INI directory
DIR_INPUTS="$PWD/w_inputs"

if [ -d $DIR_INPUTS ]
then
    echo "folder $DIR_INPUTS already exist"
else
    echo "folder $DIR_INPUTS does not exist, so I will create it"
    mkdir $DIR_INPUTS
    touch $DIR_INPUTS"/.base"
fi

# JOB script directory
DIR_JOBS="$PWD/jobs"

if [ -d $DIR_JOBS ]
then
    echo "folder $DIR_JOBS already exist"
else
    echo "folder $DIR_JOBS does not exist, so I will create it"
#    echo "mkdir $DIR_JOBS ..."
    mkdir $DIR_JOBS
    touch $DIR_JOBS"/.base"
fi

# LOGFILE directory
DIR_LOGS="$PWD/logs"

if [ -d $DIR_LOGS ]
then
    echo "folder $DIR_LOGS already exist"
else
    echo "folder $DIR_LOGS does not exist, so I will create it"
#    echo "mkdir $DIR_LOGS ..."
    mkdir $DIR_LOGS
    touch $DIR_LOGS"/.base"
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
CHAIN_DIR="$PWD/E_from_jla_wfirst_1000/chain_$n"
if [ -d $CHAIN_DIR ]
then
    echo "folder $CHAIN_DIR already exist"
else
    echo "folder $CHAIN_DIR does not exist, so I will create it"
    mkdir $CHAIN_DIR
fi

# create input parameter files (*.ini) and place them in chain_n
# 1) update the output chain root
# 2) update SN dataset file

CURRENT_INPUT="$DIR_INPUTS/bias_test_$n.ini"
CHAIN_ROOT="$CHAIN_DIR/Hz_from_JLA_WFIRST"
KEY_CHAIN_ROOT="chain_root = $CHAIN_ROOT"
KEY_WFIRST_MOCK="mock_wfirst_datafile = $PWD/SNe4Hz/mock_WFIRST_1000/WFIRST_SN_$n.txt"
KEY_JLA_MOCK="mock_snls_jla_datafile = $PWD/SNe4Hz/mock_JLA_1000/MOCK_JLA_$n.txt"
KEY_JLA_COV="mock_snls_jla_covmat = $PWD/SNe4Hz/data/JLA_cov.txt"
m4  -DKEY_CHAIN_ROOT="$KEY_CHAIN_ROOT" \
    -DKEY_WFIRST_MOCK="$KEY_WFIRST_MOCK" \
	-DKEY_JLA_MOCK="$KEY_JLA_MOCK" \
	-DKEY_JLA_COV="$KEY_JLA_COV" \
    $TEMP_INPUT > $CURRENT_INPUT

# create job script:
# 1) update KEY_JOB_ID
# 2) update KEY_LOGFILE
# 3) update KEY_INI

CURRENT_JOB="$DIR_JOBS/job.$n"

python gen_job.py "Hz_$n" "$DIR_LOGS/Hz_log_$n" "$CURRENT_INPUT" temp_job
mv temp_job $CURRENT_JOB

# Now submit the job and wait for the results
qsub -q x120.q -pe mpich 4 $CURRENT_JOB

# update counter
#echo "> n = $n"
n=`expr 1 + $n`

sleep 5 
done
