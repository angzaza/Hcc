#!/bin/bash
# Usage:
#  submitAllJobs.sh


echo -e "\npython createRunFile_new.py data_control control --run 2023C --n 30 --outName ReReco22Sep_v1Jan24_SF "
python createRunFile_new.py data_control control --run 2023C --n 100 --outName ReReco22Sep_v1Jan24_SF
sleep 1
pwd
echo -e "\nsubmit_analysis_2023C_control_ReReco22Sep_v1Jan24_SF.sh"
source submit_analysis_2023C_control_ReReco22Sep_v1Jan24_SF.sh
cd ..
sleep 1

echo -e "\npython createRunFile_new.py MC_control control --outName Summer23_v1Jan24_SF --n 30 --run 2023C --MCprocess QCD"
python createRunFile_new.py MC_control control --outName Summer23_v1Jan24_SF --n 30 --run 2023C --MCprocess QCD
sleep 1
echo -e "\nsubmit_analysis_QCD_control_Summer23_v1Jan24_SF.sh"
source submit_analysis_QCD_control_Summer23_v1Jan24_SF.sh
cd ..
sleep 1
