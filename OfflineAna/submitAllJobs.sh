#!/bin/bash
# Usage:
#  submitAllJobs.sh


echo -e "\npython createRunFile_new.py data Hcc --run 2023C --n 30 --outName ReReco22Sep_noCTagCut "
python createRunFile_new.py data Hcc --run 2023C --n 100 --outName ReReco22Sep_noCTagCut
sleep 1
pwd
echo -e "\nsubmit_analysis_2023C_Hcc_ReReco22Sep_noCTagCut.sh"
source submit_analysis_2023C_Hcc_ReReco22Sep_noCTagCut.sh
cd ..
sleep 1


#for i in QCD_PT-120-170 QCD_PT-170-300 QCD_PT-300-470 QCD_PT-470-600 QCD_PT-600-800 QCD_PT-800-1000 QCD_PT-1400-1800 QCD_PT-1800-2400 QCD_PT-2400-3200
#do
# 		echo -e "\nMC $i"
#    echo -e "\npython createRunFile_new.py MC Hcc --outName v1 --n 30 --EE preEE --MCprocess $i "
#    python createRunFile_new.py MC Hcc --outName v1 --n 30 --EE preEE --MCprocess $i
#    sleep 1
#    echo -e "\nsubmit_analysis_${i}_Hcc_noCTagCut.sh"
#    source submit_analysis_${i}_Hcc_noCTagCut.sh
#    cd ..
#    sleep 1
#done


for j in Zqq_HT-200-400 Zqq_HT-400-600 Zqq_HT-600-800 Zqq_HT-800-inf
do
    echo -e "\nMC $j"
    echo -e "\npython createRunFile_new.py MC Hcc --outName Summer23_noCTagCut --n 30 --run 2023C --MCprocess $j "
    python createRunFile_new.py MC Hcc --outName Summer23_noCTagCut --n 50 --run 2023C --MCprocess $j
    sleep 1
    echo -e "\nsubmit_analysis_${j}_Hcc_Summer23_noCTagCut.sh"
    source submit_analysis_${j}_Hcc_Summer23_noCTagCut.sh
    cd ..
    sleep 1
done


for j in Wqq_HT-200-400 Wqq_HT-400-600 Wqq_HT-600-800 Wqq_HT-800-inf
do
    echo -e "\nMC $j"
    echo -e "\npython createRunFile_new.py MC Hcc --outName Summer23_noCTagCut --n 30 --run 2023C --MCprocess $j "
    python createRunFile_new.py MC Hcc --outName Summer23_noCTagCut --n 50 --run 2023C --MCprocess $j
    sleep 1
    echo -e "\nsubmit_analysis_${j}_Hcc_Summer23_noCTagCut.sh"
    source submit_analysis_${j}_Hcc_Summer23_noCTagCut.sh
    cd ..
    sleep 1
done


echo -e "\npython createRunFile_new.py MC Hcc --outName Summer23_noCTagCut --n 30 --run 2023C --MCprocess TTto2L2Nu"
    python createRunFile_new.py MC Hcc --outName Summer23_noCTagCut --n 50 --run 2023C --MCprocess TTto2L2Nu
    sleep 1
    echo -e "\nsubmit_analysis_TTto2L2Nu_Hcc_Summer23_noCTagCut.sh"
    source submit_analysis_TTto2L2Nu_Hcc_Summer23_noCTagCut.sh
    cd ..
    sleep 1

echo -e "\npython createRunFile_new.py MC Hcc --outName Summer23_noCTagCut --n 30 --run 2023C --MCprocess VBFHCC"
    python createRunFile_new.py MC Hcc --outName Summer23_noCTagCut --n 50 --run 2023C --MCprocess VBFHCC
    sleep 1
    echo -e "\nsubmit_analysis_VBFHCC_Hcc_Summer23_noCTagCut.sh"
    source submit_analysis_VBFHCC_Hcc_Summer23_noCTagCut.sh
    cd ..
    sleep 1
