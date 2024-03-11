#!/bin/bash
# Usage:
#  submitAllJobs.sh


echo -e "\npython createRunFile_new.py data Hcc --run 2023C --n 30 --outName ReReco22Sep_v1Mar24 "
python createRunFile_new.py data Hcc --run 2023C --n 100 --outName ReReco22Sep_v1Mar24
sleep 1
pwd
echo -e "\nsubmit_analysis_2023C_Hcc_ReReco22Sep_v1Mar24.sh"
source submit_analysis_2023C_Hcc_ReReco22Sep_v1Mar24.sh
cd ..
sleep 1


#for i in QCD_PT-120-170 QCD_PT-170-300 QCD_PT-300-470 QCD_PT-470-600 QCD_PT-600-800 QCD_PT-800-1000 QCD_PT-1400-1800 QCD_PT-1800-2400 QCD_PT-2400-3200
#do
# 		echo -e "\nMC $i"
#    echo -e "\npython createRunFile_new.py MC Hcc --outName v1 --n 30 --EE preEE --MCprocess $i "
#    python createRunFile_new.py MC Hcc --outName v1 --n 30 --EE preEE --MCprocess $i
#    sleep 1
#    echo -e "\nsubmit_analysis_${i}_Hcc_v1Mar24.sh"
#    source submit_analysis_${i}_Hcc_v1Mar24.sh
#    cd ..
#    sleep 1
#done


for j in Zqq_pt-100-200 Zqq_pt-200-400 Zqq_pt-400-600 Zqq_pt-600
do
    echo -e "\nMC $j"
    echo -e "\npython createRunFile_new.py MC Hcc --outName Summer23_v1Mar24 --n 30 --run 2023C --MCprocess $j "
    python createRunFile_new.py MC Hcc --outName Summer23_v1Mar24 --n 50 --run 2023C --MCprocess $j
    sleep 1
    echo -e "\nsubmit_analysis_${j}_Hcc_Summer23_v1Mar24.sh"
    source submit_analysis_${j}_Hcc_Summer23_v1Mar24.sh
    cd ..
    sleep 1
done


for j in Wqq_pt-100-200 Wqq_pt-200-400 Wqq_pt-400-600 Wqq_pt-600
do
    echo -e "\nMC $j"
    echo -e "\npython createRunFile_new.py MC Hcc --outName Summer23_v1Mar24 --n 30 --run 2023C --MCprocess $j "
    python createRunFile_new.py MC Hcc --outName Summer23_v1Mar24 --n 50 --run 2023C --MCprocess $j
    sleep 1
    echo -e "\nsubmit_analysis_${j}_Hcc_Summer23_v1Mar24.sh"
    source submit_analysis_${j}_Hcc_Summer23_v1Mar24.sh
    cd ..
    sleep 1
done

for j in QCD-4Jets_HT-100to200 QCD-4Jets_HT-200to400 QCD-4Jets_HT-400to600 QCD-4Jets_HT-600to800  QCD-4Jets_HT-800to1000 QCD-4Jets_HT-1000to1200 QCD-4Jets_HT-1200to1500 QCD-4Jets_HT-1500to2000 QCD-4Jets_HT-2000
do
    echo -e "\nMC $j"
    echo -e "\npython createRunFile_new.py MC Hcc --outName Summer23_v1Mar24 --n 30 --run 2023C --MCprocess $j "
    python createRunFile_new.py MC Hcc --outName Summer23_v1Mar24 --n 50 --run 2023C --MCprocess $j
    sleep 1
    echo -e "\nsubmit_analysis_${j}_Hcc_Summer23_v1Mar24.sh"
    source submit_analysis_${j}_Hcc_Summer23_v1Mar24.sh
    cd ..
    sleep 1
done

echo -e "\npython createRunFile_new.py MC Hcc --outName Summer23_v1Mar24 --n 30 --run 2023C --MCprocess VBFZqq"
    python createRunFile_new.py MC Hcc --outName Summer23_v1Mar24 --n 50 --run 2023C --MCprocess VBFZqq
    sleep 1
    echo -e "\nsubmit_analysis_VBFZqq_Hcc_Summer23_v1Mar24.sh"
    source submit_analysis_VBFZqq_Hcc_Summer23_v1Mar24.sh
    cd ..
    sleep 1

echo -e "\npython createRunFile_new.py MC Hcc --outName Summer23_v1Mar24 --n 30 --run 2023C --MCprocess VBFWqq"
    python createRunFile_new.py MC Hcc --outName Summer23_v1Mar24 --n 50 --run 2023C --MCprocess VBFWqq
    sleep 1
    echo -e "\nsubmit_analysis_VBFWqq_Hcc_Summer23_v1Mar24.sh"
    source submit_analysis_VBFWqq_Hcc_Summer23_v1Mar24.sh
    cd ..
    sleep 1

echo -e "\npython createRunFile_new.py MC Hcc --outName Summer23_v1Mar24 --n 30 --run 2023C --MCprocess TTto2L2Nu"
    python createRunFile_new.py MC Hcc --outName Summer23_v1Mar24 --n 50 --run 2023C --MCprocess TTto2L2Nu
    sleep 1
    echo -e "\nsubmit_analysis_TTto2L2Nu_Hcc_Summer23_v1Mar24.sh"
    source submit_analysis_TTto2L2Nu_Hcc_Summer23_v1Mar24.sh
    cd ..
    sleep 1

echo -e "\npython createRunFile_new.py MC Hcc --outName Summer23_v1Mar24 --n 30 --run 2023C --MCprocess VBFHCC"
    python createRunFile_new.py MC Hcc --outName Summer23_v1Mar24 --n 50 --run 2023C --MCprocess VBFHCC
    sleep 1
    echo -e "\nsubmit_analysis_VBFHCC_Hcc_Summer23_v1Mar24.sh"
    source submit_analysis_VBFHCC_Hcc_Summer23_v1Mar24.sh
    cd ..
    sleep 1

echo -e "\npython createRunFile_new.py MC Hcc --outName Summer23_v1Mar24 --n 30 --run 2023C --MCprocess VBFHBB"
    python createRunFile_new.py MC Hcc --outName Summer23_v1Mar24 --n 50 --run 2023C --MCprocess VBFHBB
    sleep 1
    echo -e "\nsubmit_analysis_VBFHBB_Hcc_Summer23_v1Mar24.sh"
    source submit_analysis_VBFHBB_Hcc_Summer23_v1Mar24.sh
    cd ..
    sleep 1

echo -e "\npython createRunFile_new.py MC Hcc --outName Summer23_v1Mar24 --n 30 --run 2023C --MCprocess ggHCC"
    python createRunFile_new.py MC Hcc --outName Summer23_v1Mar24 --n 50 --run 2023C --MCprocess ggHCC
    sleep 1
    echo -e "\nsubmit_analysis_ggHCC_Hcc_Summer23_v1Mar24.sh"
    source submit_analysis_ggHCC_Hcc_Summer23_v1Mar24.sh
    cd ..
    sleep 1

echo -e "\npython createRunFile_new.py MC Hcc --outName Summer23_v1Mar24 --n 30 --run 2023C --MCprocess ggHBB"
    python createRunFile_new.py MC Hcc --outName Summer23_v1Mar24 --n 50 --run 2023C --MCprocess ggHBB
    sleep 1
    echo -e "\nsubmit_analysis_ggHBB_Hcc_Summer23_v1Mar24.sh"
    source submit_analysis_ggHBB_Hcc_Summer23_v1Mar24.sh
    cd ..
    sleep 1
