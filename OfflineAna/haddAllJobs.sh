#!/bin/bash
# Usage:
#  submitAllJobs.sh



cd 2023C_Hcc_ReReco22Sep_v1Mar24
echo -e "hadd ../mergedFiles/Hcc_v1Mar24/AnalysedTree_data_2023C_Hcc_merged.root AnalysedTree_data_2023C_Hcc*.root"
hadd ../mergedFiles/Hcc_v1Mar24/AnalysedTree_data_2023C_Hcc_merged.root AnalysedTree_data_2023C_Hcc*.root
cd ..
sleep 1

for j in Zqq_pt-100-200 Zqq_pt-200-400 Zqq_pt-400-600 Zqq_pt-600
do 
		cd ${j}_Hcc_Summer23_v1Mar24
		echo -e "hadd ../mergedFiles/Hcc_v1Mar24/AnalysedTree_MC_${j}_Hcc_merged.root AnalysedTree_MC_${j}_Hcc*.root"
		hadd ../mergedFiles/Hcc_v1Mar24/AnalysedTree_MC_${j}_Hcc_merged.root AnalysedTree_MC_${j}_Hcc*.root
		cd ..
		sleep 1
done

for j in Wqq_pt-100-200 Wqq_pt-200-400 Wqq_pt-400-600 Wqq_pt-600
do 
		cd ${j}_Hcc_Summer23_v1Mar24
		echo -e "hadd ../mergedFiles/Hcc_v1Mar24/AnalysedTree_MC_${j}_Hcc_merged.root AnalysedTree_MC_${j}_Hcc*.root"
		hadd ../mergedFiles/Hcc_v1Mar24/AnalysedTree_MC_${j}_Hcc_merged.root AnalysedTree_MC_${j}_Hcc*.root
		cd ..
		sleep 1
done

for j in QCD-4Jets_HT-100to200 QCD-4Jets_HT-200to400 QCD-4Jets_HT-400to600 QCD-4Jets_HT-600to800  QCD-4Jets_HT-800to1000 QCD-4Jets_HT-1000to1200 QCD-4Jets_HT-1200to1500 QCD-4Jets_HT-1500to2000 QCD-4Jets_HT-2000
do
    cd ${j}_Hcc_Summer23_v1Mar24
    echo -e "hadd ../mergedFiles/Hcc_v1Mar24/AnalysedTree_MC_${j}_Hcc_merged.root AnalysedTree_MC_${j}_Hcc*.root"
    hadd ../mergedFiles/Hcc_v1Mar24/AnalysedTree_MC_${j}_Hcc_merged.root AnalysedTree_MC_${j}_Hcc*.root
    cd ..
    sleep 1
done

cd VBFZqq_Hcc_Summer23_v1Mar24
echo -e "hadd ../mergedFiles/Hcc_v1Mar24/AnalysedTree_MC_VBFZqq_Hcc_merged.root AnalysedTree_MC_VBFZqq_Hcc*.root"
hadd ../mergedFiles/Hcc_v1Mar24/AnalysedTree_MC_VBFZqq_Hcc_merged.root AnalysedTree_MC_VBFZqq_Hcc*.root
cd ..
sleep 1

cd VBFWqq_Hcc_Summer23_v1Mar24
echo -e "hadd ../mergedFiles/Hcc_v1Mar24/AnalysedTree_MC_VBFWqq_Hcc_merged.root AnalysedTree_MC_VBFWqq_Hcc*.root"
hadd ../mergedFiles/Hcc_v1Mar24/AnalysedTree_MC_VBFWqq_Hcc_merged.root AnalysedTree_MC_VBFWqq_Hcc*.root
cd ..
sleep 1

cd TTto2L2Nu_Hcc_Summer23_v1Mar24
echo -e "hadd ../mergedFiles/Hcc_v1Mar24/AnalysedTree_MC_TTto2L2Nu_Hcc_merged.root AnalysedTree_MC_TTto2L2Nu_Hcc*.root"
hadd ../mergedFiles/Hcc_v1Mar24/AnalysedTree_MC_TTto2L2Nu_Hcc_merged.root AnalysedTree_MC_TTto2L2Nu_Hcc*.root
cd ..
sleep 1

cd VBFHCC_Hcc_Summer23_v1Mar24
echo -e "hadd ../mergedFiles/Hcc_v1Mar24/AnalysedTree_MC_VBFHCC_Hcc_merged.root AnalysedTree_MC_VBFHCC_Hcc*.root"
hadd ../mergedFiles/Hcc_v1Mar24/AnalysedTree_MC_VBFHCC_Hcc_merged.root AnalysedTree_MC_VBFHCC_Hcc*.root
cd ..
sleep 1

cd ggHCC_Hcc_Summer23_v1Mar24
echo -e "hadd ../mergedFiles/Hcc_v1Mar24/AnalysedTree_MC_ggHCC_Hcc_merged.root AnalysedTree_MC_ggHCC_Hcc*.root"
hadd ../mergedFiles/Hcc_v1Mar24/AnalysedTree_MC_ggHCC_Hcc_merged.root AnalysedTree_MC_ggHCC_Hcc*.root
cd ..
sleep 1

cd VBFHBB_Hcc_Summer23_v1Mar24
echo -e "hadd ../mergedFiles/Hcc_v1Mar24/AnalysedTree_MC_VBFHBB_Hcc_merged.root AnalysedTree_MC_VBFHBB_Hcc*.root"
hadd ../mergedFiles/Hcc_v1Mar24/AnalysedTree_MC_VBFHBB_Hcc_merged.root AnalysedTree_MC_VBFHBB_Hcc*.root
cd ..
sleep 1

cd ggHBB_Hcc_Summer23_v1Mar24
echo -e "hadd ../mergedFiles/Hcc_v1Mar24/AnalysedTree_MC_ggHBB_Hcc_merged.root AnalysedTree_MC_ggHBB_Hcc*.root"
hadd ../mergedFiles/Hcc_v1Mar24/AnalysedTree_MC_ggHBB_Hcc_merged.root AnalysedTree_MC_ggHBB_Hcc*.root
cd ..
sleep 1
