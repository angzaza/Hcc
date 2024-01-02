#!/bin/bash
# Usage:
#  submitAllJobs.sh



cd 2023C_Hcc_ReReco22Sep_noCTagCut
echo -e "hadd AnalysedTree_data_2023C_Hcc_merged.root AnalysedTree_data_2023C_Hcc*.root"
hadd AnalysedTree_data_2023C_Hcc_merged.root AnalysedTree_data_2023C_Hcc*.root
cd ..
sleep 1

for j in Zqq_HT-200-400 Zqq_HT-400-600 Zqq_HT-600-800 Zqq_HT-800-inf
do 
		cd ${j}_Hcc_Summer23_noCTagCut
		echo -e "hadd AnalysedTree_MC_${j}_Hcc_merged.root AnalysedTree_MC_${j}_Hcc*.root"
		hadd AnalysedTree_MC_${j}_Hcc_merged.root AnalysedTree_MC_${j}_Hcc*.root
		cd ..
		sleep 1
done

for j in Wqq_HT-200-400 Wqq_HT-400-600 Wqq_HT-600-800 Wqq_HT-800-inf
do 
		cd ${j}_Hcc_Summer23_noCTagCut
		echo -e "hadd AnalysedTree_MC_${j}_Hcc_merged.root AnalysedTree_MC_${j}_Hcc*.root"
		hadd AnalysedTree_MC_${j}_Hcc_merged.root AnalysedTree_MC_${j}_Hcc*.root
		cd ..
		sleep 1
done

cd TTto2L2Nu_Hcc_Summer23_noCTagCut
echo -e "hadd AnalysedTree_MC_TTto2L2Nu_Hcc_merged.root AnalysedTree_MC_TTto2L2Nu_Hcc*.root"
hadd AnalysedTree_MC_TTto2L2Nu_Hcc_merged.root AnalysedTree_MC_TTto2L2Nu_Hcc*.root
cd ..
sleep 1

cd VBFHCC_Hcc_Summer23_noCTagCut
echo -e "hadd AnalysedTree_MC_VBFHCC_Hcc_merged.root AnalysedTree_MC_VBFHCC_Hcc*.root"
hadd AnalysedTree_MC_VBFHCC_Hcc_merged.root AnalysedTree_MC_VBFHCC_Hcc*.root
cd ..
sleep 1
