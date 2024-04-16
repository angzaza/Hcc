#!/bin/bash
# Usage:
#  submitAllJobs.sh



cd 2023C_Hcc_ReReco22Sep_v3Mar24
echo -e "hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_data_2023C_Hcc_merged.root AnalysedTree_data_2023C_Hcc*.root"
hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_data_2023C_Hcc_merged.root AnalysedTree_data_2023C_Hcc*.root
cd ..
sleep 1

for j in Zqq1j_pt-100-200 Zqq1j_pt-200-400 Zqq1j_pt-400-600 Zqq1j_pt-600
do 
		cd ${j}_Hcc_Summer23_v3Mar24
		echo -e "hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_${j}_Hcc_merged.root AnalysedTree_MC_${j}_Hcc*.root"
		hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_${j}_Hcc_merged.root AnalysedTree_MC_${j}_Hcc*.root
		cd ..
		sleep 1
done

for j in Wqq1j_pt-100-200 Wqq1j_pt-200-400 Wqq1j_pt-400-600 Wqq1j_pt-600
do 
		cd ${j}_Hcc_Summer23_v3Mar24
		echo -e "hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_${j}_Hcc_merged.root AnalysedTree_MC_${j}_Hcc*.root"
		hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_${j}_Hcc_merged.root AnalysedTree_MC_${j}_Hcc*.root
		cd ..
		sleep 1
done

for j in Zqq2j_pt-100-200 Zqq2j_pt-200-400 Zqq2j_pt-400-600 Zqq2j_pt-600
do 
		cd ${j}_Hcc_Summer23_v3Mar24
		echo -e "hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_${j}_Hcc_merged.root AnalysedTree_MC_${j}_Hcc*.root"
		hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_${j}_Hcc_merged.root AnalysedTree_MC_${j}_Hcc*.root
		cd ..
		sleep 1
done

for j in Wqq2j_pt-100-200 Wqq2j_pt-200-400 Wqq2j_pt-400-600 Wqq2j_pt-600
do 
		cd ${j}_Hcc_Summer23_v3Mar24
		echo -e "hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_${j}_Hcc_merged.root AnalysedTree_MC_${j}_Hcc*.root"
		hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_${j}_Hcc_merged.root AnalysedTree_MC_${j}_Hcc*.root
		cd ..
		sleep 1
done

for j in QCD-4Jets_HT-100to200 QCD-4Jets_HT-200to400 QCD-4Jets_HT-400to600 QCD-4Jets_HT-600to800  QCD-4Jets_HT-800to1000 QCD-4Jets_HT-1000to1200 QCD-4Jets_HT-1200to1500 QCD-4Jets_HT-1500to2000 QCD-4Jets_HT-2000
do
    cd ${j}_Hcc_Summer23_v3Mar24
    echo -e "hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_${j}_Hcc_merged.root AnalysedTree_MC_${j}_Hcc*.root"
    hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_${j}_Hcc_merged.root AnalysedTree_MC_${j}_Hcc*.root
    cd ..
    sleep 1
done

cd VBFZqq_Hcc_Summer23_v3Mar24
echo -e "hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_VBFZqq_Hcc_merged.root AnalysedTree_MC_VBFZqq_Hcc*.root"
hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_VBFZqq_Hcc_merged.root AnalysedTree_MC_VBFZqq_Hcc*.root
cd ..
sleep 1

cd VBFWqq_Hcc_Summer23_v3Mar24
echo -e "hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_VBFWqq_Hcc_merged.root AnalysedTree_MC_VBFWqq_Hcc*.root"
hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_VBFWqq_Hcc_merged.root AnalysedTree_MC_VBFWqq_Hcc*.root
cd ..
sleep 1

cd TTto2L2Nu_Hcc_Summer23_v3Mar24
echo -e "hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_TTto2L2Nu_Hcc_merged.root AnalysedTree_MC_TTto2L2Nu_Hcc*.root"
hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_TTto2L2Nu_Hcc_merged.root AnalysedTree_MC_TTto2L2Nu_Hcc*.root
cd ..
sleep 1

cd TTto4Q_Hcc_Summer23_v3Mar24
echo -e "hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_TTto4Q_Hcc_merged.root AnalysedTree_MC_TTto4Q_Hcc*.root"
hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_TTto4Q_Hcc_merged.root AnalysedTree_MC_TTto4Q_Hcc*.root
cd ..
sleep 1

cd TTtoLNu2Q_Hcc_Summer23_v3Mar24
echo -e "hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_TTtoLNu2Q_Hcc_merged.root AnalysedTree_MC_TTtoLNu2Q_Hcc*.root"
hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_TTtoLNu2Q_Hcc_merged.root AnalysedTree_MC_TTtoLNu2Q_Hcc*.root
cd ..
sleep 1

cd TbarWplusto4Q_Hcc_Summer23_v3Mar24
echo -e "hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_TbarWplusto4Q_Hcc_merged.root AnalysedTree_MC_TbarWplusto4Q_Hcc*.root"
hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_TbarWplusto4Q_Hcc_merged.root AnalysedTree_MC_TbarWplusto4Q_Hcc*.root
cd ..
sleep 1

cd TbarWplustoLNu2Q_Hcc_Summer23_v3Mar24
echo -e "hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_TbarWplustoLNu2Q_Hcc_merged.root AnalysedTree_MC_TbarWplustoLNu2Q_Hcc*.root"
hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_TbarWplustoLNu2Q_Hcc_merged.root AnalysedTree_MC_TbarWplustoLNu2Q_Hcc*.root
cd ..
sleep 1

cd TWminusto4Q_Hcc_Summer23_v3Mar24
echo -e "hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_TWminusto4Q_Hcc_merged.root AnalysedTree_MC_TWminusto4Q_Hcc*.root"
hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_TWminusto4Q_Hcc_merged.root AnalysedTree_MC_TWminusto4Q_Hcc*.root
cd ..
sleep 1

#cd TWminustoLNu2Q_Hcc_Summer23_v3Mar24
#echo -e "hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_TWminustoLNu2Q_Hcc_merged.root AnalysedTree_MC_TWminustoLNu2Q_Hcc*.root"
#hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_TWminustoLNu2Q_Hcc_merged.root AnalysedTree_MC_TWminustoLNu2Q_Hcc*.root
#cd ..
#sleep 1

cd VBFHCC_Hcc_Summer23_v3Mar24
echo -e "hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_VBFHCC_Hcc_merged.root AnalysedTree_MC_VBFHCC_Hcc*.root"
hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_VBFHCC_Hcc_merged.root AnalysedTree_MC_VBFHCC_Hcc*.root
cd ..
sleep 1

cd ggHCC_Hcc_Summer23_v3Mar24
echo -e "hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_ggHCC_Hcc_merged.root AnalysedTree_MC_ggHCC_Hcc*.root"
hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_ggHCC_Hcc_merged.root AnalysedTree_MC_ggHCC_Hcc*.root
cd ..
sleep 1

cd VBFHBB_Hcc_Summer23_v3Mar24
echo -e "hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_VBFHBB_Hcc_merged.root AnalysedTree_MC_VBFHBB_Hcc*.root"
hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_VBFHBB_Hcc_merged.root AnalysedTree_MC_VBFHBB_Hcc*.root
cd ..
sleep 1

cd ggHBB_Hcc_Summer23_v3Mar24
echo -e "hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_ggHBB_Hcc_merged.root AnalysedTree_MC_ggHBB_Hcc*.root"
hadd ../mergedFiles/Hcc_v3Mar24/AnalysedTree_MC_ggHBB_Hcc_merged.root AnalysedTree_MC_ggHBB_Hcc*.root
cd ..
sleep 1
