import ROOT


df = ROOT.RDataFrame("FinalTree", "2023C_Hcc_v3/AnalysedTree_data_2023C_Hcc_merged_v3.root")
df.Filter("entry % 20 ==0").Snapshot("FinalTree", "2023C_Hcc_v3/AnalysedTree_forTranining_data_2023C_Hcc_v3.root",)
df.Filter("entry % 20 !=0").Snapshot("FinalTree", "2023C_Hcc_v3/AnalysedTree_forAnalysis_data_2023C_Hcc_v3.root",)

