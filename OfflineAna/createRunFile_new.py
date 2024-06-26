import sys
import os
import csv
import string
import datetime

# Define the parser
import argparse
parser = argparse.ArgumentParser(description="Options to give to the script")
# Positional arguments
parser.add_argument("dataset", type=str, choices=['data', 'data_control','data_validate', 'MC', 'MC_control', 'MC_validate'], help="Specify if data or Monte Carlo")
parser.add_argument("anatype", type=str, choices=['Hcc', 'control', 'validate'], help="Specify analysis type")
#parser.add_argument("--run", type=str, default='', choices=['2022B', '2022C_0', '2022C_1', '2022C_2', '2022C_3', '2022C_4', '2022C_5', '2022C_6', '2022C_7', '2022D-v1_0', '2022D-v1_1', '2022D-v1_2', '2022D-v1_3', '2022D-v1_4', '2022D-v1_5', '2022D-v1_6', '2022D-v1_7', '2022D-v2_0', '2022D-v2_1', '2022D-v2_2', '2022D-v2_3', '2022D-v2_4', '2022D-v2_5', '2022D-v2_6', '2022D-v2_7', '2022E_0', '2022E_1', '2022E_2', '2022E_3', '2022E_4', '2022E_5', '2022E_6', '2022E_7', '2022F_0', '2022F_1', '2022F_2', '2022F_3', '2022F_4', '2022F_5', '2022F_6', '2022F_7', '2022G_0', '2022G_1', '2022G_2', '2022G_3', '2022G_4', '2022G_5', '2022G_6', '2022G_7'], help="run in data")
parser.add_argument("--run", type=str, default='', choices=['2023D', '2023C', '2022B', '2022C', '2022D', '2022D', '2022E', '2022F', '2022G'], help="run in data")
# Optional Arguments
parser.add_argument("--outName", type=str, default="test", help="Specify name for output files")
parser.add_argument("--n", type=int, default=255, help="number of .root files per job")
#parser.add_argument("--EE",type=str, default='postEE', choices=['preEE','postEE'], help="specify if it is simulated with preEE or postEE conditions")
parser.add_argument("--MCprocess", type=str, default='', choices=['VBFHCC', 'VBFHBB', 'ggHCC','ggHBB', 'QCD','VBFZqq', 'VBFWqq','Zqq1j_pt-100-200','Zqq1j_pt-200-400','Zqq1j_pt-400-600','Zqq1j_pt-600','Zqq2j_pt-100-200','Zqq2j_pt-200-400','Zqq2j_pt-400-600','Zqq2j_pt-600', 'Wqq1j_pt-100-200','Wqq1j_pt-200-400','Wqq1j_pt-400-600','Wqq1j_pt-600', 'Wqq2j_pt-100-200','Wqq2j_pt-200-400','Wqq2j_pt-400-600','Wqq2j_pt-600', 'Zqq_HT-200-400', 'Zqq_HT-400-600', 'Zqq_HT-600-800', 'Zqq_HT-800-inf', 'Wqq_HT-200-400', 'Wqq_HT-400-600', 'Wqq_HT-600-800', 'Wqq_HT-800-inf', 'TTto2L2Nu','TTto4Q', 'TTtoLNu2Q', 'TWminusto4Q', 'TbarWplusto4Q', 'TbarWplustoLNu2Q', 'TWminustoLNu2Q', 'QCD-4Jets_HT-100to200', 'QCD-4Jets_HT-200to400', 'QCD-4Jets_HT-400to600', 'QCD-4Jets_HT-600to800', 'QCD-4Jets_HT-800to1000','QCD-4Jets_HT-1000to1200',  'QCD-4Jets_HT-1200to1500', 'QCD-4Jets_HT-1500to2000', 'QCD-4Jets_HT-2000', 'QCD_PT-120-170', 'QCD_PT-170-300', 'QCD_PT-300-470', 'QCD_PT-470-600', 'QCD_PT-600-800', 'QCD_PT-800-1000', 'QCD_PT-1000-1400', 'QCD_PT-1400-1800', 'QCD_PT-1800-2400', 'QCD_PT-2400-3200'], help="process in Monte Carlo")
args = parser.parse_args()

#prepare output filename  and option string
if args.dataset == 'data':
   out_filename = 'AnalysedTree_'+args.dataset+'_'+args.run+'_'+args.anatype
   temp = '_'+args.anatype
   option_string = ' "'+args.dataset+temp.replace("_Hcc","")+'" "'+args.run+'" "'+args.run+'"'
elif args.dataset == 'data_control':
   out_filename = 'AnalysedTree_'+args.dataset+'_'+args.run+'_'+args.anatype
   temp = '_'+args.anatype
   option_string = ' "'+args.dataset+temp.replace("_control","")+'" "'+args.run+'" "'+args.run+'"'
elif args.dataset == 'data_validate' :
   out_filename = 'AnalysedTree_'+args.dataset+'_'+args.run+'_'+args.anatype
   temp = '_'+args.anatype
   option_string = ' "'+args.dataset+temp.replace("_validate","")+'" "'+args.run+'" "'+args.run+'"'
elif args.dataset == 'MC':
   out_filename = 'AnalysedTree_'+args.dataset+'_'+args.MCprocess+'_'+args.anatype
   temp = '_'+args.anatype
   option_string = ' "'+args.dataset+temp.replace("_Hcc","")+'" "'+args.MCprocess+'" "'+args.run+'"'
elif args.dataset == 'MC_control':
   out_filename = 'AnalysedTree_'+args.dataset+'_'+args.MCprocess+'_'+args.anatype
   temp = '_'+args.anatype
   option_string = ' "'+args.dataset+temp.replace("_control","")+'" "'+args.MCprocess+'" "'+args.run+'"'
elif args.dataset == 'MC_validate':
   out_filename = 'AnalysedTree_'+args.dataset+'_'+args.MCprocess+'_'+args.anatype
   temp = '_'+args.anatype
   option_string = ' "'+args.dataset+temp.replace("_validate","")+'" "'+args.MCprocess+'" "'+args.run+'"'

#startTime = datetime.datetime.now().strftime("%Y%m%d_%H%M")

# Create target Directory if don't exist
if args.dataset == 'MC' or args.dataset == 'MC_control' or args.dataset == 'MC_validate':
   output_name = args.MCprocess+"_"+args.anatype+"_"+args.outName
   #output_name = args.MCprocess+"_"+args.anatype+"_"+args.EE+"_"+args.outName
else: 
   output_name = args.run+"_"+args.anatype+"_"+args.outName

if not os.path.exists(output_name):
    os.mkdir(output_name)
    print('Directory '+output_name+' created\n')
else:    
    print('Directory '+output_name+' already exists\n')

if args.anatype == 'Hcc':
   #### 2022
   if args.dataset == 'data' and args.run == '2022B':
      path = '' 
   if args.dataset == 'data' and args.run == '2022C':
      path = '/lustre/cms/store/user/dtroiano/JetMET/Run2022C_ntuple/230405_142008' 
   if args.dataset == 'data' and args.run == '2022D':
      path = '/lustre/cms/store/user/dtroiano/JetMET/Run2022D_ntuple/230405_143608'
   if args.dataset == 'data' and args.run == '2022E':
      path = '/lustre/cms/store/user/dtroiano/JetMET/Run2022E_ntuple/230405_144208'
   if args.dataset == 'data' and args.run == '2022F':
      path = '/lustre/cms/store/user/azaza/JetMET/DataRun3_EraF_ntuple/230405_155207'
   if args.dataset == 'data' and args.run == '2022G':
      path = '/lustre/cms/store/user/azaza/JetMET/DataRun3_EraG_ntuple/230405_155244'
   if args.dataset == 'data' and args.run == '2023C':
      path = '/lustre/cms/store/user/azaza/Data2023C_Prompt_19Mar24'
   if args.dataset == 'data' and args.run == '2023D':
      path = '/lustre/cms/store/user/azaza/Data2023_CD_Oct23_v1/ERAD'

if args.anatype == 'control':
   if args.dataset == 'data_control' and args.run == '2023C':
   #if args.dataset == 'data_control':
      path = '/lustre/cms/store/user/azaza/Data2023C_Prompt_19Mar24'
   if args.dataset == 'MC_control':
      path = '/lustre/cms/store/user/azaza/QCD_Summer23_HT'

if args.anatype == 'validate':
   if args.dataset == 'data_validate' and args.run == '2023C':
   #if args.dataset == 'data_control':
      path = '/lustre/cms/store/user/azaza/Data2023C_Prompt_19Mar24'
   #if args.dataset == 'MC_validate':
   #   path = '/lustre/cms/store/user/azaza/QCD_Summer23_HT'

if args.run == '2023C':
   if args.dataset == 'MC' and args.MCprocess == 'VBFHCC':
      path =  '/lustre/cms/store/user/azaza/VBFHto2C_M-125_TuneCP5_13p6TeV_powheg-pythia8/crab_VBFHto2C_M-125_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer23MiniAODv4/240319_144733'
   if args.dataset == 'MC' and args.MCprocess == 'VBFHBB':
      path = '/lustre/cms/store/user/azaza/VBFHto2B_M-125_TuneCP5_13p6TeV_powheg-pythia8/crab_VBFHto2B_M-125_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer23MiniAODv4/240319_144816/'
   if args.dataset == 'MC' and args.MCprocess == 'ggHCC':
      path = '/lustre/cms/store/user/azaza/GluGluHto2C_M-125_TuneCP5_13p6TeV_powheg-minlo-pythia8/crab_GluGluHto2C_M-125_TuneCP5_13p6TeV_powheg-minlo-pythia8_Run3Summer23MiniAODv4/240319_144900/'
   if args.dataset == 'MC' and args.MCprocess == 'ggHBB':
      path = '/lustre/cms/store/user/azaza/GluGluHto2B_M-125_TuneCP5_13p6TeV_powheg-minlo-pythia8/crab_GluGluHto2B_M-125_TuneCP5_13p6TeV_powheg-minlo-pythia8_Run3Summer23MiniAODv4/240319_144943/'
   if (args.dataset == 'MC' or args.dataset == 'MC_validate') and args.MCprocess == 'QCD-4Jets_HT-100to200':
      path = '/lustre/cms/store/user/azaza/QCD_Summer23_HT/QCD-4Jets_HT-100to200_TuneCP5_13p6TeV_madgraphMLM-pythia8/crab_QCD-4Jets_HT-100to200_TuneCP5_13p6TeV_madgraphMLM-pythia8_Run3Summer23MiniAODv4/240319_150226'
   if (args.dataset == 'MC' or args.dataset == 'MC_validate') and args.MCprocess == 'QCD-4Jets_HT-200to400':
      path = '/lustre/cms/store/user/azaza/QCD_Summer23_HT/QCD-4Jets_HT-200to400_TuneCP5_13p6TeV_madgraphMLM-pythia8/crab_QCD-4Jets_HT-200to400_TuneCP5_13p6TeV_madgraphMLM-pythia8_Run3Summer23MiniAODv4/240319_150307'
   if (args.dataset == 'MC' or args.dataset == 'MC_validate') and args.MCprocess == 'QCD-4Jets_HT-400to600':
      path = '/lustre/cms/store/user/azaza/QCD_Summer23_HT/QCD-4Jets_HT-400to600_TuneCP5_13p6TeV_madgraphMLM-pythia8/crab_QCD-4Jets_HT-400to600_TuneCP5_13p6TeV_madgraphMLM-pythia8_Run3Summer23MiniAODv4/240319_150348'
   if (args.dataset == 'MC' or args.dataset == 'MC_validate') and args.MCprocess == 'QCD-4Jets_HT-600to800':
      path = '/lustre/cms/store/user/azaza/QCD_Summer23_HT/QCD-4Jets_HT-600to800_TuneCP5_13p6TeV_madgraphMLM-pythia8/crab_QCD-4Jets_HT-600to800_TuneCP5_13p6TeV_madgraphMLM-pythia8_Run3Summer23MiniAODv4/240319_150430'
   if (args.dataset == 'MC' or args.dataset == 'MC_validate') and args.MCprocess == 'QCD-4Jets_HT-800to1000':
      path = '/lustre/cms/store/user/azaza/QCD_Summer23_HT/QCD-4Jets_HT-800to1000_TuneCP5_13p6TeV_madgraphMLM-pythia8/crab_QCD-4Jets_HT-800to1000_TuneCP5_13p6TeV_madgraphMLM-pythia8_Run3Summer23MiniAODv4/240319_150511'
   if (args.dataset == 'MC' or args.dataset == 'MC_validate') and args.MCprocess == 'QCD-4Jets_HT-1000to1200':
      path = '/lustre/cms/store/user/azaza/QCD_Summer23_HT/QCD-4Jets_HT-1000to1200_TuneCP5_13p6TeV_madgraphMLM-pythia8/crab_QCD-4Jets_HT-1000to1200_TuneCP5_13p6TeV_madgraphMLM-pythia8_Run3Summer23MiniAODv4/240319_150551'
   if (args.dataset == 'MC' or args.dataset == 'MC_validate') and args.MCprocess == 'QCD-4Jets_HT-1200to1500':
      path = '/lustre/cms/store/user/azaza/QCD_Summer23_HT/QCD-4Jets_HT-1200to1500_TuneCP5_13p6TeV_madgraphMLM-pythia8/crab_QCD-4Jets_HT-1200to1500_TuneCP5_13p6TeV_madgraphMLM-pythia8_Run3Summer23MiniAODv4/240319_150635'
   if (args.dataset == 'MC' or args.dataset == 'MC_validate') and args.MCprocess == 'QCD-4Jets_HT-1500to2000':
      path = '/lustre/cms/store/user/azaza/QCD_Summer23_HT/QCD-4Jets_HT-1500to2000_TuneCP5_13p6TeV_madgraphMLM-pythia8/crab_QCD-4Jets_HT-1500to2000_TuneCP5_13p6TeV_madgraphMLM-pythia8_Run3Summer23MiniAODv4/240319_150715'
   if (args.dataset == 'MC' or args.dataset == 'MC_validate') and args.MCprocess == 'QCD-4Jets_HT-2000':
      path = '/lustre/cms/store/user/azaza/QCD_Summer23_HT/QCD-4Jets_HT-2000_TuneCP5_13p6TeV_madgraphMLM-pythia8/crab_QCD-4Jets_HT-2000_TuneCP5_13p6TeV_madgraphMLM-pythia8_Run3Summer23MiniAODv4/240319_150755'
   if args.dataset == 'MC' and args.MCprocess == 'Zqq_HT-200-400':
      path = '/lustre/cms/store/user/azaza/Zto2Q-4Jets_HT-200to400_TuneCP5_13p6TeV_madgraphMLM-pythia8/crab_Zto2Q-4Jets_HT-200to400_TuneCP5_13p6TeV_madgraphMLM-pythia8_Run3Summer23MiniAODv4/'
   if args.dataset == 'MC' and args.MCprocess == 'Zqq_HT-400-600':
      path = '/lustre/cms/store/user/azaza/Zto2Q-4Jets_HT-400to600_TuneCP5_13p6TeV_madgraphMLM-pythia8/Zto2Q-4Jets_HT-400to600_TuneCP5_13p6TeV_madgraphMLM-pythia8/crab_Zto2Q-4Jets_HT-400to600_TuneCP5_13p6TeV_madgraphMLM-pythia8_Run3Summer23MiniAODv4/240118_170013/' 
   if args.dataset == 'MC' and args.MCprocess == 'Zqq_HT-600-800':
      path = '/lustre/cms/store/user/azaza/Zto2Q-4Jets_HT-600to800_TuneCP5_13p6TeV_madgraphMLM-pythia8/crab_Zto2Q-4Jets_HT-600to800_TuneCP5_13p6TeV_madgraphMLM-pythia8_Run3Summer23MiniAODv4/240118_170039/' 
   if args.dataset == 'MC' and args.MCprocess == 'Zqq_HT-800-inf':
      path = '/lustre/cms/store/user/azaza/Zto2Q-4Jets_HT-800_TuneCP5_13p6TeV_madgraphMLM-pythia8/crab_Zto2Q-4Jets_HT-800_TuneCP5_13p6TeV_madgraphMLM-pythia8_Run3Summer23MiniAODv4/240118_170103/' 
   if args.dataset == 'MC' and args.MCprocess == 'Wqq_HT-200-400':
      path = '/lustre/cms/store/user/azaza/Wto2Q-3Jets_HT-200to400_TuneCP5_13p6TeV_madgraphMLM-pythia8/crab_Wto2Q-3Jets_HT-200to400_TuneCP5_13p6TeV_madgraphMLM-pythia8_Run3Summer23MiniAODv4/240118_170150/' 
   if args.dataset == 'MC' and args.MCprocess == 'Wqq_HT-400-600':
      path = '/lustre/cms/store/user/azaza/Wto2Q-3Jets_HT-400to600_TuneCP5_13p6TeV_madgraphMLM-pythia8/crab_Wto2Q-3Jets_HT-400to600_TuneCP5_13p6TeV_madgraphMLM-pythia8_Run3Summer23MiniAODv4/240118_170214/' 
   if args.dataset == 'MC' and args.MCprocess == 'Wqq_HT-600-800':
      path = '/lustre/cms/store/user/azaza/Wto2Q-3Jets_HT-600to800_TuneCP5_13p6TeV_madgraphMLM-pythia8/crab_Wto2Q-3Jets_HT-600to800_TuneCP5_13p6TeV_madgraphMLM-pythia8_Run3Summer23MiniAODv4/240118_170238/' 
   if args.dataset == 'MC' and args.MCprocess == 'Wqq_HT-800-inf':
      path = '/lustre/cms/store/user/azaza/Wto2Q-3Jets_HT-800_TuneCP5_13p6TeV_madgraphMLM-pythia8/crab_Wto2Q-3Jets_HT-800_TuneCP5_13p6TeV_madgraphMLM-pythia8_Run3Summer23MiniAODv4/240118_170302/'
   if args.dataset == 'MC' and args.MCprocess == 'Zqq1j_pt-100-200':
      path = '/lustre/cms/store/user/azaza/Zto2Q-2Jets_PTQQ-100to200_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/crab_Zto2Q-2Jets_PTQQ-100to200_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer23MiniAODv4/240319_145028/'
   if args.dataset == 'MC' and args.MCprocess == 'Zqq1j_pt-200-400':
      path = '/lustre/cms/store/user/azaza/Zto2Q-2Jets_PTQQ-200to400_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/crab_Zto2Q-2Jets_PTQQ-200to400_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer23MiniAODv4/240319_145108/'
   if args.dataset == 'MC' and args.MCprocess == 'Zqq1j_pt-400-600':
      path = '/lustre/cms/store/user/azaza/Zto2Q-2Jets_PTQQ-400to600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/crab_Zto2Q-2Jets_PTQQ-400to600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer23MiniAODv4/240319_145149/'
   if args.dataset == 'MC' and args.MCprocess == 'Zqq1j_pt-600':
      path = '/lustre/cms/store/user/azaza/Zto2Q-2Jets_PTQQ-600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/crab_Zto2Q-2Jets_PTQQ-600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer23MiniAODv4/240319_145229/'
   if args.dataset == 'MC' and args.MCprocess == 'Wqq1j_pt-100-200':
      path = '/lustre/cms/store/user/azaza/Wto2Q-2Jets_PTQQ-100to200_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/crab_Wto2Q-2Jets_PTQQ-100to200_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer23MiniAODv4/240319_145602/'
   if args.dataset == 'MC' and args.MCprocess == 'Wqq1j_pt-200-400':
      path = '/lustre/cms/store/user/azaza/Wto2Q-2Jets_PTQQ-200to400_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/crab_Wto2Q-2Jets_PTQQ-200to400_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer23MiniAODv4/240319_145643/'
   if args.dataset == 'MC' and args.MCprocess == 'Wqq1j_pt-400-600':
      path = '/lustre/cms/store/user/azaza/Wto2Q-2Jets_PTQQ-400to600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/crab_Wto2Q-2Jets_PTQQ-400to600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer23MiniAODv4/240319_145724/'
   if args.dataset == 'MC' and args.MCprocess == 'Wqq1j_pt-600':
      path = '/lustre/cms/store/user/azaza/Wto2Q-2Jets_PTQQ-600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/crab_Wto2Q-2Jets_PTQQ-600_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer23MiniAODv4/240313_170402/'
   if args.dataset == 'MC' and args.MCprocess == 'Zqq2j_pt-100-200':
      path = '/lustre/cms/store/user/azaza/Zto2Q-2Jets_PTQQ-100to200_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/crab_Zto2Q-2Jets_PTQQ-100to200_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer23MiniAODv4/240320_150456/'
   if args.dataset == 'MC' and args.MCprocess == 'Zqq2j_pt-200-400':
      path = '/lustre/cms/store/user/azaza/Zto2Q-2Jets_PTQQ-200to400_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/crab_Zto2Q-2Jets_PTQQ-200to400_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer23MiniAODv4/240320_150541/'
   if args.dataset == 'MC' and args.MCprocess == 'Zqq2j_pt-400-600':
      path = '/lustre/cms/store/user/azaza/Zto2Q-2Jets_PTQQ-400to600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/crab_Zto2Q-2Jets_PTQQ-400to600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer23MiniAODv4/240320_150623/'
   if args.dataset == 'MC' and args.MCprocess == 'Zqq2j_pt-600':
      path = '/lustre/cms/store/user/azaza/Zto2Q-2Jets_PTQQ-600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/crab_Zto2Q-2Jets_PTQQ-600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer23MiniAODv4/240320_150706/'
   if args.dataset == 'MC' and args.MCprocess == 'Wqq2j_pt-100-200':
      path = '/lustre/cms/store/user/azaza/Wto2Q-2Jets_PTQQ-100to200_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/crab_Wto2Q-2Jets_PTQQ-100to200_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer23MiniAODv4/240320_150751/'
   if args.dataset == 'MC' and args.MCprocess == 'Wqq2j_pt-200-400':
      path = '/lustre/cms/store/user/azaza/Wto2Q-2Jets_PTQQ-200to400_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/crab_Wto2Q-2Jets_PTQQ-200to400_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer23MiniAODv4/240320_150832/'
   if args.dataset == 'MC' and args.MCprocess == 'Wqq2j_pt-400-600':
      path = '/lustre/cms/store/user/azaza/Wto2Q-2Jets_PTQQ-400to600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/crab_Wto2Q-2Jets_PTQQ-400to600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer23MiniAODv4/240320_150914/'
   if args.dataset == 'MC' and args.MCprocess == 'Wqq2j_pt-600':
      path = '/lustre/cms/store/user/azaza/Wto2Q-2Jets_PTQQ-600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/crab_Wto2Q-2Jets_PTQQ-600_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_Run3Summer23MiniAODv4/240320_150958/'
   if args.dataset == 'MC' and args.MCprocess == 'VBFZqq':
      path = '/lustre/cms/store/user/azaza/VBFZto2Q_TuneCP5_13p6TeV_madgraph-pythia8/crab_VBFZto2Q_TuneCP5_13p6TeV_madgraph-pythia8_Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v3/240319_145309/'
   if args.dataset == 'MC' and args.MCprocess == 'VBFWqq':
      path = '/lustre/cms/store/user/azaza/VBFWto2Q_TuneCP5_13p6TeV_madgraph-pythia8/crab_VBFWto2Q_TuneCP5_13p6TeV_madgraph-pythia8_Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v3/240319_145846/'
   if args.dataset == 'MC' and args.MCprocess == 'TTto2L2Nu':
      path = '/lustre/cms/store/user/azaza/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/crab_TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v2/240319_145350/'
   if args.dataset == 'MC' and args.MCprocess == 'TTto4Q':
      path = '/lustre/cms/store/user/azaza/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8/crab_TTto4Q_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v2/240319_145509/'
   if args.dataset == 'MC' and args.MCprocess == 'TTtoLNu2Q':
      path = '/lustre/cms/store/user/azaza/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/crab_TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v2/240319_145429/'
   if args.dataset == 'MC' and args.MCprocess == 'TWminusto4Q':
     path = '/lustre/cms/store/user/azaza/TWminusto4Q_TuneCP5_13p6TeV_powheg-pythia8/crab_TWminusto4Q_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v2/240319_150008/'
   if args.dataset == 'MC' and args.MCprocess == 'TbarWplusto4Q':
     path = '/lustre/cms/store/user/azaza/TbarWplusto4Q_TuneCP5_13p6TeV_powheg-pythia8/crab_TbarWplusto4Q_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer23MiniAODv4/240319_145927/'
   if args.dataset == 'MC' and args.MCprocess == 'TbarWplustoLNu2Q':
     path = '/lustre/cms/store/user/azaza/TbarWplustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/crab_TbarWplustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8_Run3Summer23MiniAODv4/240319_150049/'

  

'''if args.EE == 'preEE':
   if args.dataset == 'MC' and args.MCprocess == 'QCD_PT-120-170':
      path = '/lustre/cms/store/user/azaza/QCD_PT-120to170_TuneCP5_13p6TeV_pythia8/crab_QCD_PT-120to170_TuneCP5_13p6TeV_pythia8_Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2/230405_161508'
      #path = '/lustre/cms/store/user/azaza/QCD_PT-120to170_PROVA'
   if args.dataset == 'MC' and args.MCprocess == 'QCD_PT-170-300':
      path = '/lustre/cms/store/user/azaza/QCD_PT-170to300_EMEnriched_TuneCP5_13p6TeV_pythia8/crab_QCD_PT-170to300_EMEnriched_TuneCP5_13p6TeV_pythia8_Run3Summer22MiniAODv3/230405_184147'
   if args.dataset == 'MC' and args.MCprocess == 'QCD_PT-300-470':
      path = '/lustre/cms/store/user/azaza/QCD_PT-300to470_TuneCP5_13p6TeV_pythia8/crab_QCD_PT-300to470_TuneCP5_13p6TeV_pythia8_Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2/230405_184210'
   if args.dataset == 'MC' and args.MCprocess == 'QCD_PT-470-600':
      path = '/lustre/cms/store/user/azaza/QCD_PT-470to600_TuneCP5_13p6TeV_pythia8/crab_QCD_PT-470to600_TuneCP5_13p6TeV_pythia8_Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2/230405_184232'
   if args.dataset == 'MC' and args.MCprocess == 'QCD_PT-600-800':
      path = '/lustre/cms/store/user/azaza/QCD_PT-600to800_TuneCP5_13p6TeV_pythia8/crab_QCD_PT-600to800_TuneCP5_13p6TeV_pythia8_Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2/230405_184254'
   if args.dataset == 'MC' and args.MCprocess == 'QCD_PT-800-1000':
      path = '/lustre/cms/store/user/azaza/QCD_PT-800to1000_TuneCP5_13p6TeV_pythia8/crab_QCD_PT-800to1000_TuneCP5_13p6TeV_pythia8_Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2/230405_184317'
   if args.dataset == 'MC' and args.MCprocess == 'QCD_PT-1000-1400':
      path = ''
   if args.dataset == 'MC' and args.MCprocess == 'QCD_PT-1400-1800':
      path = '/lustre/cms/store/user/azaza/QCD_PT-1400to1800_TuneCP5_13p6TeV_pythia8/crab_QCD_PT-1400to1800_TuneCP5_13p6TeV_pythia8_Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2/230405_184339'
   if args.dataset == 'MC' and args.MCprocess == 'QCD_PT-1800-2400':
      path = '/lustre/cms/store/user/azaza/QCD_PT-1800to2400_TuneCP5_13p6TeV_pythia8/crab_QCD_PT-1800to2400_TuneCP5_13p6TeV_pythia8_Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2/230405_184401'
   if args.dataset == 'MC' and args.MCprocess == 'QCD_PT-2400-3200':
      path = '/lustre/cms/store/user/azaza/QCD_PT-2400to3200_TuneCP5_13p6TeV_pythia8/crab_QCD_PT-2400to3200_TuneCP5_13p6TeV_pythia8_Run3Summer22MiniAODv3-124X_mcRun3_2022_realistic_v12-v2/230405_184423'

   # ZToJets
   if args.dataset == 'MC' and args.MCprocess == 'Zqq_HT-200-400':
      path = '/lustre/cms/store/user/dtroiano/ZJetsToQQ_HT200to400_TuneCP5_13TeV-madgraphMLM-pythia8/Zqq_HT200to400_ntuple/230405_132349'
   if args.dataset == 'MC' and args.MCprocess == 'Zqq_HT-400-600':
      path = '/lustre/cms/store/user/dtroiano/ZJetsToQQ_HT400to600_TuneCP5_13TeV-madgraphMLM-pythia8/Zqq_HT400to600_ntuple/230405_133832'   
   if args.dataset == 'MC' and args.MCprocess == 'Zqq_HT-600-800':
      path = '/lustre/cms/store/user/dtroiano/ZJetsToQQ_HT600to800_TuneCP5_13TeV-madgraphMLM-pythia8/Zqq_HT600to800_ntuple/230405_134357'
   if args.dataset == 'MC' and args.MCprocess == 'Zqq_HT-800-inf':
      path = '/lustre/cms/store/user/dtroiano/ZJetsToQQ_HT800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/Zqq_HT800toInf_ntuple/230405_121440'



if args.EE == 'postEE':
   if args.dataset == 'MC' and args.MCprocess == 'VBFHCC':
      path =  '/lustre/cms/store/user/azaza/VBFHToCC_M-125_TuneCP5_13p6TeV-powheg-pythia8_Run3/VBFHToCC_130X_muRun3_2023_realistic_v8_ntuple/230613_142000'

   if args.dataset == 'MC' and args.MCprocess == 'QCD_PT-120-170':
      path = '/lustre/cms/store/user/azaza/QCD_PT-120to170_TuneCP5_13p6TeV_pythia8/crab_QCD_PT-120to170_TuneCP5_13p6TeV_pythia8_Run3Summer22EEMiniAODv3/230405_183617'
   if args.dataset == 'MC' and args.MCprocess == 'QCD_PT-170-300':
      path = '/lustre/cms/store/user/azaza/QCD_PT-170to300_TuneCP5_13p6TeV_pythia8/crab_QCD_PT-170to300_TuneCP5_13p6TeV_pythia8_Run3Summer22EEMiniAODv3/230405_183640'
   if args.dataset == 'MC' and args.MCprocess == 'QCD_PT-300-470':
      path = '/lustre/cms/store/user/azaza/QCD_PT-300to470_TuneCP5_13p6TeV_pythia8/crab_QCD_PT-300to470_TuneCP5_13p6TeV_pythia8_Run3Summer22EEMiniAODv3/230405_183703'
   if args.dataset == 'MC' and args.MCprocess == 'QCD_PT-470-600':
      path = '/lustre/cms/store/user/azaza/QCD_PT-470to600_TuneCP5_13p6TeV_pythia8/crab_QCD_PT-470to600_TuneCP5_13p6TeV_pythia8_Run3Summer22EEMiniAODv3/230405_183725'
   if args.dataset == 'MC' and args.MCprocess == 'QCD_PT-600-800':
      path = '/lustre/cms/store/user/azaza/QCD_PT-600to800_TuneCP5_13p6TeV_pythia8/crab_QCD_PT-600to800_TuneCP5_13p6TeV_pythia8_Run3Summer22EEMiniAODv3/230405_183748'
   if args.dataset == 'MC' and args.MCprocess == 'QCD_PT-800-1000':
      path = '/lustre/cms/store/user/azaza/QCD_PT-800to1000_TuneCP5_13p6TeV_pythia8/crab_QCD_PT-800to1000_TuneCP5_13p6TeV_pythia8_Run3Summer22EEMiniAODv3/230405_183810'
   if args.dataset == 'MC' and args.MCprocess == 'QCD_PT-1000-1400':
      path = '/lustre/cms/store/user/azaza/QCD_PT-1000to1400_TuneCP5_13p6TeV_pythia8/crab_QCD_PT-1000to1400_TuneCP5_13p6TeV_pythia8_Run3Summer22EEMiniAODv3/230405_183833'
   if args.dataset == 'MC' and args.MCprocess == 'QCD_PT-1400-1800':
      path = '/lustre/cms/store/user/azaza/QCD_PT-1400to1800_TuneCP5_13p6TeV_pythia8/crab_QCD_PT-1400to1800_TuneCP5_13p6TeV_pythia8_Run3Summer22EEMiniAODv3/230405_183855'
   if args.dataset == 'MC' and args.MCprocess == 'QCD_PT-1800-2400':
      path = '/lustre/cms/store/user/azaza/QCD_PT-1800to2400_TuneCP5_13p6TeV_pythia8/crab_QCD_PT-1800to2400_TuneCP5_13p6TeV_pythia8_Run3Summer22EEMiniAODv3/230405_183917'
   if args.dataset == 'MC' and args.MCprocess == 'QCD_PT-2400-3200':'''

#generating the list of all .root files in given directory and subdirectories
fileList = []
for r, d, f in os.walk(path): # r=root, d=directories, f = files
    for file in f:
        if '.root' in file:
            fileList.append(os.path.join(r, file))
            print file

#prepare final script
#final_script = open("submit_analysis_"+startTime+".sh", "w")
final_script = open("submit_analysis_"+output_name+".sh", "w")
final_script.write("#!/bin/bash\n")
final_script.write("chmod 777 -R *\n")
final_script.write("cd "+output_name+"\n")

#loop to generate one .cpp+executable+batch system conf file for each group of "n" files
n_chunk = len(fileList)//args.n
print('Number of files is {0:2d}'.format(len(fileList)))
print('Number of jobs is {0:2d}'.format(n_chunk+1))
for file_index in range(n_chunk+1):
      chunk = '' 
      for idx, l in enumerate(fileList):
         if idx < args.n*(file_index+1) and idx >= args.n*file_index:
             l = l.rstrip()
             l = '        chain->AddFile("{}");\n'.format(l)
             chunk = chunk + l

      #analysis.cpp template
      with open("templates/Analysis_template.cpp", "r") as in_file:
          buf = in_file.readlines()

      cpp_filename = "Analysis_"+args.dataset+"_"+args.run+args.MCprocess+"_"+args.anatype+"_chunk"+str(file_index)+".cpp"
      with open(cpp_filename, "w") as out_file:
          for lb in buf:
              if lb == '        //AddFile_'+args.dataset+'\n':
                  #write group of files
                  out_file.write(chunk)
              elif lb == '        //OutFile_'+args.dataset+'\n':
                  #write output file name
                  out_file.write('        fileout = "'+out_filename+str(file_index)+'.root";\n')
              else: out_file.write(lb)

              #elif lb == '            TString fileout = "AddOutput_'+args.dataset+args.MCprocess+'_'+args.anatype+'.root";\n':
                  #write output file name
               #   out_file.write('        TString fileout = "'+out_filename+str(file_index)+'.root";\n')
              #else: out_file.write(lb)

      #executable template
      with open("templates/launch_analysis_template.job", "r") as launch_infile:
          buf2 = launch_infile.readlines()

      launch_filename = "launch_analysis_"+args.dataset+"_"+args.run+args.MCprocess+"_"+args.anatype+"_"+str(file_index)+".job"
      with open(output_name+"/"+launch_filename, "w") as launch_outfile:
          for lb2 in buf2:
              if lb2 == "#compile\n":
                  launch_outfile.write("cd "+output_name+"\n")
                  launch_outfile.write("g++ -I $ROOTSYS/include ../"+cpp_filename+" `root-config --glibs` `root-config --libs` `root-config --cflags` -lTMVA -L $ROOTSYS/lib -o executable"+str(file_index)+"\n")
              elif lb2 == "#execute\n":
                  launch_outfile.write('./executable'+str(file_index)+option_string+'\n')
              else: launch_outfile.write(lb2)

      #myCondor template
      with open("templates/my_HTCondor_template.job", "r") as myCondor_infile:
          buf3 = myCondor_infile.readlines()

      condor_filename = "my_HTCondor_"+args.dataset+"_"+args.run+args.MCprocess+"_"+args.anatype+"_"+str(file_index)+".job"
      with open(output_name+"/"+condor_filename, "w") as myCondor_outfile:
          for lb3 in buf3:
              if lb3 == "Executable = launch_analysis_template.job\n":
                  myCondor_outfile.write("Executable = "+launch_filename+"\n")
              else: myCondor_outfile.write(lb3)

      #add lines to final script
      final_script.write("echo condor_submit "+condor_filename+" -name ettore\n")
      final_script.write("condor_submit "+condor_filename+" -name ettore\n")

final_script.close()
#submitName = "submit_analysis_"+startTime+".sh"
#source submitName
