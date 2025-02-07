Hcc, Zqq Analyzer for CMS Run3

CMSSW: 12.4.3 

------

# Install

follow install.sh instruction

# access to cms vo
    voms-proxy-init  --voms cms -valid 192:00

# Run the ntuplizer
- To run on Data 
        
        cmsRun python/templateData_Run3_Hcc_cfg.py

- To run on MC
        
        cmsRun python/templateMC_Zcc_Run3.py

!IMPORTANT: change the global tag accordingly to the dataset you use in input 
(https://github.com/BariGEMJetTau/Hcc/blob/main/HccAna/python/templateMC_Zcc_Run3.py#L17)

- To submit on CRAB use the script https://github.com/BariGEMJetTau/Hcc/blob/main/HccAna/python/crab/crab_Zcc.py

        crab submit -c crab_Zcc.py


- To submit on CRAB on multiple dataset, you need to define a txt file with the list of dataset you want to ntuplize, e.g. https://github.com/angzaza/Hcc/blob/main/Utilities/crab/dataset_MC_Summer23.txt

	then you can use a python script as https://github.com/angzaza/Hcc/blob/main/Utilities/crab/SubmitCrabJobs_MCSummer23.py
	this takes as options: -t (tag), -d (datasets), -c (crab config file)
	the crab config file for MC Summer23 is HccAna/python/crab/crabConfig_TEMPLATE_MCSummer23.py
	the crab config file takes as input the python config file for the ntuplizer https://github.com/angzaza/Hcc/blob/main/HccAna/python/templateMC_QCDZqq_Summer23_JECdB.py


Note: if you run on a private dataset, you should change config.Data.inputDBS = 'global' in config.Data.inputDBS = 'phys03' in the crab config file
