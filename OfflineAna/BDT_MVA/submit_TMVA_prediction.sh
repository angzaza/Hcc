#!/bin/bash

#export CVSROOT=:pserver:anonymous@cmscvs.cern.ch:/cvs_server/repositories/CMSSW
#export CVSROOT=":ext:azaza@lxplus.cern.ch:/afs/cern.ch/work/c/cmsbuild/public/cvs/CMSSW"
#export CVS_RSH=ssh
#source /cvmfs/cms.cern.ch/cmsset_default.sh

#source /afs/cern.ch/cms/LCG/LCG-2/UI/cms_ui_env_3_2.sh
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/exp_soft/crab/pbs_python-3.5.0/lib

#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib
#export STORE=/lustre/cms/store/user/azaza/

#cd /lustrehome/azaza/HccAnalysis/CMSSW_12_4_3/src/BDT_MVA
#eval `scramv1 runtime -sh`
#set env. variables if needed
#export HOME=/cvmfs/cms.cern.ch/slc7_amd64_gcc630/lcg/root/6.10.08/etc/
#cd /lustrehome/azaza/HccAnalysis/CMSSW_12_4_3/src/BDT_MVA

# Define paths and variables
#EXECUTABLE="./your_executable_name"
OUTPUT_DIR="/lustrehome/azaza/HccAnalysis/CMSSW_12_4_3/src/BDT_MVA/Output_files"
SAMPLE_DIR="/lustrehome/azaza/HccAnalysis/CMSSW_12_4_3/src/Analysis/mergedFiles/Hcc_v3Mar24" 

cd $OUTPUT_DIR

# Find all sample files in the directory
SAMPLES=($(find "$SAMPLE_DIR" -name "*.root"))

# Loop over each sample
#for sample in "${SAMPLES[@]}"; do
for sample_path in "${SAMPLES[@]}"; do
    # Extract sample name without path and extension
    #sample_name=$(basename "$sample_path" | sed 's/.root//')
    sample_name=$(basename "$sample_path" | sed 's/.root//; s/AnalysedTree_//g; s/_merged//g')
    # Construct output file path based on sample name
    OUTPUT_FILE="${OUTPUT_DIR}/rootFiles/output_${sample_name}.root"
    # create config file
    mkdir BDT_${sample_name}
    cd BDT_${sample_name}
    echo -e "creo il file launch_${sample_name}.job"
    cat > "launch_${sample_name}.job" << EOF
#!/bin/bash
export CVSROOT=:pserver:anonymous@cmscvs.cern.ch:/cvs_server/repositories/CMSSW
export CVSROOT=":ext:azaza@lxplus.cern.ch:/afs/cern.ch/work/c/cmsbuild/public/cvs/CMSSW"
export CVS_RSH=ssh
source /cvmfs/cms.cern.ch/cmsset_default.sh

export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/opt/exp_soft/crab/pbs_python-3.5.0/lib
export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/usr/lib
export STORE=/lustre/cms/store/user/azaza/


cd /lustrehome/azaza/HccAnalysis/CMSSW_12_4_3/src/BDT_MVA/templates
eval \`scramv1 runtime -sh\`
#set env. variables if needed
export HOME=/cvmfs/cms.cern.ch/slc7_amd64_gcc630/lcg/root/6.10.08/etc/
cd /lustrehome/azaza/HccAnalysis/CMSSW_12_4_3/src/BDT_MVA
pwd
echo -e "creando l'eseguibile"

cd Output_files/BDT_${sample_name}
pwd
g++ -I \$ROOTSYS/include ../../templates/TMVA_prediction.cpp \`root-config --glibs\` \`root-config --libs\` \`root-config --cflags\` -lTMVA -L \$ROOTSYS/lib -o executable0
sleep 1

./executable0 "$sample_path" "$OUTPUT_FILE"
EOF
    chmod +x launch_${sample_name}.job

    # Create submission file
    cat > "submit_${sample_name}.job" << EOF
universe = vanilla
Executable = launch_${sample_name}.job
Requirements = OpSys == "LINUX" && (Arch != "DUMMY")
Output = sleep_\$(Cluster)_\$(Process).stdout
Error = sleep_\$(Cluster)_\$(Process).stderr
Log = sleep_\$(Cluster)_\$(Process).log
request_memory = 500
notify_user = angela.zaza@ba.infn.it
Queue

EOF
 
    # Submit the job to HTCondor
    condor_submit "submit_${sample_name}.job"
    cd ..
done

cd ..
