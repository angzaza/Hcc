Code to train and apply a BDTG for signal vs background discrimination based on the TMVA package provided by ROOT

Execute the BDT classification on root

	.x TMVA_Classification.C


Apply the model previously trained on all the samples

	source submit_TMVA_prediction.sh


Get a single output to be used as input for combine on ROOT

	.x singleOutFile_maker.C

