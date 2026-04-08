# NBES
NBES
How to run the model:
1. go to "scripts"-folder
2. open "NBESOperateModel.R"
3. install missing packages if necessary
4. set parameter/simulation settings in the script (currently set to default values)
5. run the model at the end of the script

The model output will be stored in the output folder. The script "NBESProcessOutput.R" can be used to process and summarise the simulation results for further analysis. 
The file names for reading/writing have to be adjusted accordingly.


# Data
Simulated data used for the analysis are available on Zenodo https://doi.org/10.5281/zenodo.15274625 and need to be stored in output if you wish to replicate the analysis. 

# Analysis
Scripts for analysis is stored in the analysis folder, please follow the order of the rscripts starting with the slope calulcations