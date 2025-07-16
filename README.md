This is Gene Expression-based Chromatin State Imputation (GECSI) model. Please follow the below instructions for installing and using this model.

## Standard Usage ##

### Step 1: Installation
GECSI is a tool based on R version 4.1.0 with the requirement of installing multiple packages. In order to avoid potential environment conflicts, we recommend that you install the tool in a clean environment by running:
`chmod +x install_GECSI.sh`
`install_GECSI.sh`

Since it installs the R environment and all related packages from scratch, it takes ~25min to finish.

This will automatically create an environment called "gecsi-r" for you. 

After the environment is installed, whenever you want to run GECSI, please activate the environment using:
`conda activate gecsi-r`

### Step 2: Run GECSI

In order to make sure that GECSI is correctly installed, we recommend that you 

After the environment is set up, run GECSI.sh. Check "-h" or "--help" for help with the options.
`chmod +x GECSI.sh`
`./GECSI.sh -h`
`./GECSI.sh --help`

Standard steps of performing GECSI is:
`./GECSI.sh -c compute_dist`
`./GECSI.sh -c find_nearest_ct`
`./GECSI.sh -c train`
`./GECSI.sh -c apply`

Or, if you want to use our pre-trained models, you can skip the "train" step and follow the below steps (TODO)

## Use Snakemake ##

Snakemake 

## (Optional Step): Parameter Tuning

