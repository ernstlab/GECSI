This is Gene Expression-based Chromatin State Imputation (GECSI) program. This program is used for generating chromatin state annotations from biosamples that have gene expression data. Here are the basic instructions for installing and using this program.

## Quick Guide ##

### Step 0: Clone This Repository
To download and use this project locally, run the following command in your terminal:

`git clone https://github.com/jingyuanf/GECSI.git`

This will create a folder named GECSI in your current working directory. You can then navigate into it:

`cd GECSI`



### Step 1: Installation
GECSI is a tool based on R version 4.1.0 with the requirement of installing multiple packages. In order to avoid potential environment conflicts, we recommend that you install the tool in a clean environment by running:

`chmod +x install_GECSI.sh`
`install_GECSI.sh`

Since it installs the R environment and all related packages from scratch, it takes ~30min to finish.

This will automatically create an environment called "gecsi-r" for you. 

After the environment is installed, whenever you want to run GECSI, please activate the environment using:

`conda activate gecsi-r`

### Step 2: Check that GECSI is correctly installed ###

After the environment is set up, to correctly execute the program, you will need to first use the following commands to make certain files executable:

`chmod +x GECSI.sh`

`chmod -R +x ./scripts/*.R`

`chmod -R +x ./preprocessing/*.sh`

In order to make sure that GECSI is correctly installed, we recommend that you run the following script that goes through the whole process using the example data provided.

`chmod +x run_GECSI_example.sh`

`./run_GECSI_example.sh`

This example script will run for about 30 minutes in an intel 32G CPU core. It will train toy models in 50 epigenomes and apply them in 10 epigenomes. Output of this example script will be stored in `example_results/` and log output will be stored in `log/`

## How to use the program? ##

### 1. Use previously trained models to apply to your gene expression data


After making sure that GECSI is correctly set up, you can use the program to train and apply your model now! Check "-h" or "--help" for help with the options. 


`chmod +x GECSI.sh`

`./GECSI.sh -h`

`./GECSI.sh --help`

Alternatively, you can check a more detailed manual here [provide a link to the manual].


Standard steps of performing GECSI is:
`./GECSI.sh -c compute_dist`

`./GECSI.sh -c find_nearest_ct`

`./GECSI.sh -c train`

`./GECSI.sh -c apply`



