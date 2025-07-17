# Gene Expression-based Chromatin State Imputation (GECSI) #

This is Gene Expression-based Chromatin State Imputation (GECSI) program. This program is used for generating chromatin state annotations from biosamples that have gene expression data. Here are the basic instructions for installing and using this program.

## Quick Guide ##

### Step 0: Clone This Repository
To download and use this project locally, run the following command in your terminal:

`git clone https://github.com/jingyuanf/GECSI.git`

This will create a folder named GECSI in your current working directory. You can then navigate into it:

`cd GECSI`



### Step 1: Installation
GECSI is a tool based on R version 4.1.0 with the requirement of installing multiple packages. In order to avoid potential environment conflicts, you will need to install the tool in a clean environment by running:

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

Once this is finished, it means that GECSI is correctly set up!

## Guide for using previously trained models to apply to your gene expression data ##

Please follow the below steps sequentially and closely.

### Step 1: Download pre-trained models and training gene expression data ###

You need to download all models we trained using the International Human Epigenome Consortium (IHEC) into your working directory. 

**Important**: Before downloading, make sure that you set up a **project folder** in the **working directory** with a name you like. This will be the same working directory and project folder where GECSI generates outputs. In this guide, we use `your_wd` to specify the working directory and `your_proj_name` to specify the project folder.

`mkdir -p your_wd/your_proj_name/`

Navigate to the project folder

`cd your_wd/your_proj_name/`

Now obtain the pre-trained models using

`wget "placeholder_for_zenodo/Train_0.zip" -O Train_0.zip` or `curl -L "placeholder_for_zenodo/Train_0.zip" -o Train_0.zip`

`unzip Train_0.zip`

Similarly, you can download the training sample gene expression data using

`wget placeholder_for_zenodo/training_sample_gene_expression.tsv -O training_sample_gene_expression.tsv`

and the names of training samples using

`wget placeholder_for_zenodo/training_sample_names.txt -O training_sample_names.txt`

### Step 2: Prepare your gene expression data where GECSI needs to be applied ###
You will need to prepare your gene expression data in the following format: (you may check `data/gene_exp_test.tsv` for the correct format):

A tab-delimited file where each row is a gene and each column is a sample. The values need to be log2 transformed gene expression TPM data. The first row must be the sample names, and the first column must be **gene symbols**.

**Note**: pre-trained GECSI models currently only support human gene expression data and the identifier needs to be gene symbols (instead of Ensembl IDs) to make it compatible with the pre-trained data. The predicted annotations will be in Genome Assembly hg38.

Example:

```
Sample1	Sample2	Sample3
TP53	5.2	    6.1	    4.8
BRCA1	2.3	    2.9	    3.0
GAPDH	10.1	11.0	9.8
```

### Step 3: Compute distance between new samples and training samples ###

Use the following command to compute the distance between new samples and training samples. **Note**: Make sure to stick with the same working directory and project name as in Step 1 (so that the output will be generated in the same folder as where you downloaded the models). 

`./GECSI.sh -c compute_dist -d --gexp-ref "/path/to/training_sample_gene_expression.tsv" --gexp-apply "/path/to/your_sample_gene_expression.tsv" --proj-name "your_proj_name" -o "your_wd"`

### Step 4: Find nearest samples between new samples and training samples ###

Use the following command to find the nearest training samples for the new samples. **Note**: Make sure that "your_sample_names.txt" contain the exact same sample names as in the gene expression data file.

`./GECSI.sh -c find_nearest --ref-list "/path/to/training_sample_names.txt" --apply-list "/path/to/your_sample_names.txt" --proj-name "your_proj_name" -o "your_wd" --nref 414 `

### Step 5: Process training sample chromatin states and apply pre-trained models to your new data ###

**Caution**: Preprocessing steps will take up to **30G** of storage. Make sure you preserve enough space for running these commands.

1. You will need to first download all training sample chromatin state annotations to generate features for prediction. 

    `wget placeholder_for_chromhmm/training_sample_chromatin_states.zip -O training_sample_chromatin_states.zip`

2. Unzip the zip file. You will obtain a folder named `training_sample_chromatin_states/` that contain all training sample chromatin state bed files.

    `unzip training_sample_chromatin_states.zip`

3. Bin the chromatin state files using the following command into 200 basepair into a separate folder called `training_sample_chromatin_states_binned/`:

    `./preprocessing/bin_bed.sh "/path/to/training_sample_chromatin_states/" 200 "/path/to/training_sample_chromatin_states_binned/" hg38`

4. Now you can apply the pre-trained GECSI models on your new samples. 

    **Important**: The following command performs application for each new sample and for each chromosome (specified in `--apply-sam` and `--apply-chr` command). You may want to parallelize this process if you want to generate outputs for multiple samples and chromosomes.


    `./GECSI.sh -c apply --chrstate "/path/to/training_sample_chromatin_states_binned/" --apply-sam "sample_name" --ref-list "/path/to/training_sample_names.txt" --proj-name "your_proj_name" -o "your_wd" -k "5" --ref-chr "all" --apply-chr "chr" --sample-size "100000" --nref "5" --lambda "0.0001"`

Once it's finished, you will see the predicted chromatin state annotations in `your_wd/your_proj_name/Apply_0/predictions/chr/`. The results will be stored in ".rds" format in this folder and ".bed.gz" format in the subfolder `bed_files/`.

****

## Guide for training your own model using your gene expression and chromatin state data ##
#TODO
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



