# Gene Expression-based Chromatin State Imputation (GECSI) #

**Gene Expression-based Chromatin State Imputation (GECSI)** is a tool for generating chromatin state annotations from biosamples with gene expression data. Here are the basic instructions for installing and using this program.

## Table of Contents

### [Installation](#how-do-i-install-it)
- [Step 0: Clone This Repository](#step-0-clone-this-repository)
- [Step 1: Installation](#step-1-installation)
- [Step 2: Check that GECSI is Correctly Installed](#step-2-check-that-gecsi-is-correctly-installed)

### [Using Pre-Trained Models](#how-do-i-use-previously-trained-models-to-apply-to-my-own-gene-expression-data)
- [Step 1: Download Pre-Trained Data](#step-1-download-pre-trained-data-and-training-gene-expression-data)
- [Step 2: Prepare Your Gene Expression Data](#step-2-prepare-your-gene-expression-data)
- [Step 3: Compute Sample Distance](#step-3-compute-sample-distance)
- [Step 4: Find Nearest Samples](#step-4-find-nearest-samples)
- [Step 5: Apply Pre-Trained Models](#step-5-apply-pre-trained-models)

### [Train Your Own Model](#how-do-i-train-my-own-model)
- [Step 1: Prepare Data](#step-1-prepare-data)
  - [Chromatin States](#training-samples-chromatin-states-data)
  - [Gene Expression](#training-samples-and-applying-samples-gene-expression-data)
- [Step 2: Compute Sample Distance](#step-2-compute-sample-distance)
- [Step 3: Find Nearest Samples](#step-3-find-nearest-samples)
- [Step 4: Train GECSI Model](#step-4-train-gecsi-model)
- [Step 5: Apply GECSI Model](#step-5-apply-gecsi-model)

### [Troubleshooting and Support](#questions-or-issues)

# How do I install it? #

## Step 0: Clone This Repository
To download and use this project locally, run the following command in your terminal:

    git clone https://github.com/jingyuanf/GECSI.git

This will create a folder named GECSI in your current working directory. You can then navigate into it:

    cd GECSI


## Step 1: Installation
GECSI is a tool based on R version 4.1.0 with the requirement of installing multiple packages. In order to avoid potential environment conflicts, you will need to install the tool in a clean environment by running:

    chmod +x ./install_GECSI.sh

    ./install_GECSI.sh

Note: Currently this tool is supported on linux environment. The environment set up requires mamba, so if your environment doesn't have mamba, this code will automatically install it for you. If mamba install fails, you may need to manually install mamba and rerun this command. This code takes ~30min to finish.

This will automatically create an environment called "gecsi-r" for you. 

After the environment is installed, whenever you want to run GECSI, please activate the environment using:

    conda activate gecsi-r

## Step 2: Check that GECSI is correctly installed ##

In order to make sure that GECSI is correctly installed, we recommend that you run the following script that goes through the whole process using the example data provided.

    ./run_GECSI_example.sh

This example script will run for about 30 minutes in an intel 32G CPU core. It will train toy models in 50 epigenomes and apply them in 10 epigenomes. Output of this example script will be stored in `example_results/` and log output will be stored in `log/`

Once this is finished, it means that GECSI is correctly set up!

# How do I use previously trained models to apply to my own gene expression data? #

The models can be downloaded from `./models/Train_0/models/lr-model-multinom/all`. However, currently the International Human Epigenome Consortium (IHEC) EpiAtlas data has not been released so using previously trained models is not yet supported as it requires using unreleased data. Once the data is released at https://ihec-epigenomes.org/epiatlas/data/, we will release all files required for applying previously trained models. 

Once all files are released, you can use previously trained models by following the steps below:

## Step 1: Download pre-trained data and training gene expression data ##

* You need to download all models we trained using the International Human Epigenome Consortium (IHEC) into your working directory. 

* **Important**: Before downloading, make sure that you set up a **project folder** in the **working directory** with a name you like. This will be the same working directory and project folder where GECSI generates outputs. In this guide, we use `your_wd` to specify the working directory and `your_proj_name` to specify the project folder.

        mkdir -p your_wd/your_proj_name/

* Navigate to the project folder

        cd your_wd/your_proj_name/

* Now obtain the pre-trained models using

        wget "placeholder_link/Pretrained.zip" -O Pretrained.zip
    
    **or** 

        curl -L "placeholder_link/Pretrained.zip" -o Pretrained.zip

* Unzip the zip file

        unzip Pretrained.zip

If the above downloading fails, you can download them manually.

## Step 2: Prepare Your Gene Expression Data ##
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

## Step 3: Compute Sample Distance ##

Use the following command to compute the distance between new samples and training samples. **Note**: Make sure to stick with the same working directory and project name as in Step 1 (so that the output will be generated in the same folder as where you downloaded the models). 

    ./GECSI.sh -c compute_dist -d --gexp-ref "your_wd/your_proj_name/training_sample_gene_expression.tsv" --gexp-apply "/path/to/your_sample_gene_expression.tsv" --proj-name "your_proj_name" -o "your_wd"

## Step 4: Find Nearest Samples ##

Use the following command to find the nearest training samples for the new samples. **Note**: Make sure that "your_sample_names.txt" contain the exact same sample names as in the gene expression data file.

    ./GECSI.sh -c find_nearest --ref-list "your_wd/your_proj_name/training_sample_names.txt" --apply-list "/path/to/your_sample_names.txt" --proj-name "your_proj_name" -o "your_wd" --nref 5

## Step 5: Apply Pre-trained Models ##

**Caution**: These steps will take up to **30G** of storage. Make sure you preserve enough space for running these commands.


Bin the chromatin state files using the following command into 200 basepair into a separate folder called `training_sample_chromatin_states_binned/`:

```
./preprocessing/bin_bed.sh "your_wd/your_proj_name/training_sample_chromatin_states/" 200 "your_wd/your_proj_name/training_sample_chromatin_states_binned/" hg38 4
```

Now you can apply the pre-trained GECSI models on your new samples. 

**Important**: The following command performs application for each new sample and for each chromosome (specified in `--apply-sam` and `--apply-chr` command). You may want to parallelize this process if you want to generate outputs for multiple samples and chromosomes using GNU Parallel or a job array.


```
./GECSI.sh -c apply --chrstate "your_wd/your_proj_name/training_sample_chromatin_states_binned/" --apply-sam "<apply-sam>" --ref-list "your_wd/your_proj_name/training_sample_names.txt" --proj-name "your_proj_name" -o "your_wd" -k "5" --ref-chr "all" --apply-chr "<chr>" --sample-size "100000" --nref "5" --lambda "0.0001"
```

Once it's finished, you will see the predicted chromatin state annotations in `your_wd/your_proj_name/Apply_0/predictions/<chr>/`. The results will be stored in ".rds" format in this folder and ".bed.gz" format in the subfolder `bed_files/`.




# How do I train my own model? #

Follow the steps below for a standard setting of the training and application procedure. Check `./GECSI.sh -h` for more configuration options. 

## Step 1: Prepare data ##

### Training Samples Chromatin States data ###

Chromatin states files must be in **.bed** or **.bed.gz** format and stored in an input directory (for example `"/path/to/training_sample_chromatin_states/"`). 

To preprocess the chromatin state data, run the following command. The binned chromatin states will be saved to `"/path/to/training_sample_chromatin_states_binned/"`

    ./preprocessing/bin_bed.sh "/path/to/training_sample_chromatin_states/" $BINNING_RESOLUTION "/path/to/training_sample_chromatin_states_binned/" $GENOME_ASSEMBLY #CHROMATIN_STATE_COLUMN

For example, to use a binning resolution of 200 bp, a genome assembly of hg38, and to extract chromatin states from the 4th column (1-based indexing):

    ./preprocessing/bin_bed.sh "/path/to/training_sample_chromatin_states/" 200 "/path/to/training_sample_chromatin_states_binned/" hg38 4

**Note**: All chromatin state files should span the same genomic regions across chromosomes (i.e., consistent chromosome coverage and coordinates).

Optional Check: Ensure all chromatin state BED files cover the same genomic regions:

    bash preprocessing/check_bed_range_consistency.sh "/path/to/training_sample_chromatin_states/"

### Training Samples and Applying Samples Gene Expression data ###

Follow the same procedure as in [Step 2: Prepare Your Gene Expression Data](#step-2-prepare-your-gene-expression-data) in the previous section.

## Step 2: Compute Sample Distance ##

`"your_wd"` is the working directory you will place your project folder in. `"your_proj_name"` is the name of project folder where all input and output data will be stored.

    ./GECSI.sh -c compute_dist -d --gexp-ref "/path/to/training_sample_gene_expression.tsv" --gexp-apply "/path/to/applying_sample_gene_expression.tsv" --proj-name "your_proj_name" -o "your_wd"

## Step 3: Find Nearest Samples ##

**Note**: Make sure that `"training_sample_names.txt"` and `"applying_sample_names.txt"` contain sample names exactly matching the headers in the gene expression files. `NUM_REF` specifies the number of nearest reference samples used for model assembly.

    ./GECSI.sh -c find_nearest --ref-list "/path/to/training_sample_names.txt" --apply-list "/path/to/applying_sample_names.txt" --proj-name "your_proj_name" -o "your_wd" --nref NUM_REF

## Step 4: Train GECSI model ##

Check `./GECSI.sh -h` for a full list of configurable options.

**Note**: The example below assumes you are using the **Roadmap Epigenomics 18-state model**.
If you are using a different model, update the following parameters accordingly:
`--states-list`, `--num-states`, and `--quies-state`.

    ./GECSI.sh -c train 
                --chr "all" \
                --chrstate "/path/to/training_sample_chromatin_states_binned/" \
                --train-sam "TRAINING-SAM" <training is done one sample at a time> \
                --ref-list  "/path/to/training_sample_names.txt" \
                --proj-name "your_proj" 
                -o "your_wd" \
                --task-id 0 \ # Specify any integer if you would like to separate for different tasks
                -k 5 \ # Number of features in the model
                --sample-size 100000 \ # Number of positions sampled for model training 
                --states-list "./data/chr_state_list.txt" \ # A file containing chromatin state names separated by lines 
                --num-states 18 \ # Total number of states
                --quies-state "18_Quies" \ # The name for quiescent state

## Step 5: Apply GECSI model ##

Check [Step 5: Apply Pre-trained Models](#step-5-apply-pre-trained-models) for more details.

Check `./GECSI.sh -h` for a full list of configurable options.

```
./GECSI.sh -c apply \
            --chrstate "/path/to/training_sample_chromatin_states_binned/" \
            --apply-sam "<apply-sam>" \
            --ref-list "your_wd/your_proj_name/training_sample_names.txt" \
            --proj-name "your_proj_name" \
            -o "your_wd" 
            -k "5" \
            --ref-chr "all" \
            --apply-chr "<chr>" \ # Specify a chromosome
            --sample-size "100000" 
            --nref "5" \
            --lambda "0.0001" \
            --num-states 18 \ # Total number of states
            --states-list "./data/chr_state_list.txt" \ # A file containing chromatin state names separated by lines 
```

Once it's finished, you will see the predicted chromatin state annotations in `your_wd/your_proj_name/Apply_<task_id>/predictions/<chr>/`. The results will be stored in ".rds" format in this folder and ".bed.gz" format in the subfolder `bed_files/`.


# Questions or Issues? #

If you encounter problems or have questions, please open an issue on the GitHub repository or contact the maintainer.


