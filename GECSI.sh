#!/bin/bash

# Function to print usage
help() {
    echo "Usage: $0 --command {compute_dist|find_nearest|train|apply|help} [options...]"
    echo "Options:"
    echo "  -c | --command <command>       Command to execute (compute_dist, find_nearest, train, apply, help)"
    echo "General Commands: "
    echo "  -v | --verbose            Enable verbose output"
    echo "  -d | --debug              Enable debug mode"
    echo "  --chr <chr>               Chromosome to use (e.g., 'chr1', 'chrX', 'all') (defalt: 'all')"
    echo "  --proj-name <name>        Project name. Distance metrics will be stored in a folder named by this project name, inside a larger folder defined by output directory. (default: 'new_proj')"
    echo "  -o | --outdir <dir>       Output directory (Required)"
    echo "  --seed <number>           Seed for randomization (default:42)"
    echo "  --task-id <string>        A unique identifier for a sequence of steps. Separates different experiments in the same project folder.(default:'1')"
    echo "  --overwrite               Whether to overwrite existing models and predictions (default:false)"
    echo "  --parallel                Number of threads to parallelly run some steps (default:3)"
    echo "  --num-states <number>     Number of chromatin states to use in the model (default:18)"
    echo "  --quies-state <name>      Name of the quiescent state (default: '18_Quies')"
    echo "  --states-list <file>      A file providing a list of available chromatin state names in order, separated in lines. If not provided, will default to the Roadmap 18-state chromatin state ordering."

    echo "Commands for compute_dist: "
    echo "  --dist <method>           Distance Method to calculate nearest samples. 'pcc' for Pearson Correlation, 'scc' for Spearman Correlation, and 'euc' for Euclidean distance) (default: 'pcc')"
    echo "  --gexp-ref <file>         File path for Gene Expression data used for reference. In tab-delimited format. (Required)"
    echo "  --gexp-apply <file>       File path for Gene Expression data used for application. In tab-delimited format. (Required)"
    echo "  -v, -d, --proj-name, -o, --task-id as described in 'General Commands'"

    echo "Commands for find_nearest: "
    echo "  --ref-list <file>         File listing sample names that can be used for reference. "
    echo "  --nref <nref>             Number of nearest reference samples to find and train the model. (default:100)"
    echo "  --apply-list <file/name>   File (txt) containing a list of sample names for which GECSI will be applied. Or a single sample name. (Required)"
    echo "  -v, -d, --proj-name, -o, --task-id as described in 'General Commands'"

    echo "Commands for train: "
    echo "  --ref-list <file>         File listing sample names that can be used for reference. "
    echo "  --train-sam <file/name>   A sample name (recommended) or a file (txt) containing a list of sample names for which GECSI will be trained on. (Required)"
    echo "  --chrstate <dir>          A directory containing reference chromatin state annotation data. File names should contain sample names for reference purpose. (Required)"
    echo "  --sample_size <integer>   Sample size (number of positions) used for training. (default: 500000)"
    echo "  --dist <method>           Distance Method to calculate nearest samples. 'pcc' for Pearson Correlation, 'scc' for Spearman Correlation, and 'euc' for Euclidean distance) (default: 'pcc')"
    echo "  -k | --k-nearest <integer> Number of nearest reference celltypes to use as features (default: 10)"
    echo "  --save-train              Whether to save predictions for training cell types. (default: false)"
    echo "  -v, -d, --chr, --overwrite, --parallel, --ref-list,--proj-name, -o, --task-id, --seed, --num_states, --quies_state, --states_list as described in 'General Commands'"
    
    echo "Commands for apply: "
    echo "  --ref-chr <chr>           Chromosome used in training step (e.g., 'chr1', 'chrX', 'all') (defalt: 'all')"
    echo "  --apply-chr <chr>         Chromosome to use in application step (e.g., 'chr1', 'chrX', 'all') (defalt: 'chr1') (Note: "all" is not recommended. Parallelization is preferred.)"
    echo "  --nref <integer>          Number of nearest reference samples to find and train the model. (default:100)"
    echo "  --chrstate <dir>          A directory containing reference chromatin state annotation data. File names should contain sample names for reference purpose. (Required)"
    echo "  --apply-sam <nref>        A sample name (recommended) or a file (txt) containing a list of sample names for which GECSI will be applied to. (Required)"
    echo "  --sample-size <integer>   Sample size (number of positions) that was used for training. (default:500000)"
    echo "  --ignore-missing-model    When encountering models that are missing, ignore and use the other models instead of reporting error"
    echo "  --lambda <double>          Lasso regularization parameter (lambda)"
    echo "  -k | --k-nearest <integer> Number of nearest reference celltypes used as features (default: 10)"
    echo "  -v, -d, --overwrite, --parallel, --ref-list, --proj-name, -o, --task-id, --seed as described in 'General Commands'"

    exit 0
}

usage() {
    help
    exit 1
}

# Initialize variables
command=""
verbose=false
debug=false
chr="all"
ref_list="default"
proj_name="new_proj"
seed=42
task_id="1"
dist="scc"
nref=100
gexp_ref=""
gexp_apply=""
outdir=""
apply_list=""
train_sam=""
chrstate=""
sample_size=500000
k=10
save_train=false
overwrite=false
ref_chr="all"
apply_chr="chr1"
apply_sam=""
ignore_missing_model=false
parallel=3
lambda=0.001
nstates=18
quies_state="18_Quies"
states_list=""



mydir="$(dirname "$(realpath "$0")")"
cd $mydir

# Parse options

# Parse options
while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
        -c|--command)
        command="$2"
        shift # past argument
        shift # past value
        ;;
        --chr)
        chr="$2"
        shift # past argument
        shift # past value
        ;;
        --ref-list)
        ref_list="$2"
        shift # past argument
        shift # past value
        ;;
        --gexp-ref)
        gexp_ref="$2"
        shift # past argument
        shift # past value
        ;;
        --gexp-apply)
        gexp_apply="$2"
        shift # past argument
        shift # past value
        ;;
        --proj-name)
        proj_name="$2"
        shift # past argument
        shift # past value
        ;;
        -o|--outdir)
        outdir="$2"
        shift # past argument
        shift # past value
        ;;
        --seed)
        seed="$2"
        shift # past argument
        shift # past value
        ;;
        --task-id)
        task_id="$2"
        shift # past argument
        shift # past value
        ;;
        --dist)
        dist="$2"
        shift # past argument
        shift # past value
        ;;
        --nref)
        nref="$2"
        shift
        shift
        ;;
        --apply-list)
        apply_list="$2"
        shift
        shift
        ;;
        --train-sam)
        train_sam="$2"
        shift
        shift
        ;;
        -k|--k_nearest)
        k="$2"
        shift
        shift
        ;;
        --lambda)
        lambda="$2"
        shift
        shift
        ;;
        --ref-chr)
        ref_chr="$2"
        shift
        shift
        ;;
        --apply-chr)
        apply_chr="$2"
        shift
        shift
        ;;
        --apply-sam)
        apply_sam="$2"
        shift
        shift
        ;;
        --sample-size)
        sample_size="$2"
        shift
        shift
        ;;
        --num-states)
        nstates="$2"
        shift
        shift
        ;;
        --quies-state)
        quies_state="$2"
        shift
        shift
        ;;
        --states-list)
        states_list="$2"
        shift
        shift
        ;;
        --ignore-missing-model)
        ignore_missing_model=true
        shift
        ;;
        --chrstate)
        chrstate="$2"
        shift
        shift
        ;;
        --parallel)
        parallel="$2"
        shift
        shift
        ;;        
        -v|--verbose)
        verbose=true
        shift # past argument
        ;;
        -d|--debug)
        debug=true
        shift # past argument
        ;;
        --save_train)
        save_train=true
        shift
        ;;
        --overwrite)
        overwrite=true
        shift
        ;;
        -h|--help)
        help
        ;;
        *)
        echo "Unknown option: $key"
        help
        ;;
    esac
done

# Check that a command was provided
if [ -z "$command" ]; then
    echo "Error: A command must be specified with -c."
    usage
fi

# Enable verbose output if the flag is set
if [ "$verbose" = true ]; then
    echo "Verbose mode enabled."
fi


if [ "$verbose" = true ]; then
    echo "command: $command"
    echo "ref_list: $ref_list"
    echo "gexp_ref: $gexp_ref"
    echo "gexp_apply: $gexp_apply"
    echo "proj_name: $proj_name"
    echo "outdir: $outdir"
    echo "seed: $seed"
    echo "task_id: $task_id"
    echo "dist: $dist"
    echo "verbose: $verbose"
    echo "debug: $debug"
    echo "nref: $nref"
    echo "apply_list: $apply_list"
    echo "train_sam: $train_sam"
    echo "chrstate: $chrstate"
    echo "sample_size: $sample_size"
    echo "k: $k"
    echo "lambda: $lambda"
    echo "save_train: $save_train"
    echo "overwrite: $overwrite"
    echo "ignore_missing_model: $ignore_missing_model"
    echo "apply_sam: $apply_sam"
    echo "ref_chr: $ref_chr"
    echo "parallel: $parallel"
fi

# Enable debug mode if the flag is set
if [ "$debug" = true ]; then
    echo "Debug mode enabled."
    set -x  # Enable debugging output
fi


# Define functions
compute_dist() {
    if [ "$verbose" = true ]; then
        echo "gene_exp_ref: $gexp_ref"
        echo "gene_exp_apply: $gexp_apply"
        echo "outdir: $outdir"
    fi

    if [ -z "$gexp_ref" ]; then
        echo "Error: The gene expression for training command requires a string of file path (-gt)."
        usage
    fi
    if [ -z "$gexp_apply" ]; then
        echo "Error: The gene expression for applying command requires a string of file path (-ga)."
        usage
    fi
    if [ -z "$outdir" ]; then
        echo "Error: The output directory command requires a string of file path (-o)."
        usage
    fi    
    compute_dist_code_name="./scripts/compute_dist.R"
    if [ "$verbose" = true ]; then
        ${compute_dist_code_name} -v \
                            --wd $mydir \
                            --gene_train ${gexp_ref} \
                            --gene_test ${gexp_apply} \
                            --proj_name ${proj_name} \
                            -r ${outdir} \
                            --dist ${dist} \
                            --task_id ${task_id} \

    else
        ${compute_dist_code_name} -q \
                            --wd $mydir \
                            --gene_train ${gexp_ref} \
                            --gene_test ${gexp_apply} \
                            --proj_name ${proj_name} \
                            -r ${outdir} \
                            --dist ${dist} \
                            --task_id ${task_id} \

    fi
    
}

find_nearest_ct() {
    if [ "$verbose" = true ]; then
        echo "ref_list: $ref_list"
        echo "apply_list: $apply_list"
        echo "outdir: $outdir"
    fi

    if [ -z "$ref_list" ]; then
        echo "Error: The samples for training command requires a string of file path or sample name (--ref-list)."
        usage
    fi
    if [ -z "$apply_list" ]; then
        echo "Error: The samples for application command requires a string of file path or sample name (--apply-list)."
        usage
    fi
    if [ -z "$outdir" ]; then
        echo "Error: The output directory command requires a string of file path (-o)."
        usage
    fi   

    generate_ct_name="./scripts/find_nearest_ct.R"
    if [ "$verbose" = true ]; then
        ${generate_ct_name} -v \
                            --training_ct ${ref_list} \
                            --wd $mydir \
                            --holdout ${apply_list} \
                            --proj_name ${proj_name} \
                            -r ${outdir} \
                            --num_ct_to_train ${nref} \
                            --task_id ${task_id} \
                            --dist ${dist} \

    else
        ${generate_ct_name} -q \
                            --training_ct ${ref_list} \
                            --wd $mydir \
                            --holdout ${apply_list} \
                            --proj_name ${proj_name} \
                            -r ${outdir} \
                            --num_ct_to_train ${nref} \
                            --task_id ${task_id} \
                            --dist ${dist} \

    fi
}

train() {
    if [ "$verbose" = true ]; then
        echo "train_sam: $train_sam"
        echo "chrstate: $chrstate"
        echo "sample_size: $sample_size"
        echo "k: $k"
        echo "save_train: $save_train"
        echo "overwrite: $overwrite"
        echo "ref_list: $ref_list"
        echo "outdir: $outdir"
        echo "num_states: $nstates"
        echo "quies_state: $quies_state"
        echo "states_list: $states_list"
    fi

    if [ -z "$ref_list" ]; then
        echo "Error: The reference list command requires a string of file path (--ref-list)."
        usage
    fi
    if [ -z "$train_sam" ]; then
        echo "Error: The sample(s) for training command requires a string of file path or sample name (--train-sam)."
        usage
    fi
    if [ -z "$chrstate" ]; then
        echo "Error: The chromatin state directory command requires a string of file path (--chrstate)."
        usage
    fi    

    train_name="./scripts/train_GECSI.R"

    if [ -z "$states_list" ]; then
        /usr/bin/time -v ${train_name} -v ${verbose} \
                            --chr ${chr} \
                            --training_ct ${train_sam} \
                            --wd $mydir \
                            --all_training_samples ${ref_list} \
                            --proj_name ${proj_name} \
                            -r ${outdir} \
                            -c ${chrstate} \
                            --sample_size ${sample_size} \
                            --seed ${seed} \
                            --task_id ${task_id} \
                            -k ${k} \
                            --save_train ${save_train} \
                            --dist ${dist} \
                            --overwrite ${overwrite} \
                            --num_states ${nstates} \
                            --quies_state ${quies_state} \
                            --parallel ${parallel} \

    else
        /usr/bin/time -v ${train_name} -v ${verbose} \
                            --chr ${chr} \
                            --training_ct ${train_sam} \
                            --wd $mydir \
                            --all_training_samples ${ref_list} \
                            --proj_name ${proj_name} \
                            -r ${outdir} \
                            -c ${chrstate} \
                            --sample_size ${sample_size} \
                            --seed ${seed} \
                            --task_id ${task_id} \
                            -k ${k} \
                            --save_train ${save_train} \
                            --dist ${dist} \
                            --overwrite ${overwrite} \
                            --num_states ${nstates} \
                            --quies_state ${quies_state} \
                            --states_list ${states_list} \
                            --parallel ${parallel} \

    fi
}

apply() {

    if [ "$verbose" = true ]; then
        echo "ref_chr: $ref_chr"
        echo "nref: $nref"
        echo "apply_sam: $apply_sam"
        echo "k: $k"
        echo "sample_size: $sample_size"
        echo "overwrite: $overwrite"
        echo "ref_list: $ref_list"
        echo "outdir: $outdir"
        echo "chrstate: $chrstate"
        echo "lambda: $lambda"
        echo "num_states: $nstates"
        echo "states_list: $states_list"
    fi

    if [ -z "$ref_list" ]; then
        echo "Error: The reference list command requires a string of file path (--ref-list)."
        usage
    fi
    if [ -z "$apply_sam" ]; then
        echo "Error: The sample(s) for application command requires a string of file path or sample name (--apply-sam)."
        usage
    fi
    if [ -z "$chrstate" ]; then
        echo "Error: The chromatin state directory command requires a string of file path (--chrstate)."
        usage
    fi

    apply_name="./scripts/apply_GECSI.R"

    if [ -z "$states_list" ]; then
        /usr/bin/time -v ${apply_name} -v ${verbose}\
                            --train_chr ${ref_chr} \
                            --test_chr ${apply_chr} \
                            --num_train_ct ${nref} \
                            --wd $mydir \
                            --holdout ${apply_sam} \
                            --proj_name ${proj_name} \
                            -r ${outdir} \
                            -c ${chrstate} \
                            --sample_size ${sample_size} \
                            --seed ${seed} \
                            --task_id ${task_id} \
                            -k ${k} \
                            --lambda ${lambda} \
                            --ignore_missing_model ${ignore_missing_model} \
                            --dist ${dist} \
                            --training_samples ${ref_list} \
                            --parallel ${parallel} \
                            --overwrite ${overwrite} \
                            --num_states ${nstates} \

    else
        /usr/bin/time -v ${apply_name} -v ${verbose}\
                            --train_chr ${ref_chr} \
                            --test_chr ${apply_chr} \
                            --num_train_ct ${nref} \
                            --wd $mydir \
                            --holdout ${apply_sam} \
                            --proj_name ${proj_name} \
                            -r ${outdir} \
                            -c ${chrstate} \
                            --sample_size ${sample_size} \
                            --seed ${seed} \
                            --task_id ${task_id} \
                            -k ${k} \
                            --lambda ${lambda} \
                            --ignore_missing_model ${ignore_missing_model} \
                            --dist ${dist} \
                            --training_samples ${ref_list} \
                            --parallel ${parallel} \
                            --overwrite ${overwrite} \
                            --num_states ${nstates} \
                            --states_list ${states_list} \
    
    fi

}


# Execute the appropriate function based on the command
case $command in
    compute_dist)
        compute_dist
        ;;
    find_nearest)
        find_nearest_ct
        ;;
    train)
        train
        ;;
    apply)
        apply
        ;;
    help)
        help
        ;;
    *)
        echo "Unknown command: $command"
        usage
        ;;
esac

# Disable debug mode after execution
if [ "$debug" = true ]; then
    set +x  # Disable debugging output
fi

