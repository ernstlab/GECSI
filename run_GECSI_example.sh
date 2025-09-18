#!/bin/bash

# set -e  # Exit on error
# trap 'echo -e "\n[ERROR] Line $LINENO: $BASH_COMMAND failed. See logs for details." >&2' ERR


# Set the working directory to the directory of the script
SCRIPT_DIR="$(dirname "$0")"

# Change to that directory
cd "$SCRIPT_DIR" || exit

# Input paths and parameters
gene_exp_train="./data/gene_exp_train.tsv"
gene_exp_test="./data/gene_exp_test.tsv"
proj_name="GECSI_example_project"
outdir="./example_results"
logdir="./log"
dist_method="scc"
seed=42
task_id=0
ref_list="./data/training_ct.txt"
apply_list="./data/testing_ct.txt"
chrstate_dir="./data/training_ct_chromhmm_bins"
state_list="./data/chr_state_list.txt"
k=5
nref=5
train_chr="chr1"
apply_chr="chr1"

# Create output directory
mkdir -p "$outdir"
mkdir -p "$logdir"

# Count dynamic steps (train + apply)
num_train=$(wc -l < "$ref_list")
num_apply=$(wc -l < "$apply_list")

# Define total number of steps (static + dynamic)
total_steps=$((1 + 1 + 1 + num_train + num_apply))  # preprocessing + compute_dist + find_nearest + train + apply
current_step=1

# Progress bar function
# Progress bar function
print_progress() {
    bar_length=50  # total number of characters in the bar
    percent=$((100 * current_step / total_steps))
    filled=$((bar_length * current_step / total_steps))
    empty=$((bar_length - filled))

    bar=$(printf "%0.s#" $(seq 1 $filled))
    space=$(printf "%0.s " $(seq 1 $empty))

    printf "\rProgress: [%s%s] %3d%% (%d/%d)" "$bar" "$space" "$percent" "$current_step" "$total_steps"
    ((current_step++))
}

print_progress
./preprocessing.sh > ${logdir}/preprocessing.log 2>&1
status=$?
if [[ $status -ne 0 ]]; then
    echo "[ERROR] Preprocessing failed (exit code $status). Check ${logdir}/preprocessing.log" >&2
    exit $status
fi

chmod +x ./GECSI.sh
print_progress
./GECSI.sh -c compute_dist -d --gexp-ref "${gene_exp_train}" --gexp-apply "${gene_exp_test}" --proj-name "${proj_name}" --seed ${seed} --dist "${dist_method}" -o "$outdir" > ${logdir}/compute_dist.log 2>&1
status=$?
if [[ $status -ne 0 ]]; then
    echo "[ERROR] GECSI compute_dist failed (exit code $status). Check ${logdir}/compute_dist.log" >&2
    exit $status
fi

print_progress
./GECSI.sh -c find_nearest --ref-list "${ref_list}" --apply-list "${apply_list}" --proj-name "${proj_name}" --dist "${dist_method}" -o "$outdir" --task-id "$task_id" --nref "$nref" > ${logdir}/find_nearest.log 2>&1
status=$?
if [[ $status -ne 0 ]]; then
    echo "[ERROR] GECSI find_nearest failed (exit code $status). Check ${logdir}/find_nearest.log" >&2
    exit $status
fi

for training_ct in $(cat "./data/training_ct.txt"); do
    echo "Training on cell type: $training_ct" > ${logdir}/train_${training_ct}.log
    ./GECSI.sh -c train --chr "chr1" --chrstate "${chrstate_dir}" --train-sam "${training_ct}" --ref-list "${ref_list}" --proj-name "${proj_name}" --dist "${dist_method}" -o "$outdir" --task-id "$task_id" -k "$k" --sample-size "1000" --states-list $state_list --num-states 18 --quies-state "18_Quies" --overwrite >> ${logdir}/train_${training_ct}.log 2>&1
    status=$?
    if [[ $status -ne 0 ]]; then
        echo "[ERROR] GECSI train failed for ${training_ct} (exit code $status). Check ${logdir}/train_${training_ct}.log" >&2
        exit $status
    fi
    echo "Training completed for cell type: $training_ct" >> ${logdir}/train_${training_ct}.log
    print_progress
done

for testing_ct in $(cat "./data/testing_ct.txt"); do
    echo "Applying model to cell type: $testing_ct" > ${logdir}/apply_${testing_ct}.log
    ./GECSI.sh -c apply --chrstate "${chrstate_dir}" --apply-sam "${testing_ct}" --ref-list "${ref_list}" --proj-name "${proj_name}" --dist "${dist_method}" -o "$outdir" --task-id "$task_id" -k "$k" --ref-chr $train_chr --apply-chr $apply_chr --sample-size "1000" --nref "$nref" --lambda "0.0001" --num-states 18 --states-list $state_list --overwrite --probability >> ${logdir}/apply_${testing_ct}.log 2>&1
    status=$?
    if [[ $status -ne 0 ]]; then
        echo "[ERROR] GECSI apply failed for ${testing_ct} (exit code $status). Check ${logdir}/apply_${testing_ct}.log" >&2
        exit $status
    fi
    echo "Application completed for cell type: $testing_ct" >> ${logdir}/apply_${testing_ct}.log
    print_progress
done




