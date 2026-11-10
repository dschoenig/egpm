#!/bin/bash
# Shared submission logic for rerun_block_*.sh. Not meant to be run directly.
#
# The calling script defines:
#   scripts       cluster scripts to submit in this block
#   prev_scripts  cluster scripts of the previous block (must have left the queue)
#   inter_files   stale 3_marginal_egp.R .inter.rds files to remove before submitting
#
# Every check runs before anything is submitted, so a block is either submitted
# completely or not at all. Job names default to the script file name, which is
# how jobs are matched against the queue.

# Per-user limit on queued jobs, where each array task counts as one job.
MAX_JOBS="${MAX_JOBS:-1000}"

n_tasks() {
  local range
  range=$(grep -oP '(?<=^#SBATCH --array=)[0-9]+-[0-9]+' "$1" || true)
  if [ -z "$range" ]; then
    echo 1
  else
    echo $(( ${range#*-} - ${range%-*} + 1 ))
  fi
}

submit_block() {
  local s f queued_names queued_tasks needed=0

  for s in "${scripts[@]}"; do
    [ -f "$s" ] || { echo "Missing script: $s" >&2; exit 1; }
    needed=$(( needed + $(n_tasks "$s") ))
  done

  queued_names=$(squeue -u "$USER" -h -o "%j")

  for s in "${prev_scripts[@]}"; do
    if grep -qxF "$s" <<< "$queued_names"; then
      echo "Previous block still in the queue ($s). Wait until it has finished." >&2
      exit 1
    fi
  done

  for s in "${scripts[@]}"; do
    if grep -qxF "$s" <<< "$queued_names"; then
      echo "Already in the queue: $s. Not submitting this block twice." >&2
      exit 1
    fi
  done

  queued_tasks=$(squeue -u "$USER" -h -r | wc -l)
  if [ $(( queued_tasks + needed )) -gt "$MAX_JOBS" ]; then
    echo "This block needs $needed jobs, $queued_tasks are queued, limit is $MAX_JOBS." >&2
    echo "Wait until more jobs have finished." >&2
    exit 1
  fi

  if [ "${#prev_scripts[@]}" -gt 0 ]; then
    echo "Reminder: check that the previous block finished without failures, e.g."
    echo "  sacct -u $USER -S <date> --format=JobName%45,State | grep -v COMPLETED"
  fi

  # 3_marginal_egp.R resumes from an existing .inter.rds file, which would mix
  # old results into the rerun.
  for f in "${inter_files[@]}"; do
    if [ -f "$f" ]; then
      echo "Removing $f"
      rm -f "$f"
    fi
  done

  echo "Submitting ${#scripts[@]} scripts ($needed jobs) …"
  for s in "${scripts[@]}"; do
    echo "Submitting $s"
    sbatch "$s"
  done
}
