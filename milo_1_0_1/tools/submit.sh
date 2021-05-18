#!/bin/bash
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=24
#SBATCH --mem=24G
#SBATCH -t 16:00:00
#SBATCH -C 'avx2'

export JOB_NAME=

export TEMPORARY_DIR=/tmp/$SLURM_JOB_ID
export JOB_SOURCE_DIR="$SLURM_SUBMIT_DIR"

function cleanup
{
  echo "---"
  echo "Starting clean up"

  for file in "$TEMPORARY_DIR"/*.{out,xyz,txt}; do
    [ -e $file ] || continue
    cp -vp $file "$JOB_SOURCE_DIR"
  done

  mkdir -p "$JOB_SOURCE_DIR/com_log_files/"
  echo "Archiving .com and .log files"
  tar -cf ${JOB_NAME}_${SLURM_JOB_ID}_com_log_files.tar *.{com,log}
  cp -v ${JOB_NAME}_${SLURM_JOB_ID}_com_log_files.tar "${JOB_SOURCE_DIR}/com_log_files/"

  cd "$JOB_SOURCE_DIR"
  rm -fr "$TEMPORARY_DIR"

  echo "Clean up finished at `date`"
}

echo "---"
echo "Milo SLURM Diagnostic Information"
echo "---"
echo "Start date and time: `date`"
echo "---"
echo "Nodes assigned:"
cat `/fslapps/fslutils/generate_pbs_nodefile`
echo "---"
echo "Job Source Directory: $JOB_SOURCE_DIR"
echo "Temporary Directory: $TEMPORARY_DIR"
mkdir $TEMPORARY_DIR
cp -v "$JOB_SOURCE_DIR/$JOB_NAME.in" "$TEMPORARY_DIR"
cd "$TEMPORARY_DIR"
echo "---"
echo "Starting Milo run"

module load python/3.8
#module load g09
module load g16
python -m milo_1_0_1 < "$JOB_NAME.in" > "${JOB_NAME}_${SLURM_JOB_ID}.out" &
pid=$!
# Associate the function "cleanup" with the TERM signal, which is usually
# how jobs get killed
trap "kill $pid; cleanup; exit 1" TERM SIGTERM KILL SIGKILL
wait $pid

cleanup
