#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create an input file with opposite phase for each output file in the directory.

This script will create and input file (.in) and a SLURM submission script
(.sh) for each .out file in the directory this script is called from. The input
files will have opposite phase, meaning they will have push_apart instead of
bring_together or bring_together instead of push_apart. This is only helpful
for trajectories that are initiated from a transition state. This only works if
the .out files had phase specified as push_apart or bring_together, not as
random.

You can use command line options to determine max_steps for the input file and
memory, number of CPUs, and time limit for the submission script.
"""

import argparse
import os


def main():
    """Serve as main."""
    parser = argparse.ArgumentParser(description="Make Milo input files "
                                     "with opposite phase from .out files.\n")

    parser.add_argument("-t", "--time", required=True,
                        help="Amount of time the submission script should "
                        "request, formatted as HH:MM:SS.")
    parser.add_argument("-m", "--memory", type=int, required=True,
                        help="Amount of memory, in GB, the submission script "
                        "should request. Also modifies the input file to "
                        "request 1 fewer GBs than the submission script. ")
    parser.add_argument("-p", "--processors", type=int, required=True,
                        help="Number of processors the submission script "
                        "should request. Also modifies the input file to match"
                        " the submission script.")
    parser.add_argument("-s", "--max_steps", default=None,
                        help="The maximum number of trajectory steps. Either "
                        "an integer or 'no_limit'. Default: leave the same as"
                        "the .out file.")
    parser.add_argument("--no_script", action="store_true",
                        help="Stops submission scripts from being output.")
    parser.add_argument("--gaussian09", action="store_true",
                        help="Uses 'module load g09' in the submission script,"
                        " instead of 'module load g16' (the default).")

    args = parser.parse_args()

    out_files = [f for f in os.listdir('.') if os.path.isfile(f) and
                 f.endswith('.out')]

    for out_file in out_files:
        with open(out_file, mode="r") as out_reader:
            input = list()
            random_seed = None
            job_name = out_file.replace(".out", "") + "_rev"

            in_input_section = False
            in_random_seed_section = False
            for line in out_reader:
                if "### Input File ----------------------------------" in line:
                    in_input_section = True
                elif "### Default Parameters Being Used -------------" in line:
                    in_input_section = False
                elif in_input_section:
                    input.append(line)
                elif "### Random Seed -------------------------------" in line:
                    in_random_seed_section = True
                elif in_random_seed_section:
                    random_seed = line.strip()
                    break

            make_input_file(job_name, input, random_seed, args.max_steps,
                            args.memory, args.processors)
            if not args.no_script:
                make_submission_script(job_name, args.time, args.memory,
                                       args.processors, args.gaussian09)


def make_input_file(job_name, old_input, random_seed, max_steps, memory,
                    procs):
    """Make an input file that matches, but with phase switched."""
    new_input = list()

    in_job_section = False
    random_seed_set = False
    for line in old_input:
        if "$job" in line.casefold():
            in_job_section = True
        elif in_job_section and "$end" in line.casefold():
            in_job_section = False
            if not random_seed_set:
                new_input.append(f"    random_seed             {random_seed}"
                                 f"\n")
                random_seed_set = True
        elif in_job_section and "random_seed" in line.casefold():
            line = f"    random_seed             {random_seed}\n"
            random_seed_set = True
        elif (in_job_section and max_steps is not None and
                "max_steps" in line.casefold()):
            line = f"    max_steps               {max_steps}\n"
        elif in_job_section and "phase" in line.casefold():
            if "bring_together" in line.casefold():
                line = line.replace("bring_together", "push_apart")
            elif "push_apart" in line.casefold():
                line = line.replace("push_apart", "bring_together")
            else:
                print(f"{job_name[:-4]}.out -- Warning: This script only "
                      "works as intended if 'push_apart' or 'bring_together' "
                      "was specified for phase in the original input file.")
        elif in_job_section and "memory" in line.casefold():
            line = f"    memory                  {memory - 1}\n"
        elif in_job_section and "processors" in line.casefold():
            line = f"    processors              {procs}\n"
        new_input.append(line)

    with open(job_name + ".in", "w") as file:
        for line in new_input:
            file.write(line)


def make_submission_script(job_name, time_string, memory, procs, g09_flag):
    """Make a submission script given filename and time/job parameters."""
    with open(job_name + ".sh", "w") as file:
        file.write('#!/bin/bash\n')
        file.write(f'#SBATCH --nodes=1 --ntasks=1 --cpus-per-task={procs}\n')
        file.write(f'#SBATCH --mem={memory}G\n')
        file.write(f'#SBATCH -t {time_string}\n')
        file.write("#SBATCH -C 'avx2'\n")
        file.write('\n')
        file.write(f'export JOB_NAME={job_name}\n')
        file.write('\n')
        file.write('export TEMPORARY_DIR=/tmp/$SLURM_JOB_ID\n')
        file.write('export JOB_SOURCE_DIR="$SLURM_SUBMIT_DIR"\n')
        file.write('\n')
        file.write('function cleanup\n')
        file.write('{\n')
        file.write('  echo "---"\n')
        file.write('  echo "Starting clean up"\n')
        file.write('\n')
        file.write('  for file in "$TEMPORARY_DIR"/*.{out,xyz,txt}; do\n')
        file.write('    [ -e $file ] || continue\n')
        file.write('    cp -vp $file "$JOB_SOURCE_DIR"\n')
        file.write('  done\n')
        file.write('\n')
        file.write('mkdir -p "$JOB_SOURCE_DIR/com_log_files/"\n')
        file.write('echo "Archiving .com and .log files"\n')
        file.write('tar -cf ${JOB_NAME}_${SLURM_JOB_ID}_com_log_files.tar '
                   '*.{com,log}\n')
        file.write('cp -v ${JOB_NAME}_${SLURM_JOB_ID}_com_log_files.tar '
                   '"${JOB_SOURCE_DIR}/com_log_files/"\n')
        file.write('\n')
        file.write('  cd "$JOB_SOURCE_DIR"\n')
        file.write('  rm -fr "$TEMPORARY_DIR"\n')
        file.write('\n')
        file.write('  echo "Clean up finished at `date`"\n')
        file.write('}\n')
        file.write('\n')
        file.write('echo "---"\n')
        file.write('echo "Milo SLURM Diagnostic Information"\n')
        file.write('echo "---"\n')
        file.write('echo "Start date and time: `date`"\n')
        file.write('echo "---"\n')
        file.write('echo "Nodes assigned:"\n')
        file.write('cat `/fslapps/fslutils/generate_pbs_nodefile`\n')
        file.write('echo "---"\n')
        file.write('echo "Job Source Directory: $JOB_SOURCE_DIR"\n')
        file.write('echo "Temporary Directory: $TEMPORARY_DIR"\n')
        file.write('mkdir $TEMPORARY_DIR\n')
        file.write('cp -v "$JOB_SOURCE_DIR/$JOB_NAME.in" "$TEMPORARY_DIR"\n')
        file.write('cd "$TEMPORARY_DIR"\n')
        file.write('echo "---"\n')
        file.write('echo "Starting Milo run"\n')
        file.write('\n')
        file.write('module load python/3.8\n')
        if g09_flag:
            file.write('module load g09\n')
        else:
            file.write('module load g16\n')
        file.write('python -m milo_1_0_2 < "$JOB_NAME.in" > ')
        file.write('"$JOB_NAME.out" &\n')
        file.write('pid=$!\n')
        file.write('# Associate the function "cleanup" with the TERM signal, '
                   'which is usually\n')
        file.write('# how jobs get killed\n')
        file.write('trap "kill $pid; cleanup; exit 1" TERM SIGTERM KILL '
                   'SIGKILL\n')
        file.write('wait $pid\n')
        file.write('\n')
        file.write('cleanup\n')
        file.write('\n')


if __name__ == "__main__":
    main()
