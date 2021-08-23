#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create multiple input files with different seeds to form an ensemble.

This script can setup input files in two ways:
 - from a single input file, corresponding to only one direction, or
 - from two input files, one forward and one backward.

The script will also create SLURM submission scripts, unless turned off with
the --no_script flag.
"""

import argparse
from random import randrange


def main():
    """Serve as main."""
    parser = argparse.ArgumentParser(description="Make Milo input files for an"
                                     " ensemble of trajectories.")

    parser.add_argument("-n", "--number_trajectories", type=int, required=True,
                        help="The number of trajectories to make input files "
                        "for.")

    parser.add_argument("-f", "--forward", type=argparse.FileType("r"),
                        required=True, help="A template input file for the "
                        "forward direction.")
    parser.add_argument("-b", "--backward", type=argparse.FileType("r"),
                        default=None, help="A template input file for the "
                        "backward direction.")

    parser.add_argument("-t", "--time", default="time_placeholder",
                        help="Amount of time the submission script should "
                        "request for forward trajectoryies, formatted as "
                        "HH:MM:SS.")
    parser.add_argument("--time_backward", default=None,
                        help="Amount of time the submission script should "
                        "request for backward trajectories, formatted as "
                        "HH:MM:SS. If not provided, defaults to same time as "
                        "forward trajectories.")

    parser.add_argument("--no_script", action="store_true",
                        help="Stops submission scripts from being output.")

    args = parser.parse_args()

    gaussian_string = "g16"

    forward_template = list()
    in_job_section = False
    random_seed_set = False
    forward_memory = "memory_placeholder"
    forward_procs = "processors_placeholder"
    for line in args.forward:
        if "$job" in line.casefold():
            in_job_section = True
        elif in_job_section and "$end" in line.casefold():
            in_job_section = False
            if not random_seed_set:
                forward_template.append("    random_seed             "
                                        "random_seed_placeholder\n")
                random_seed_set = True
        elif in_job_section and "random_seed" in line.casefold():
            line = "    random_seed             random_seed_placeholder\n"
            random_seed_set = True
        elif in_job_section and "memory" in line.casefold():
            forward_memory = int(line.split()[1]) + 1
        elif in_job_section and "processors" in line.casefold():
            forward_procs = int(line.split()[1])
        elif in_job_section and "program" in line.casefold() and \
                "gaussian09" in line.casefold():
            gaussian_string = "g09"
        forward_template.append(line)

    if not args.no_script and forward_memory == "memory_placeholder":
        print("Warning: Could not find memory specification in forward "
              "template file. Replace 'memory_placeholder' in submission "
              "script before submitting jobs.")
    if not args.no_script and forward_procs == "processors_placeholder":
        print("Warning: Could not find processors specification in forward "
              "template file. Replace 'processors_placeholder' in submission "
              "script before submitting jobs.")

    if args.backward is not None:
        backward_template = list()
        in_job_section = False
        random_seed_set = False
        backward_memory = "memory_placeholder"
        backward_procs = "processors_placeholder"
        for line in args.backward:
            if "$job" in line.casefold():
                in_job_section = True
            elif in_job_section and "$end" in line.casefold():
                in_job_section = False
                if not random_seed_set:
                    backward_template.append("    random_seed             "
                                            "random_seed_placeholder\n")
                    random_seed_set = True
            elif in_job_section and "random_seed" in line.casefold():
                line = "    random_seed             random_seed_placeholder\n"
                random_seed_set = True
            elif in_job_section and "memory" in line.casefold():
                backward_memory = int(line.split()[1]) + 1
            elif in_job_section and "processors" in line.casefold():
                backward_procs = int(line.split()[1])
            backward_template.append(line)
        if not args.no_script and forward_memory == "memory_placeholder":
            print("Warning: Could not find memory specification in backward "
                  "template file. Replace 'memory_placeholder' in submission "
                  "script before submitting jobs.")
        if not args.no_script and forward_procs == "processors_placeholder":
            print("Warning: Could not find processors specification in forward"
                  " template file. Replace 'processors_placeholder' in "
                  "submission script before submitting jobs.")

    if not args.no_script and args.time == "time_placeholder":
        print("Warning: Time not specified. Replace 'time_placeholder' in "
              "submission script before submitting jobs.")

    input_file_basename = args.forward.name.split(".")[0]

    for _ in range(args.number_trajectories):
        random_seed = randrange(10000000000, 99999999999)
        if args.backward is not None:
            forward_name = (input_file_basename + "_" +
                           str(random_seed) + "_fwd.in")
            make_input_file(forward_name, forward_template, random_seed)
            backward_name = (input_file_basename + "_" +
                           str(random_seed) + "_bwd.in")
            make_input_file(backward_name, backward_template, random_seed)
            if not args.no_script:
                make_submission_script(forward_name.replace(".in", ".sh"),
                    args.time, forward_memory, forward_procs, gaussian_string)
                make_submission_script(backward_name.replace(".in", ".sh"),
                    args.time_backward if args.time_backward is not None else
                    args.time, backward_memory, backward_procs,
                    gaussian_string)
        else:
            name = (input_file_basename + "_" +
                    str(random_seed) + ".in")
            make_input_file(name, forward_template, random_seed)
            if not args.no_script:
                make_submission_script(name.replace(".in", ".sh"), args.time,
                                       forward_memory, forward_procs,
                                       gaussian_string)


def make_input_file(file_name, template, random_seed):
    """Make an input file that matches, but with phase switched."""
    with open(file_name, "w") as file:
        for line in template:
            if "random_seed_placeholder" in line:
                line = line.replace("random_seed_placeholder",
                                    str(random_seed))
            file.write(line)


def make_submission_script(file_name, time_string, memory, procs,
                           gaussian_string):
    """Make a submission script given filename and time/job parameters."""
    with open(file_name, "w") as file:
        file.write('#!/bin/bash\n')
        file.write(f'#SBATCH --nodes=1 --ntasks=1 --cpus-per-task={procs}\n')
        file.write(f'#SBATCH --mem={memory}G\n')
        file.write(f'#SBATCH -t {time_string}\n')
        file.write("#SBATCH -C 'avx2'\n")
        file.write('\n')
        file.write(f'export JOB_NAME={file_name.split(".")[0]}\n')
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
        file.write(f'module load {gaussian_string}\n')
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
