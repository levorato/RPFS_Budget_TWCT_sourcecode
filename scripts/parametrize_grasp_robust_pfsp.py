#!/usr/bin/env python
# coding: utf-8

import subprocess
import sys
import itertools


def format_list_param(param_string):
    param_string = param_string.replace("(", "")
    param_string = param_string.replace(")", "")
    param_string = param_string.replace(",", "")
    return param_string


def main(argv):

    t_values = [20, 40, 60, 80, 100]
    n_jobs = ['10 jobs']  #, '20 jobs', '50 jobs'] #, '100 jobs', '150 jobs', '200 jobs']
    perc_dev_list = ['10%']  #, '20%'] #, '30%', '40%', '50%']
    base_dir = '/home/ubuntu/RPFS_Budget_TWCT/instances/robust/ying/data'
    output_folder_param = '--output-folder=/home/ubuntu/pfsp-experiments/parametrize-robust-grasp'
    executable_path = './bin/pfsp'
    mh = 'grasp'
    total_vnd_neighborhoods = 8
    vnd_list = [x for x in range(1, total_vnd_neighborhoods + 1)]  # generate all possible VND neighborhood numbers
    error_count = 0
    for job_num in n_jobs:
        for perc_dev in perc_dev_list:
            for T1 in t_values:
                for T2 in t_values:
                    for first_improvement in ['true', 'false']:
                        for vnd_size in range(1, len(vnd_list) + 1):
                            # Generate all permutations of size 'vnd_size'
                            for vnd_permutation in itertools.permutations(vnd_list, vnd_size):
                                #print('Solving robust PFSP with T1=' + str(T1) + ', T2=' + str(T2) + ' and mh=' + str(mh) + '...\n')
                                rob_ub_type_param = '--rob-ub-type=' + mh
                                input_folder_param = '--input-file-dir=' + base_dir + '/' + job_num + '/' + perc_dev
                                budget_param = '--budget-T=' + str(T1) + ' ' + str(T2)
                                tf_param = '--grasp-tf=5'
                                if job_num == '50 jobs':
                                    tf_param = '--grasp-tf=30'
                                vnd_fi_param = '--first-improvement=' + str(first_improvement)
                                vnd_permutation = format_list_param(str(vnd_permutation))
                                vnd_perm_param = '--vnd-permutation=' + str(vnd_permutation)
                                vnd_size_param = '--vnd-size=' + str(vnd_size)
                                exec_list = [executable_path, '--time-limit=5', tf_param, '--model-version=104',
                                             rob_ub_type_param, budget_param, output_folder_param, input_folder_param,
                                             vnd_size_param, vnd_perm_param, vnd_fi_param]
                                print(str(exec_list))
                                run_count = 0
                                try:
                                    result = subprocess.check_call(exec_list, stderr=subprocess.PIPE)
                                    run_count += 1
                                except subprocess.CalledProcessError as e:
                                    error_count += 1
                                    print("ERROR (1) invoking process! " + str(e))
                                    try:
                                        result = subprocess.check_call(exec_list, stderr=subprocess.PIPE)
                                        run_count += 1
                                    except subprocess.CalledProcessError as e:
                                        error_count += 1
                                        print("ERROR (2) invoking process! " + str(e))
                                        result = subprocess.check_call(exec_list, stderr=subprocess.PIPE)
                                #print("Return: " + result)
                            # end for vnd_permutation

    print("EXPERIMENTS EXECUTION DONE. Total number of errors: " + str(error_count))


if __name__ == "__main__":
    main(sys.argv[1:])
