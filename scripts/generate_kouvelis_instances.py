#!/usr/bin/env python
# coding: utf-8

import sys, getopt
import csv
import os
import os.path
import argparse
import time
import numpy as np


class Range(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __eq__(self, other):
        return self.start <= other <= self.end

    def __str__(self):
        return "range(" + str(self.start) + ", " + str(self.end) + ")"


def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def main(argv):
    csv.field_size_limit(1000000000)
    print("Kouvelis et al.(2000) instance generator")
    print("=========================================\n")

    parser = argparse.ArgumentParser(description='Generate instance file (.in) for the 2-machine robust flowshop problem.')
    parser.add_argument("--discrete", type=str2bool, nargs='?',
                        const=True, default=False,
                        help="Generate discrete scenario instances.")
    args = parser.parse_args()

    # measure elapsed time during instance generation
    start = time.time()
    np.random.seed(99)
    if args.discrete:
        print('Generating discrete scenarios instance files...')
        generateAllInstanceFiles_DiscreteCase()
    else:
        print('Generating continuous intervals instance files...')
        generateAllInstanceFiles_IntervalCase()

    end = time.time()
    elapsed = end - start
    print("Instance generation took {0:.2f} seconds.".format(elapsed))
    print('Done.')


def generateAllInstanceFiles_IntervalCase():
    num_replicas = 10
    n_list = [9, 12, 15]  # number of jobs
    alpha1_list = [0.2, 0.6, 1.0]
    alpha2_list = [0.2, 0.6, 1.0]
    beta_list = [[1.0, 1.0], [1.2, 1.0], [1.0, 1.2]]
    m = 2
    current_directory = '../instances/robust/kouvelis/test/interval'
    if not os.path.exists(current_directory):
        os.makedirs(current_directory)

    # Total of 810 test problems
    instance_count = 0
    for n in n_list:  # number of jobs
        for alpha1 in alpha1_list:
            for alpha2 in alpha2_list:
                for beta in beta_list:  # pair of values (beta_1, beta_2)
                    for r in range(num_replicas):  # one full instance is generated for each replica
                        # create a new instance file for combination n x alpha1 x alpha2 x beta x replica_number
                        filename = 'm' + str(m) + 'n' + str(n) + 'alphaone' + str(alpha1) + 'alphatwo' + str(alpha2) \
                                   + 'beta' + str(beta).replace(' ', '') \
                                   + 'r' + str(r + 1) + ".txt"
                        output_file = open(os.path.join(current_directory, filename), 'w')
                        matrices = generateProcessingTimeMatrixInterval(m, n, alpha1, alpha2, beta)
                        try:
                            output_file.write("# m n rep Lambda alpha1 alpha2 beta\r\n")
                            output_file.write("{0} {1} {2} {3} {4} {5} {6}\r\n".format(m, n, r + 1, 0, alpha1, alpha2,
                                                                                       str(beta).replace(' ', '')))
                            for matrix in matrices:
                                output_file.write("{0}\r\n".format(matrix))
                            print("Created output instance file {0}".format(filename))
                            instance_count += 1
                        finally:
                            output_file.close()
                    # end for r
                # end for beta
            # end for alpha 2
        # end for alpha 1
    # end for n
    print('Generated ' + str(instance_count) + ' processing interval scenario instances.')


def generateAllInstanceFiles_DiscreteCase():
    num_replicas = 10
    n_list = [9, 12, 15]  # number of jobs
    l_list = [4, 8, 12]  # number of scenarios
    alpha_list = [0.2, 0.6, 1.0]
    beta_list = [[1.0, 1.0], [1.2, 1.0], [1.0, 1.2]]
    m = 2
    current_directory = '../instances/robust/kouvelis/test/discrete'
    if not os.path.exists(current_directory):
        os.makedirs(current_directory)

    # Total of 810 test problems
    instance_count = 0
    for n in n_list:  # number of jobs
        for alpha in alpha_list:
            for beta in beta_list:  # pair of values (beta_1, beta_2)
                for r in range(num_replicas):  # one full instance is generated for each replica
                    for lambda_size in l_list:
                        # create a new instance file for combination n x alpha x beta x lambda_size
                        filename = 'm' + str(m) + 'n' + str(n) + 'alpha' + str(alpha) + 'beta' \
                                   + str(beta).replace(' ', '') \
                                   + 's' + str(lambda_size) + 'r' + str(r + 1) + ".txt"
                        output_file = open(os.path.join(current_directory, filename), 'w')
                        matrices = []
                        for l in range(lambda_size):  # for each scenario l : gen list of scenarios {1, ..., l}
                            matrices.append(generateProcessingTimeMatrixDiscrete(m, n, alpha, beta))
                        # end for l
                        try:
                            output_file.write("# m n rep Lambda alpha beta\r\n")
                            output_file.write("{0} {1} {2} {3} {4} {5}\r\n".format(m, n, r + 1, lambda_size, alpha,
                                                                                   str(beta).replace(' ', '')))
                            scen_count = 1
                            for matrix in matrices:
                                output_file.write("Scenario {0}\n{1}\r\n".format(scen_count, matrix))
                                scen_count += 1
                            print("Created output instance file {0}".format(filename))
                            instance_count += 1
                        finally:
                            output_file.close()
                    # end for lambda_size
                # end for r
            # end for beta
        # end for alpha
    # end for n
    print('Generated ' + str(instance_count) + ' discrete scenario instances.')


def generateProcessingTimeMatrixInterval(m, n, alpha1, alpha2, beta):
    # sample the lower range of processing time (p_low_ij) from the uniform distribution of integers on the interval
    # p_low => uniform[10*beta[j], (10+40*alpha1)*beta[j]]
    # the generated p_low matrix of random numbers has dimensions (n x m)
    p_low_str = str('')
    p_low = []
    for j in range(m):  # for each machine j
        r_vector = np.random.uniform(10*beta[j], (10+40*alpha1)*beta[j] + 1, size=n)
        r_vector = np.around(r_vector, 3)
        r_str = str(r_vector)
        r_str = r_str.replace('[', '')
        r_str = r_str.replace(']', '')
        r_str = r_str.replace(',', '')
        r_str = r_str.replace('\r', '')
        r_str = r_str.replace('\n', '')
        p_low_str += r_str + "\n"
        p_low.append(r_vector)
    # end for j
    # sample the upper end of processing time (p_high_ij) from the uniform distribution of integers on the interval
    # p_high => uniform[p_low(i, j), p_low(i, j) * (1 + alpha2)]
    # the generated p_high matrix of random numbers has dimensions (n x m)
    p_high_str = str('')
    for j in range(m):  # for each machine j
        for i in range(n):  # for each job i
            x = round(np.random.uniform(p_low[j][i], (p_low[j][i] * (1 + alpha2)) + 1, size=1)[0], 3)
            p_high_str += str(x) + ' '
        # end for i
        p_high_str += '\n'
    # end for j
    return [p_low_str, p_high_str]


def generateProcessingTimeMatrixDiscrete(m, n, alpha, beta):
    # sample the the processing time from the uniform distribution of integers on the interval
    # [10*beta[j], (10+40*alpha)*beta[j]]
    # the generated matrix of random numbers has dimensions (n x m)
    matrix = str('')
    for j in range(m):  # for each machine j
        r_vector = np.random.uniform(10*beta[j], (10+40*alpha)*beta[j] + 1, size=n)
        r_vector = np.around(r_vector, 3)
        r_str = str(r_vector)
        r_str = r_str.replace('[', '')
        r_str = r_str.replace(']', '')
        r_str = r_str.replace(',', '')
        r_str = r_str.replace('\r', '')
        r_str = r_str.replace('\n', '')
        matrix += r_str + "\n"
    # end for j
    # print(matrix)
    return matrix


if __name__ == "__main__":
    main(sys.argv[1:])

