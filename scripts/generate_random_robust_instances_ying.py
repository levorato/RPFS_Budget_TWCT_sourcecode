# =============================================================================
#   generate_random_robust_instances_ying.py
# =============================================================================
#   Python script to generate random instances for the Robust
#    Permutation Flow Shop Problem, processing time intervals, with
#    an arbitraty number of jobs and machines.
#   Based on Ying's (2005) recipe for generating random instances.
#   Reference:
#   @article{Ying2015,
#      author = {Ying, Kuo Ching},
#      doi = {10.1057/jors.2014.100},
#      issn = {14769360},
#      journal = {Journal of the Operational Research Society},
#      number = {9},
#      pages = {1413--1425},
#      title = {{Scheduling the two-machine flowshop to hedge against processing time uncertainty}},
#      volume = {66},
#      year = {2015}
#   }
# =============================================================================

import sys, getopt
import glob
import os
import os.path
import argparse
import time
import math
import numpy as np
import numpy.random
import re
import datetime


def main(argv):
    parser = argparse.ArgumentParser(description='Generate Robust PFSP \
        random instance files.')
    parser.add_argument('-n', type=int, required=True, help="number of jobs")
    parser.add_argument('-m', type=int, required=True, help="number of machines")
    parser.add_argument('--outputfolder', required=True,
                        help='output folder to write the new instances.')
    args = parser.parse_args()
    m = args.m
    n = args.n
    output_folder = args.outputfolder
    args = parser.parse_args()
    print('Generating Robust PFSP random instance files for m = {0} and \
        n = {1}...'.format(m, n))
    print('Output folder is ', output_folder)
    processInstanceFiles(m, n, output_folder)


# The expected / nominal processing time p_bar_{i,j} (j = 1,.., n; i = 1,.., m)
# is an integer to be generated from the uniform distribution [10, 50].
# The largest processing time deviation (p_hat{i,j}) is set as a ratio of the
# expected processing time (i.e. p_hat{i,j} = \alpha \times p_bar_{i,j}), where
# \alpha = 10, 20, 30, 40 and 50 \%.
# Ten instances are generated for each combination of m, n and \alpha,
# resulting in a total of 10 \times \alpha test instances.
# With particular interest in examining the impact of the uncertain budget
# parameters \Gamma_1, \Gamma_2, ..., \Gamma_m on scheduling performance,
# five ratios (i.e. 20, 40, 60, 80 and 100%) of jobs with uncertain processing
# times in M_1, M_2, ..., M_m, respectively, were evaluated for each test instance.
# (e.g. number of jobs n = 10, 20, 50,100,150 and 200)
def generate_robust_instance(m, n, output_dir):
    n_str = '{0:03d}'.format(n)
    m_str = '{0:03d}'.format(m)
    for gen_count in range(10):
        gen_count_str = '{0:02d}'.format(gen_count + 1)
        # Generate random matrix using uniform distribution [10, 50]
        P_bar = np.random.randint(10, 51, size=(m, n))
        for alpha in [0.1, 0.2, 0.3, 0.4, 0.5]:
            alpha_str = '{0:02d}'.format(int(alpha * 100))
            n_str = '{0:03d}'.format(n)
            output_filepath = os.path.join(output_dir, alpha_str + '%')
            output_filename = "RB" + n_str + alpha_str + gen_count_str + '_' + n_str + "_" + m_str + "_" \
                + alpha_str + "%_" + gen_count_str + ".txt"
            if not os.path.exists(output_filepath):
                os.makedirs(output_filepath)
            output_filepath = os.path.join(output_filepath, output_filename)
            print("Generating output file " + output_filepath)
            with open(output_filepath, 'w') as output_file:
                # writes the header info
                output_file.write("! n m alpha seq\r\n")
                output_file.write("{0} {1} {2} {3}\r\n".format(n, m, alpha, gen_count+1))
                P_hat = np.zeros((m,n), dtype=float)
                P_hat = P_bar * alpha
                # Export to text file: is each line a job, each column is a machine
                for i in range(n):  # For each job i
                    line_str = list()
                    for r in range(m):  # For each machine r
                        line_str.append('{0:.2f}\t'.format(P_bar[r][i]))
                    for r in range(m):  # For each machine r
                        line_str.append('{0:.2f}'.format(P_hat[r][i]))
                        if r < m - 1:
                            line_str.append('\t')
                    output_file.write('{0}\r\n'.format(''.join(line_str)))
            print("Created output robust instance file {0}\n".format(output_filepath))
        # end for alpha
    # end for gen_count

def processInstanceFiles(m, n, output_folder):
    # measure elapsed time during instance generation
    start = time.time()
    datetime_str='{date:%Y_%m_%d-%H_%M_%S}'.format(date=datetime.datetime.now())
    output_folder = os.path.join(output_folder, datetime_str)
    output_dir = os.path.join(output_folder, '{0}x{1}'.format(n, m))
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    generate_robust_instance(m, n, output_dir)
    end = time.time()
    elapsed = end - start
    print("Instance generation took {0:.2f} seconds.".format(elapsed))
    print("\nDone.\n")


if __name__ == "__main__":
    main(sys.argv[1:])
