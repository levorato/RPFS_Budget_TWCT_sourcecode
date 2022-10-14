# =============================================================================
#   generate_taillard_robust_instances.py
#    Python script to generate Taillard-based instances for the Robust
#    Permutation Flow Shop Problem, processing time intervals.
# =============================================================================

import sys, getopt
import glob
import os
import os.path
import argparse
import time
import math
import numpy as np
import re
import datetime


def main(argv):
    parser = argparse.ArgumentParser(description='Convert Taillard instance \
        files (.in) to robust instance files.')
    parser.add_argument('folders', nargs='+',
                        help='the folders containing the instance files')
    parser.add_argument('--filefilter', default='.in', required=False,
                        help='the filename for instance files (default: .in)')
    parser.add_argument('--outputfolder', required=True,
                        help='the output folder to write the new instances.')
    args = parser.parse_args()
    folders = args.folders
    file_filter = args.filefilter
    output_folder = args.outputfolder
    args = parser.parse_args()
    print('Input folders are ', folders)
    print('File filter is ', file_filter)
    print('Output folder is ', output_folder)
    processInstanceFiles(folders, file_filter, output_folder)


def read_input_file(filepath):
    with open(filepath, 'r') as content_file:
        # reads the lines of the instance file
        row_count = 1
        m = n = r = 0
        matrix = np.zeros((1,1), dtype=float)
        for row in content_file:
            linestring = ''.join(row)
            column = re.split('\t| |\r|\n', linestring)
            column = [x for x in column if len(x) > 0]  # Remove empty elements
            if len(row) == 0 or ("!" in row):  # ignore comments and empty lines
                continue
            #print("Reading row {0}".format(column))
            if row_count == 1:  # line 1 => # of jobs (n) and # of machines (m)
                n = int(column[0])
                m = int(column[1])
                matrix = np.zeros((m,n), dtype=float)
            else: # valores correspondentes aos tempos das operacoes.
                for i in range(n):
                    matrix[r][i] = float(column[i])
                r += 1
            # end if
            row_count += 1
        # end read lines
        content_file.close()
        print("Successfully read Taillard input file.")
        return m, n, matrix
# end read_input_file()


# For each Taillard instance file (given by m, n, P_bar) in each group
# (number of jobs n = 10, 20, 50,100,150 and 200):
# The expected / nominal processing time p_bar_{i,j} (j = 1,.., n; i = 1,.., m)
# is the original processing time value (p_{i,j}) from the Taillard instance.
# The largest processing time deviation is set as a ratio of the
# expected processing time (i.e. p_hat{i,j} = \alpha \times p_bar_{i,j}), where
# \alpha = 10, 20, 30, 40 and 50 \%. One instance is generated for each value of
# \alpha, resulting in a total of 5 robust test instances for each original
# Taillard instance.
# With particular interest in examining the impact of the uncertain budget
# parameters \Gamma_1, \Gamma_2, ..., \Gamma_m on scheduling performance,
# five ratios (i.e. 20, 40, 60, 80 and 100%) of jobs with uncertain processing
# times in M_1, M_2, ..., M_m, respectively, were evaluated for each test instance.
def generate_robust_instance(m, n, P_bar, output_dir, taillard_instance_name):
    for alpha in [0.1, 0.2, 0.3, 0.4, 0.5]:
        alpha_str = '{0:d}%'.format(int(alpha * 100))
        output_filepath = os.path.join(output_dir, alpha_str)
        output_filename = taillard_instance_name + "_" + alpha_str + ".txt"
        if not os.path.exists(output_filepath):
            os.makedirs(output_filepath)
        output_filepath = os.path.join(output_filepath, output_filename)
        print("Generating output file " + output_filepath)
        with open(output_filepath, 'w') as output_file:
            # writes the header info
            output_file.write("! n m alpha\r\n")
            output_file.write("{0} {1} {2}\r\n".format(n, m, alpha))
            P_hat = np.zeros((m,n), dtype=float)
            P_hat = P_bar * alpha
            # Export to text file: is each line a job, each column is a machine
            for i in range(n):  # For each job i
                line_str = list()
                for r in range(m):  # For each machine r
                    line_str.append('{0:f}\t'.format(P_bar[r][i]))
                for r in range(m):  # For each machine r
                    line_str.append('{0:f}'.format(P_hat[r][i]))
                    if r < m - 1:
                        line_str.append('\t')
                output_file.write('{0}\r\n'.format(''.join(line_str)))
        print("Created output robust instance file {0}\n".format(output_filepath))
    # end for alpha


def processInstanceFiles(folders, file_filter, output_folder):
    # measure elapsed time during instance generation
    start = time.time()
    # datetime_str = '{date:%Y_%m_%d-%H_%M_%S}'.format( date=datetime.datetime.now() )
    output_folder = os.path.join(output_folder)
    for folder in folders:
        print("Processing folder " + ''.join(folder))
        for root, subFolders, files in os.walk(folder):
            # sort dirs and files
            subFolders.sort()
            files.sort()
            if len(files):
                file_list = [f for f in files if file_filter in f]
                for file in file_list:
                    # output_dir = os.path.join(output_folder, )
                    output_dir = output_folder
                    if not os.path.exists(output_dir):
                        os.makedirs(output_dir)
                    input_filepath = os.path.join(root, file)
                    print("Processing input file " + input_filepath)
                    m, n, P_bar = read_input_file(input_filepath)
                    taillard_instance_name = input_filepath[ \
                       input_filepath.rfind(os.sep)+1:input_filepath.rfind('.')]
                    generate_robust_instance(m, n, P_bar, output_dir, \
                        taillard_instance_name)
    # end for
    end = time.time()
    elapsed = end - start
    print("Instance generation took {0:.2f} seconds.".format(elapsed))
    print("\nDone.\n")


if __name__ == "__main__":
    main(sys.argv[1:])
