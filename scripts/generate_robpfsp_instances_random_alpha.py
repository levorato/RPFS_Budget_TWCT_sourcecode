'''
    Script to generate Robust PFSP instance files with random alpha,
    i.e. random proc. time oscillation, for each processing time operation.
    For each original Ying instance, generates 2 new instance files, for both
    makespan and Weighted Completion Time (WCT) objectives.
    Alpha is randomly drawn from a uniform distribution in the interval
    [0, MAX_ALPHA] where MAX_ALPHA can be equal to different values
    (e.g. 30, 50, 100, 200, 400). 

    Usage: python3 generate_robpfsp_instances_random_alpha.py --output ./output
                ../instances/robust/ying/rob-pfsp-wct/10jobs
'''

import sys, getopt
import csv
import glob
import os
import os.path
import argparse
import time
import numpy as np
seed = 42
MAX_JOBS = 400  # the maximum number of jobs allowed
MAX_MACHINES = 400

def main(argv):
    csv.field_size_limit(1000000000)

    parser = argparse.ArgumentParser(description='Generate Rob-PFSP instance files with random alpha.')
    parser.add_argument('folders', nargs='+',
                        help='the folders containing the instance files')
    parser.add_argument('--filefilter', default='.txt', required=False,
                        help='the filename for instance files (default: .txt)')
    parser.add_argument('--output', required=True,
                        help='the output folder for writing the converted instance files.')
    args = parser.parse_args()
    folders = args.folders
    filter = args.filefilter
    output = args.output

    args = parser.parse_args()
    print('Input folders are ', folders)
    print('File filter is ', filter)
    print('Output folder is ', output)
    processInstanceFiles(folders, filter, output)


def read_robust_input_file(x):
    file = open(x, 'r')
    m = 0
    n = 0
    i = j = 0
    P_bar = np.zeros((MAX_MACHINES, MAX_JOBS))
    P_hat = np.zeros((MAX_MACHINES, MAX_JOBS))
    w = np.zeros((MAX_JOBS))
    step = 0
    for ln in file.readlines():
        line = ln.strip()
        if len(line) > 0:  # ignore empty lines
            if "nJobs" in line:
                print("Step 1")
                step = 1
                continue
            #end
            if "Weights" in line:
                step = 2
                j = 0
                continue
            #end
            if "P_bar" in line:
                step = 3
                j = 0
                continue
            #end
            if "P_hat" in line:
                step = 4
                j = 0
                continue
            #end
            line_array = np.fromstring(line, dtype=float, sep=' ')
            # println("Read header line $(j): $(line) ; $(line_array)")
            # determine the number of machines
            if step == 1:  # n m
                n = int(line_array[0])
                m = int(line_array[1])
                print("The number of jobs is $(n).")
                print("The number of machines is $(m).")
                P_bar = np.zeros((m, n))
                P_hat = np.zeros((m, n))
                w = np.zeros((n))
                continue
            elif step == 2:  # job weights
                w[j] = line_array[0]
                j += 1
            elif step == 3:  # P_bar
                for i in range(0, m):
                    P_bar[i, j] = line_array[i]
                #end
                j += 1
            elif step == 4:  # P_hat
                for i in range(0, m):
                    P_hat[i, j] = line_array[i]
                #end
                j += 1
            #end
        #end
    #end
    file.close()
    # println("M = $(m), N = $(n), P_bar = $(P_bar), P_hat = $(P_hat), w = $(w).")
    print("Flow shop robust instance file read OK : $(x).")
    return m, n, P_bar, P_hat, w


def processInstanceFiles(folders, filter, outputpath):
    np.random.seed(seed)
    if not os.path.exists(outputpath):
        os.makedirs(outputpath)
    for folder in folders:
        print("Processing folder " + ''.join(folder))
        for root, subFolders, files in os.walk(folder):
            # sort dirs and files
            subFolders.sort()
            files.sort()
            if len(files):
                file_list = [f for f in files if filter in f]
                for file in file_list:
                    idx = file.find("_wct_inputs.txt")
                    alpha = file[idx-2:idx]
                    print("alpha = ", alpha)
                    if alpha == '10' or alpha.find('_') >= 0:  # we will generate the instances based on alpha=10 instances
                        # measure elapsed time during instance generation
                        start = time.time()
                        file = os.path.join(root, file)
                        filename = file
                        print("Processing file " + filename)
                        m, n, P_bar, P_hat, w = read_robust_input_file(filename)
                        print('Instance size: n = {}, m = {}'.format(n, m))
                        print(P_bar, P_hat, w)
                        # randomly generate matrix of alpha values
                        # uniform distribution in the interval [0, 1)
                        for alpha_dev in [0.3, 0.5, 1.0, 2.0, 4.0]:
                            print("alpha_dev = " + str(alpha_dev))
                            alpha_values = np.random.rand(m, n)
                            print("alpha_values", alpha_values)
                            for exportweights in ['True', 'False']:
                                # output file
                                filename = file[file.rfind(os.sep) + 1:file.find("_")]
                                filename += '_' + str(n) + '_' + str(m) + '_R' + '%d' % (alpha_dev*100)
                                if exportweights == 'True':
                                    filename += '_wct_inputs.txt'
                                else:
                                    filename += '_cmax_inputs.txt'
                                output_file = open(os.path.join(outputpath, filename), 'w')
                                try:
                                    output_file.write("# nJobs | nMachines\r\n")
                                    output_file.write("{0} {1}\r\n".format(n, m))
                                    if exportweights == 'True':
                                        output_file.write("# Job Weights\r\n")
                                        np.savetxt(output_file, w, fmt='%d')
                                    mlist = []
                                    for x in range(1, m):
                                        mlist.append(' m{} |'.format(x))
                                    mlist.append(' m{}\r\n'.format(m))
                                    output_file.write('# P_bar :' + ''.join(mlist))
                                    # Ying instance files always have integer processing times
                                    np.savetxt(output_file, np.transpose(P_bar), fmt='%.2f')
                                    # Generate processing time deviations, based on alpha value
                                    matrix_dev = np.zeros((m, n))
                                    for i in range(0, m):
                                        for j in range(0, n):
                                            matrix_dev[i, j] = alpha_dev * P_bar[i, j] * alpha_values[i, j]
                                    # processing time deviations may be float
                                    output_file.write('# P_hat :' + ''.join(mlist))
                                    np.savetxt(output_file, np.transpose(matrix_dev), fmt='%.2f')
                                finally:
                                    output_file.close()
                                print("\nCreated output instance file {0}".format(filename))
                            # end for
                        # end for
                        end = time.time()
                        elapsed = end - start
                        print("Instance generation took {0:.2f} seconds.".format(elapsed))
                        #break
                    # end if
                # end loop
                # process last file
        print("\nDone.\n")


if __name__ == "__main__":
    main(sys.argv[1:])
