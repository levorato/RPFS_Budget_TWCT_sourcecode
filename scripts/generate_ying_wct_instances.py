'''
    Script to generate Ying's robust PFSP instance files for the Weighted Completion Time (WCT) objective.
'''

import sys, getopt
import csv
import glob
import os
import os.path
import argparse
import time
import numpy as np


def main(argv):
    csv.field_size_limit(1000000000)

    parser = argparse.ArgumentParser(description='Generate Ying instance files for Rob-PFSP-WCT problem.')
    parser.add_argument('folders', nargs='+',
                        help='the folders containing the instance files')
    parser.add_argument('--filefilter', default='.txt', required=False,
                        help='the filename for instance files (default: .txt)')
    parser.add_argument('--output', required=True,
                        help='the output folder for writing the converted instance files.')
    parser.add_argument('--exportweights', default='True', required=False,
                        help='export job weights to instance file?')
    args = parser.parse_args()
    folders = args.folders
    filter = args.filefilter
    output = args.output
    exportweights = args.exportweights

    args = parser.parse_args()
    print('Input folders are ', folders)
    print('File filter is ', filter)
    print('Output folder is ', output)
    print('Export job weights to instance file? ', exportweights)
    processInstanceFiles(folders, filter, output, exportweights)


def processInstanceFiles(folders, filter, outputpath, exportweights):
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
                    alpha = file[5:7]
                    if alpha == '10':  # we will manually generate all alpha values for a given instance
                        # measure elapsed time during instance generation
                        start = time.time()
                        file = os.path.join(root, file)
                        filename = file
                        print("Processing file " + filename)
                        matrix = np.loadtxt(file, skiprows=2)
                        n = int(np.size(matrix, 0))
                        m = int(np.size(matrix, 1) / 2)
                        print('Instance size: n = {}, m = {}'.format(n, m))
                        #print('Before: ' + str(matrix))
                        # resize matrix and remove the uncertain part
                        matrix2 = np.resize(matrix, (n, m))
                        for i in range(0, n):
                            for j in range(0, m):
                                matrix2[i, j] = matrix[i, j]
                        # randomly generate list of job weights (uniform distribution)
                        jobweights = np.random.randint(1, 101, n)
                        #print('After: ' + str(matrix2))
                        for alpha_dev in [0.1, 0.2, 0.3, 0.4, 0.5]:
                            # output file
                            filename = file[file.rfind(os.sep) + 1:file.rfind(os.sep) + 6]
                            filename += '%d' % (alpha_dev*100)
                            filename += file[file.rfind(os.sep) + 8:file.rfind(os.sep) + 10]
                            filename += '_' + str(n) + '_' + str(m) + '_' + '%d' % (alpha_dev*100)
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
                                    np.savetxt(output_file, jobweights, fmt='%d')
                                mlist = []
                                for x in range(1, m):
                                    mlist.append(' m{} |'.format(x))
                                mlist.append(' m{}\r\n'.format(m))
                                output_file.write('# P_bar :' + ''.join(mlist))
                                # Ying instance files always have integer processing times
                                np.savetxt(output_file, matrix2, fmt='%d')
                                # Generate processing time deviations, based on alpha value
                                matrix_dev = alpha_dev * matrix2
                                # processing time deviations may be float
                                output_file.write('# P_hat :' + ''.join(mlist))
                                np.savetxt(output_file, matrix_dev, fmt='%.2f')
                            finally:
                                output_file.close()
                            print("\nCreated output instance file {0}".format(filename))
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
