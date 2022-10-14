# Script to convert Ying's robust PFSP instance file (.txt) to Taillard-format instance file (.txt)

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

    parser = argparse.ArgumentParser(description='Convert Ying instance file to Taillard-format file.')
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


def processInstanceFiles(folders, filter, outputpath):
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
                    # measure elapsed time during instance generation
                    start = time.time()
                    file = os.path.join(root, file)
                    filename = file
                    print("Processing file " + filename)
                    matrix = np.loadtxt(file)
                    n = int(np.size(matrix, 0))
                    m = int(np.size(matrix, 1) / 2)
                    print('Instance size: n = {}, m = {}'.format(n, m))
                    #print('Before: ' + str(matrix))
                    # resize matrix and remove the uncertain part
                    matrix2 = np.resize(matrix, (n, m))
                    for i in range(0, n):
                        for j in range(0, m):
                            matrix2[i, j] = matrix[i, j]
                    #print('After: ' + str(matrix2))
                    # output file
                    filename = filename[filename.rfind(os.sep) + 1:filename.rfind('.')]
                    filename += '_' + str(n) + '_' + str(m) + '_inputs.txt'
                    output_file = open(os.path.join(outputpath, filename), 'w')
                    try:
                        output_file.write("# nJobs | nMachines\r\n")
                        output_file.write("{0} {1}\r\n".format(n, m))
                        mlist = []
                        for x in range(1, m):
                            mlist.append(' m{} |'.format(x))
                        mlist.append(' m{}\r\n'.format(m))
                        output_file.write('#' + ''.join(mlist))
                        # Ying instance files always have integer processing times
                        np.savetxt(output_file, matrix2, fmt='%d')
                    finally:
                        output_file.close()
                    print("\nCreated output instance file {0}".format(filename))
                    end = time.time()
                    elapsed = end - start
                    print("Instance conversion took {0:.2f} seconds.".format(elapsed))
                    #break
                # end loop
                # process last file
        print("\nDone.\n")


if __name__ == "__main__":
    main(sys.argv[1:])
