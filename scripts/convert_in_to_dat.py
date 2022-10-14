# Script to convert generic instance file (.in) to CPLEX dat file (.dat)

import sys, getopt
import csv
import StringIO
import glob
import os
import os.path
import argparse
import time
import math

def main(argv):
    csv.field_size_limit(1000000000)

    parser = argparse.ArgumentParser(description='Convert instance file (.in) to CPLEX dat file.')
    parser.add_argument('folders', nargs='+',
                        help='the folders containing the instance files')
    parser.add_argument('--filefilter', default='.in', required=False,
                        help='the filename for instance files (default: .in)')
    args = parser.parse_args()
    folders = args.folders
    filter = args.filefilter

    args = parser.parse_args()

    print 'Input folders are ', folders
    print 'File filter is ', filter

    processInstanceFiles(folders, filter, True)
    processInstanceFiles(folders, filter, False)


def processInstanceFiles(folders, filter, old_format):
    for folder in folders:
        print "Processing folder " + ''.join(folder)

        for root, subFolders, files in os.walk(folder):
            # sort dirs and files
            subFolders.sort()
            files.sort()

            if len(files):
                file_list = [f for f in files if filter in f]

                for file in file_list:
                    directory = root + os.sep + 'dat'
                    if not os.path.exists(directory):
                        os.makedirs(directory)

                    # measure elapsed time during instance generation
                    start = time.time()

                    file = os.path.join(root, file)
                    filename = file
                    print "Processing file " + filename
                    filename = filename[filename.rfind(os.sep)+1:filename.rfind('.')]
                    content_file = open(file, 'r')
                    if old_format:
                        filename += '-v1'
                    else:
                        filename += '-v2'
                    filename += ".dat"
                    output_file = open(os.path.join(directory, filename), 'w')

                    try:
                        # detect the dialect (separator used in the csv file
                        #dialect = csv.Sniffer().sniff(content_file.read(1024), delimiters=" ")

                        # reads the contents of csv to pandas dataframe
                        content_file2 = open(file, 'r')
                        text_content = content_file2.read()
                        output = StringIO.StringIO(text_content)
                        reader = csv.reader(output, delimiter =' ')

                        # reads the lines of the instance file
                        row_count = 1
                        V = []
                        C = []
                        p = []
                        s = []
                        upper_bound = 0
                        num_operacoes = 0
                        for row in reader:
                            linestring = ''.join(row)
                            column = []
                            for col in row:
                                column.append(col)
                            if len(row) == 0:
                                continue
                            print "Reading row {0}".format(column)

                            if row_count == 1:  # linha 1 => numero de pocos (n)  e numero de operacoes (m).
                                n = long(column[0])
                                m = long(column[1])
                            elif row_count == 2: # linha 2 => m valores: 1 se i-esima operacao exige recurso critico e 0 caso contrario.
                                for i in xrange(m):
                                    C.append(int(column[i]))
                                print "Read C = {0}".format(C)
                            elif row_count <= n + 2: # linhas 3 a n+2 => vazoes dos pocos (m^3/hora), 1 valor por linha.
                                # converte os valores de vazao de 'por hora' para 'por meia hora'
                                V.append(2 * float(column[0]))
                            else: # linhas n+2 a 2n+1 => m valores correspondentes aos tempos das operacoes.
                                Pline = []
                                Sline = []
                                setup_time_sum = 0
                                num_operacoes = 0
                                for i in xrange(m):
                                    time_value = float(column[i])
                                    # discretiza os intervalos de tempo em numeros de periodos de meia hora
                                    time_in_half_hour_elements = int(2 * time_value)
                                    if C[i] == 1 or old_format:  # se a operacao i eh critica, insere na matrix de tempos D
                                        Pline.append(time_in_half_hour_elements)
                                        Sline.append(setup_time_sum)
                                        setup_time_sum = 0
                                        num_operacoes += 1
                                    else:
                                        setup_time_sum += time_in_half_hour_elements
                                    upper_bound += time_in_half_hour_elements
                                p.append(Pline)
                                s.append(Sline)

                            row_count += 1
                        # end read lines
                        content_file2.close()
                        print "Successfully read input file, generating dat file."

                        # writes the header info
                        output_file.write("numPocos = {0};\r\n".format(n))
                        output_file.write("numOperacoes = {0};\r\n".format(num_operacoes))
                        # the following line was replaced by setup times matrix
                        if old_format:
                            output_file.write("numPeriodos = {0};\r\n".format(upper_bound))
                            output_file.write("C = {0};\r\n".format(C))
                        output_file.write("V = {0};\r\n".format(V))
                        output_file.write("p = {0};\r\n".format(p))
                        if not old_format:
                            output_file.write("s = {0};\r\n".format(s))
                    finally:
                        content_file.close()
                        output_file.close()
                    print "\nCreated output dat instance file {0}".format(filename)
                    end = time.time()
                    elapsed = end - start
                    print "Instance conversion took {0:.2f} seconds.".format(elapsed)
                # end loop
                # process last file


        print "\nDone.\n"


if __name__ == "__main__":
    main(sys.argv[1:])
