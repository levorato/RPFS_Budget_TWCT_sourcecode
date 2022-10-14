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
import random

class Range(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __eq__(self, other):
        return self.start <= other <= self.end

    def __str__(self):
        return "range(" + str(self.start) + ", " + str(self.end) + ")"


def main(argv):
    csv.field_size_limit(1000000000)

    parser = argparse.ArgumentParser(description='Generate instance file (.in) for the flowshop problem.')
    parser.add_argument('-n', type=int, required=True, help="number of machines (pocos)")
    parser.add_argument('-m', type=int, required=True, help="number of jobs (operacoes)")
    parser.add_argument('-c', type=float, required=True, choices=[Range(0.0, 1.0)],
                        help="percentage of critical operations")
    args = parser.parse_args()

    print 'Number of machines is ', str(args.n)
    print 'Number of jobs is ', str(args.m)
    print 'Percentage of critical operations is ', str(args.c)

    processInstanceFiles(args.n, args.m, args.c)


def processInstanceFiles(n, m, c):

    # measure elapsed time during instance generation
    start = time.time()

    filename = 'n' + str(n) + 'm' + str(m) + 'c' + str(c) + ".in"
    current_directory = 'instances'
    output_file = open(os.path.join(current_directory, filename), 'w')

    try:
        # linha 1 => numero de pocos (n)  e numero de operacoes (m).
        output_file.write("{0} {1}\r\n".format(n, m))

        # linha 2 => m valores: 1 se i-esima operacao exige recurso critico e 0 caso contrario.
        C = [0 for i in xrange(m)]
        critical_list = range(m)
        random.shuffle(critical_list)
        num_critical = int(c * m)
        for i in xrange(num_critical):
            C[critical_list[i]] = 1
        print "Number of critical operations is {0}".format(num_critical)
        print "Generated C = {0}".format(C)
        Cstr = str(C)
        Cstr = Cstr.replace('[', '')
        Cstr = Cstr.replace(']', '')
        Cstr = Cstr.replace(',', '')
        output_file.write("{0}\r\n".format(Cstr))

        # linhas 3 a n+2 => vazoes dos pocos (m^3/hora), 1 valor por linha.
        # n linhas, cada uma com a vazao de um poco: um numero aleatorio entre 10.00 e 200.00
        V = [0 for i in xrange(n)]
        for i in xrange(n):
            V[i] = round(100 * round(random.uniform(0.1, 2), 4), 2)
            output_file.write("{0}\r\n".format(V[i]))
        print "Generated V = {0}".format(V)

        # n linhas: n+2 a 2n+1 => m valores correspondentes aos tempos das operacoes.
        for i in xrange(n):
            Pline = []
            for i in xrange(m):
                time_times_2 = round(random.uniform(0, 24), 0)
                # gera tempos de operacoes entre 0 e 12 horas (com intervalos de meia hora, ie, 0.5)
                Pline.append(float(time_times_2 / 2))
            Pstr = str(Pline)
            Pstr = Pstr.replace('[', '')
            Pstr = Pstr.replace(']', '')
            Pstr = Pstr.replace(',', '')
            print str(Pline) + '\n'
            output_file.write("{0}\r\n".format(Pstr))
    finally:
        output_file.close()
    print "\nCreated output dat instance file {0}".format(filename)
    end = time.time()
    elapsed = end - start
    print "Instance generation took {0:.2f} seconds.".format(elapsed)
    # end loop
    print "\nDone.\n"


if __name__ == "__main__":
    main(sys.argv[1:])
