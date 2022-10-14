#!/usr/bin/env python
# coding: utf-8

# to run this script, please install:
# sudo apt-get install python-numpy python-scipy python-matplotlib ipython ipython-notebook python-pandas
# python-sympy python-nose python-tk
# OR pip install pandas <etc.>

import sys, getopt
import csv
import os
import os.path
import argparse
import pandas as pd
import scipy.stats as stats
import numpy as np
import matplotlib as mpl
from statsmodels.stats.multicomp import MultiComparison
from matplotlib.backends.backend_pdf import PdfPages
import statsmodels.api as sm
from statsmodels.formula.api import ols
#import researchpy as rp
import fnmatch
import copy

## agg backend is used to create plot as a .png file
mpl.use('agg')

import matplotlib.pyplot as plt


def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def main(argv):

    csv.field_size_limit(100000)
    pd.set_option('display.max_columns', None)
    mpl.rcParams['figure.figsize'] = (16.0, 12.0)
    mpl.style.use('ggplot')

    parser = argparse.ArgumentParser(description='Process PFSP result files.')
    parser.add_argument('--folders', nargs='+',
                        help='the folders containing the result files (one for each experiment)')
    parser.add_argument('--filefilter', default='*-solution_info.csv', required=False,
                        help='the file extension for result files (default: *.csv)')
    args = parser.parse_args()
    folders = args.folders
    filter = args.filefilter
    args = parser.parse_args()

    print('Input folders are ', folders)
    print('File filter is ', filter)
    processResult(folders, filter)


def extract_instance_name(filename):
    instance_name = filename[:filename.find('.')]
    if filename.find('reversed') >= 0:
        instance_name += '_rev'
    return instance_name


def processResult(folders, filter):
    for folder in folders:
        print("Processing folder " + ''.join(folder))
        output_path = os.path.join(folder, "summary")
        if not os.path.exists(output_path):
            os.makedirs(output_path)

        t1_list = []
        t2_list = []
        count = 1
        for root, subFolders, files in os.walk(folder):
            # sort dirs and files
            files.sort()
            for filename in fnmatch.filter(files, filter):
                column_list, value_list = processTraceFile(folder, filename)
                phatom_filename = filename[:filename.rfind("-")] + "-phatomed_count.csv"
                column_list2, value_list2 = processPhatomedInfoFile(folder, phatom_filename)
                t1_list.append([extract_instance_name(filename)] + list(value_list))
                t2_list.append([extract_instance_name(filename)] + list(value_list2))
                count += 1
                #if count > 40:
                #    break
            # end for
        # end for
        df = pd.DataFrame(t1_list, columns = ['instance'] + list(column_list))
        df2 = pd.DataFrame(t2_list, columns=['instance'] + list(column_list2))
        df = df.join(df2.set_index('instance'), on='instance')
        filename = os.path.join(output_path, "PFSP_Summary_All_Instances.csv")
        df.to_csv(filename)
        # end process all result files of an instance / filename


def print_df_head(df, name):
    print('\n=========================================\n' + name + '\n=========================================\n')
    print("Grouped dataframe deterministic:\n" + str(df[df['RTCS_Type'] == "Deterministic"].head()))
    print("Grouped dataframe robust:\n" + str(df[df['RTCS_Type'] == "Robust"].head()))

def print_full_df(df, name):
    print('\n=========================================\n' + name + '\n=========================================\n')
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        print(df)


def processTraceFile(folder, filename):
    print("\nProcessing trace file : " + filename + "...\n")
    df = pd.read_csv(os.path.join(folder, filename), sep = ':', header = None, names = ['key', 'value'])
    column_list = df['key'].values
    value_list = df['value'].values
    return column_list, value_list


def processPhatomedInfoFile(folder, filename):
    print("Processing phatom info file : " + filename + "...\n")
    df = pd.read_csv(os.path.join(folder, filename), sep=',')  #, index_col=["LB_num"])
    df = df.rename(columns=lambda x: x.strip())  # strip whitespace from column names
    column_list = list(df["LB_num"])
    df = df[df["LB_num"] != 0]
    df = df.astype(int)
    #print(df)
    nn_values = list(df.NN.values)
    inv_values = list(df['invocations'])
    column_list = [x for x in column_list if x != 0]
    column_list = [("Inv_LB_" + str(x)) for x in column_list] + [("NN_LB_" + str(x)) for x in column_list]
    #print("Columns: " + str(column_list))
    value_list = inv_values + nn_values
    #print("Values: " + str(value_list))
    return column_list, value_list


if __name__ == "__main__":
    main(sys.argv[1:])
