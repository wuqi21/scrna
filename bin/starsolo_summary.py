#!/usr/bin/env python

import pandas as pd
import argparse
import utils
import json
import csv
from collections import defaultdict

def parse_summary(f):
    parsed_data = {}
    reader = csv.reader(f)
    for row in reader:
        parsed_data[row[0]] = row[1]
    return parsed_data


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Starsolo summary')
    parser.add_argument('--read_stats', help='cellReadsStats file')
    parser.add_argument('--barcodes', help='barcode file')
    parser.add_argument('--summary', help='summary file')
    parser.add_argument('--sample', help='sample name')
    args = parser.parse_args()
    dtypes = defaultdict(lambda :'Int32')
    dtypes['CB'] = 'object'
    df = pd.read_csv(args.read_stats, sep='\t', header=0, index_col=0, skiprows=[1], dtype=dtypes) # skip first line cb not pass whitelist
    df_count = df.loc[:,['nUMIunique','countedU']] # keep dataframe format
    df = df.loc[:,['cbMatch', 'cbPerfect','genomeU', 'genomeM', 'exonic', 'intronic','exonicAS','intronicAS','countedU']]
    s = df.sum()
    # json does not recognize NumPy data types. TypeError: Object of type int64 is not JSON serializable
    valid = int(s['cbMatch'])
    perfect = int(s['cbPerfect'])
    corrected = valid - perfect
    genomeU= int(s['genomeU'])
    genomeM = int(s['genomeM'])
    mapped = genomeU + genomeM
    exonic = int(s['exonic'])
    intronic = int(s['intronic'])
    antisense = int(s['exonicAS'] + s['intronicAS'])
    intergenic = mapped - exonic - intronic - antisense
    countedU = int(s['countedU'])
    data_dict = {
        'sample': args.sample,
        'valid': valid,
        'perfect': perfect,
        'corrected': corrected,
        'genomeU': genomeU,
        'genomeM': genomeM,
        'mapped': mapped,
        'exonic': exonic,
        'intronic': intronic,
        'antisense': antisense,
        'intergenic': intergenic,
        'countedU': countedU
    }
    read_stats_file = args.sample + '.read_stats.json'
    with open(read_stats_file, 'w') as f:
        json.dump(data_dict, f)
    
    # summary
    parsed_data = parse_summary(open(args.summary))
    summary_file = args.sample + '.summary.json'
    with open(summary_file, 'w') as f:
        json.dump(parsed_data, f)

    # UMI count 
    df_count.sort_values(by='nUMIunique', ascending=False, inplace=True)
    cbs = utils.read_one_col(args.barcodes)
    df_count['mark'] = 'BG'
    for cb in cbs:
        df_count.loc[cb, 'mark'] = 'CB'
    UMI_file = args.sample + '.UMI_count.tsv' 
    df_count.to_csv(UMI_file, index=True)