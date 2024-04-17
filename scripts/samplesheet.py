#!/usr/bin/env python
"""
samplesheet.py folder1,folder2
read all paired fastq files in folders and create a samplesheet.csv file.
"""
import csv
import os
import argparse

def get_read(folder):
    read1_list = [f'_1', f'R1', f'R1_001']
    read2_list = [f'_2', f'R2', f'R2_001']
    fq_list = ['fq', 'fastq']
    suffix1_list = [
        f'{x}.{y}.gz'
        for x in read1_list
        for y in fq_list
    ]
    suffix2_list = [
        f'{x}.{y}.gz'
        for x in read2_list
        for y in fq_list
    ]

    for file in os.listdir(folder):
        r1 = os.path.join(folder, file)
        if os.path.isfile(r1):
            for suffix1, suffix2 in zip(suffix1_list, suffix2_list):
                if file.endswith(suffix1):
                    prefix = file[:-len(suffix1)]
                    r2 = os.path.join(folder, prefix+suffix2)
                    if not os.path.exists(r2):
                        print("warning: {r1} exists, but {r2} not exists. continue")
                        continue
                    prefix = prefix.split("_")[0]
                    print(prefix, r1, r2)
                    yield prefix,r1,r2
                    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('folders')
    args = parser.parse_args()
    with open('samplesheet.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['sample', 'fastq_1', 'fastq_2'])
        for folder in args.folders.split(','):
            folder = os.path.abspath(folder)
            for prefix, r1, r2 in get_read(folder):
                writer.writerow([prefix, r1, r2])
    print("samplesheet.csv created")


