#!/usr/bin/env python

import argparse
import json
import sys
import re
import os
from collections import defaultdict


def parse_pattern(pattern):
    """
    >>> pattern_dict = parse_pattern("C8L16C8L16C8L1U12T18")
    >>> pattern_dict['C']
    [[0, 8], [24, 32], [48, 56]]
    >>> pattern_dict['L']
    [[8, 24], [32, 48], [56, 57]]
    """
    pattern_dict = defaultdict(list)
    p = re.compile(r'([CLUNT])(\d+)')
    tmp = p.findall(pattern)
    if not tmp:
        sys.exit()(f'Invalid pattern: {pattern}')
    start = 0
    for p, length in tmp:
        end = start + int(length)
        pattern_dict[p].append([start, end])
        start = end
    return pattern_dict

def get_solo_pattern(pattern) -> str:
    """
    Returns:
        starsolo_cb_umi_args
    """
    pattern_dict = parse_pattern(pattern)
    if len(pattern_dict['U']) != 1:
        sys.exit(f'Error: Wrong pattern:{pattern}. \n Solution: fix pattern so that UMI only have 1 position.\n')
    ul, ur = pattern_dict["U"][0]
    umi_len = ur - ul

    if len(pattern_dict["C"]) == 1:
        solo_type = 'CB_UMI_Simple'
        l, r = pattern_dict["C"][0]
        cb_start = l + 1
        cb_len = r - l
        umi_start = ul + 1
        cb_str = f'--soloCBstart {cb_start} --soloCBlen {cb_len} --soloCBmatchWLtype 1MM'
        umi_str = f'--soloUMIstart {umi_start} --soloUMIlen {umi_len} '
    else:
        solo_type = 'CB_UMI_Complex'
        cb_pos = ' '.join([f'0_{l}_0_{r-1}' for l, r in pattern_dict["C"]])
        umi_pos = f'0_{ul}_0_{ur-1}'
        cb_str = f'--soloCBposition {cb_pos} --soloCBmatchWLtype EditDist_2 '
        umi_str = f'--soloUMIposition {umi_pos} --soloUMIlen {umi_len} '

    starsolo_cb_umi_args = f'--soloType {solo_type} ' + cb_str + umi_str
    return starsolo_cb_umi_args


def write_output(parsed_protocol, starsolo_cb_umi_args, whitelist, sample):
    """
    write files:
        path "*starsolo_cb_umi_args.txt", emit: starsolo_cb_umi_args
        path "*parsed_protocol.txt", emit: parsed_protocol
        path "*whitelist.txt", emit: whitelist
    """
    with open(sample + '.parsed_protocol.txt', 'w') as f:
        f.write(parsed_protocol)
    with open(sample + '.starsolo_cb_umi_args.txt', 'w') as f:
        f.write(starsolo_cb_umi_args)
    with open(sample + '.whitelist.txt', 'w') as f:
        f.write(whitelist)


if __name__ == '__main__':
    """
    args:
        --sample ${meta.id} \\
        --R1_fastq ${reads[0]} \\
        --json_file ${json_file} \\
        --protocol ${protocol} \\
        --pattern ${pattern} \\
    output:
        path "*starsolo_cb_umi_args.txt", emit: starsolo_cb_umi_args
        path "*parsed_protocol.txt", emit: parsed_protocol
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=True)
    parser.add_argument('--R1_fastq', required=True)
    parser.add_argument('--json_file', required=True)
    parser.add_argument('--protocol', required=True)
    parser.add_argument('--pattern')
    # add version
    parser.add_argument('--version', action='version', version='1.0')

    args = parser.parse_args()

    if args.protocol == 'new':
        # TODO
        pass

    elif args.protocol == 'auto':
        # TODO
        pass

    else:
        # get the protocol from protocols.json
        with open(args.json_file) as f:
            p = json.load(f)[args.protocol]
        parsed_protocol = args.protocol
        pattern = p['pattern']
        whitelist = p['whitelist']

    starsolo_cb_umi_args = get_solo_pattern(pattern)
    write_output(parsed_protocol, starsolo_cb_umi_args, whitelist, args.sample)







