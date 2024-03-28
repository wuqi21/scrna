#!/usr/bin/env python

import argparse
import itertools
import json
import math
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict

import pyfastx

SAM_attributes = 'NH HI nM AS CR UR CB UB GX GN '

def get_seq_str(seq, sub_pattern):
    """
    join seq slices.

    Args:
        seq: usually R1 read
        sub_pattern: [slice(0,8),slice(16,24)]

    Returns:
        joined intervals seq

    Raises:
        IndexError: if seq length is not enough.

    >>> sub_pattern_dict = [slice(0,2)]
    >>> seq = "A" * 2 + "T" * 2
    >>> get_seq_str(seq, sub_pattern_dict)
    'AA'
    """
    seq_len = len(seq)
    expect_len = sub_pattern[-1].stop
    if seq_len < expect_len:
        raise IndexError(f"read length({seq_len} bp) less than expected length({expect_len} bp) in read: {seq}")
    return "".join([seq[x] for x in sub_pattern])


def findall_mismatch(seq, n_mismatch=1, bases='ACGTN'):
    """
    choose locations where there's going to be a mismatch using combinations
    and then construct all satisfying lists using product

    Return:
    all mismatch <= n_mismatch set.

    >>> answer = set(["TCG", "AAG", "ACC", "ATG", "ACT", "ACN", "GCG", "ANG", "ACA", "ACG", "CCG", "AGG", "NCG"])
    >>> seq_set = findall_mismatch("ACG")
    >>> seq_set == answer
    True
    """
    seq_set = set()
    seq_len = len(seq)
    if n_mismatch > seq_len:
        n_mismatch = seq_len
    for locs in itertools.combinations(range(seq_len), n_mismatch):
        seq_locs = [[base] for base in seq]
        for loc in locs:
            seq_locs[loc] = list(bases)
        for poss in itertools.product(*seq_locs):
            seq_set.add(''.join(poss))
    return seq_set

def get_mismatch_dict(seq_list, n_mismatch=1):
    """
    Return:
    mismatch dict. Key: mismatch seq, value: seq in seq_list

    >>> seq_list = ["AACGTGAT", "AAACATCG"]
    >>> mismatch_dict = get_mismatch_dict(seq_list)
    >>> mismatch_dict["AACGTGAA"] == "AACGTGAT"
    True
    """
    mismatch_dict = {}
    for seq in seq_list:
        seq = seq.strip()
        if seq == '':
            continue
        for mismatch_seq in findall_mismatch(seq, n_mismatch):
            mismatch_dict[mismatch_seq] = seq
    return mismatch_dict


def read_one_col(fn):
    """read one column file into list"""
    with open(fn) as f:
        return [x.strip() for x in f]


def get_protocol_dict(assets_dir):
    json_file = os.path.join(assets_dir, "protocols.json")
    protocol_dict = json.load(open(json_file))
    whitelist_dir = os.path.join(assets_dir, "whitelist")
    # add folder prefix
    for protocol in protocol_dict:
        cur = protocol_dict[protocol]
        bc = cur.get("bc", [])
        linker = cur.get("linker", [])
        if bc:
            cur["bc"] = [os.path.join(whitelist_dir, protocol, x) for x in bc]
        if linker:
            cur["linker"] = [os.path.join(whitelist_dir, protocol, x) for x in linker]
        cur["pattern_dict"] = parse_pattern(cur["pattern"])
    return protocol_dict


class Protocol:
    """
    Auto detect singleron protocols from R1-read
    GEXSCOPE-MicroBead
    GEXSCOPE-V1
    GEXSCOPE-V2
    sCircle-V1
    """

    def __init__(self, fq1_list, sample, assets_dir='assets/', max_read=10000):
        '''
        Args:
            assets_dir: Expects file 'protocols.json' and 'whitelist/{protocol}' folder under assets_dir
        '''
        self.fq1_list = fq1_list
        self.max_read = max_read
        self.sample = sample
        self.protocol_dict = get_protocol_dict(assets_dir)
        self.scircle_R1_LEN = self.protocol_dict["sCircle-V1"]["pattern_dict"]["U"][-1].stop

    def check_protocol(self):
        """check protocol in the fq1_list"""
        fq_protocol = {}
        for fastq1 in self.fq1_list:
            protocol = self.get_protocol(fastq1)
            fq_protocol[fastq1] = protocol
        if len(set(fq_protocol.vals())) != 1:
            sys.exit(f'Error: multiple protocols are not allowed for one sample: {self.sample}! \n' + str(fq_protocol))
        protocol = list(fq_protocol.vals())[0]
        return protocol

    def check_linker(self, seq, protocol, i):
        """check if seq matches the linker i of protocol"""
        linker_fn = self.protocol_dict[protocol]["linker"]
        linker = read_one_col(linker_fn[i-1])
        linker_mismatch_dict = get_mismatch_dict(linker)

        pattern = self.protocol_dict[protocol]["pattern"]
        pattern_dict = parse_pattern(pattern)
        seq_linker = seq[pattern_dict["L"][i-1]]

        return seq_linker in linker_mismatch_dict


    def seq_protocol(self, seq):
        """
        Returns: protocol or None

        >>> runner = Protocol("fake_fq1_string", "fake_sample")
        >>> seq = "TCGACTGTC" + "ATCCACGTGCTTGAGA" + "TTCTAGGAT" + "TCAGCATGCGGCTACG" + "TGCACGAGA" + "C" + "CATATCAATGGG" + "TTTTTTTTTT"
        >>> runner.seq_protocol(seq)
        'GEXSCOPE-V2'
        >>> seq = "NCAGATTC" + "TCGGTGACAGCCATAT" + "GTACGCAA" + "CGTAGTCAGAAGCTGA" + "CTGAGCCA" + "C" + "TCCGAAGCCCAT" + "TTTTTTTTTT"
        >>> runner.seq_protocol(seq)
        'GEXSCOPE-V1'
        >>> seq = "NCAGATTC" + "TCGGTGACAGCCATAT" + "GTACGCAA" + "CGTAGTCAGAAGCTGA" + "CTGAGCCA"  + "TCCGAAGCC" + "CTGTCT"
        >>> runner.seq_protocol(seq)
        'GEXSCOPE-sCircle'
        >>> seq = "NCAGATTC" + "TCGGTGACAGCCATAT" + "GTACGCAA" + "CGTAGTCAGAAGCTGA" + "CTGAGCCA"  + "TCCGAAGCC"
        >>> runner.seq_protocol(seq)
        'GEXSCOPE-sCircle'
        >>> seq = "ATCGATCGATCG" + "ATCGATCG" + "C" + "TTTTTTTTTT"
        >>> runner.seq_protocol(seq)
        'GEXSCOPE-MicroBead'
        """
        # check "GEXSCOPE-sCircle" first as it is similar to "GEXSCOPE-V1"
        protocol = "GEXSCOPE-sCircle"
        if self.check_linker(seq, protocol, 1) and (len(seq) == self.scircle_R1_LEN or self.check_linker(seq, protocol, 3)):
            return protocol

        for protocol in ["GEXSCOPE-V1", "GEXSCOPE-V2"]:
            if self.check_linker(seq, protocol, 1):
                return protocol

        # check if it is MicroBead
        if seq[16:20] != "TTTT" and seq[22:26] == "TTTT":
            return "GEXSCOPE-MicroBead"


    def get_protocol(self, fq1):
        results = defaultdict(int)

        with pyfastx.Fastx(fq1) as fh:
            for _ in range(self.max_read):
                entry = fh.__next__()
                seq = entry.sequence
                protocol = self.seq_protocol(seq)
                if protocol:
                    results[protocol] += 1
        # if it is 0, then no other linker types
        if results["scopeV1"] != 0:
            results["scopeV1"] += results["scopeV2.1.1"]
        sorted_counts = sorted(results.items(), key=lambda x: x[1], reverse=True)
        self.get_protocol.logger.info(sorted_counts)

        protocol, read_counts = sorted_counts[0][0], sorted_counts[0][1]
        percent = float(read_counts) / self.max_read
        if percent < 0.5:
            self.get_protocol.logger.warning("Valid protocol read counts percent < 0.5")
        if percent < 0.1:
            self.get_protocol.logger.error("Valid protocol read counts percent < 0.1")
            raise Exception(
                'Auto protocol detection failed! '
            )
        protocol.get_protocol.logger.info(f'protocol: {protocol}')

        return protocol



def parse_pattern(pattern):
    """
    >>> pattern_dict = parse_pattern("C8L16C8L16C8L1U12T18")
    >>> pattern_dict['C']
    [slice(0, 8, None), slice(24, 32, None), slice(48, 56, None)]
    >>> pattern_dict['L']
    [slice(8, 24, None), slice(32, 48, None), slice(56, 57, None)]
    """
    pattern_dict = defaultdict(list)
    p = re.compile(r'([CLUNT])(\d+)')
    tmp = p.findall(pattern)
    if not tmp:
        sys.exit(f'Invalid pattern: {pattern}')
    start = 0
    for p, length in tmp:
        end = start + int(length)
        pattern_dict[p].append(slice(start,end))
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
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', required=True)
    parser.add_argument('--R1_fastq', required=True)
    parser.add_argument('--R2_fastq', required=True)
    parser.add_argument('--asset_dir', required=True)
    parser.add_argument('--protocol', required=True)
    parser.add_argument('--pattern')
    parser.add_argument('--whitelist')
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







