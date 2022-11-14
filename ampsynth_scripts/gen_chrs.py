#!/usr/bin/env python3

import argparse as ap
import copy
import gzip
import multiprocessing
import os
import random

## python script to generate a given number of chromosomes containing a specified number of inversions, translocations, SNPs etc.


parser = ap.ArgumentParser(description="A tool to generate synthetic chromosome sequences containing defined SVs")
parser.add_argument("-i", dest='infasta', required=True, type=str, help="A FASTA file used as base genome sequence.")
parser.add_argument("-c", dest="cores", help="Number of cores to use for parallel computation. Defaults to 4.", type=int, default=4)
parser.add_argument("-n", dest='genomes', default=5, type=int, help="How many genomes to simulate. Default 5.")

parser.add_argument("-s", dest='snp_rate', default=1, type=int, help="SNP rate to simulate, in SNPs/kbp. Default 1.")
parser.add_argument("--indel", dest='indel_rate', default=10, type=int, help="Indel rate to simulate, in indels/Mbp. Default 10.")

parser.add_argument("--max", dest='max', default=100000, type=int, help="Max size of SVs to simulate. Default 100 kbp.")
parser.add_argument("--min", dest='min', default=50, type=int, help="Min size of SVs to simulate. Default 50 bp.")

parser.add_argument("--invs", dest='invs', default=20, type=int, help="No of inversions to simulate for each generated chromosome. Default 20.")
parser.add_argument("--dels", dest='dels', default=10, type=int, help="No of deletions to simulate for each generated chromosome. Default 10.")
parser.add_argument("--inss", dest='inss', default=10, type=int, help="No of insertions to simulate for each generated chromosome. Default 10.")
parser.add_argument("--hdrs", dest='hdrs', default=10, type=int, help="No of highly-diverged regions to simulate for each simulated chromosome. Default 10.")
parser.add_argument("--dups", dest='dups', default=5, type=int, help="No of interspersed duplications to simulate  for each generated chromosome. Default 5.")
parser.add_argument("--tands", dest='tands', default=10, type=int, help="No of tandem duplications to simulate  for each generated chromosome. Default 10.")
parser.add_argument("--tposs", dest='tposs', default=10, type=int, help="No of transpositions within chromosomes to simulate for each generated chromosome. Default 10.")
parser.add_argument("--tlocs", dest='tlocs', default=20, type=int, help="No of translocations across chromosomes to simulate for each generated genome. Default 20.")


ALPHABET=['A', 'C', 'G', 'T']
INDEL_MIN=10
INDEL_MAX=50

def read_fasta(f):
    """Helper function to read a fasta file given as a string path into a dictionary of chromosomes.
    """
    ret = {}
    with gzip.open(f, 'rt') as fin:
        key = ''
        strbuf = ''
        for line in fin:
            if line[0] == '>':
                if key:
                    ret[key] = strbuf
                    strbuf = ''
                key = line[1:].strip()
            else:
                strbuf += line.strip()
        ret[key] = strbuf # add last chr

    return ret

args = parser.parse_args()

ref = read_fasta(args.infasta)

def rand_seq(n):
    ## generates a random dna sequence of length n
    return ''.join(random.choices(ALPHABET, k=n))

def genl(min=args.min, max=args.max):
    return random.randrange(min, max)

def sim_chr(ch):
    """
    Receives a reference chromosome sequence, simulates variants in it, returns it containing the variants and the amount of expected synteny to the reference.
    """
    lch = len(ch)
    syn = lch # expected amount of synteny to reference on this chromosome.
    # calculated assuming infinite sites, i.e. no overlapping variants

    # generate snps, do not save as they should not affect synteny
    ch = list(ch) # temporarily store as a list, as strs are static in python
    for _ in range(int(lch*args.snp_rate/1000)):
        pos = random.randrange(lch)
        ch[pos] = random.choice(ALPHABET) # randomly assign, may not always produce snps
    ch = ''.join(ch)

    # generate indels, do not save as they should not affect synteny?
    for _ in range(int(lch*args.indel_rate/1000000)):
        l = genl(min=INDEL_MIN, max=INDEL_MAX)
        #syn -= l
        pos = random.randrange(lch-l) # ensure there is sufficient space
        if random.choice([True, False]): # randomly choose whether to insert or delete
            ch = ch[:pos] + rand_seq(l) + ch[pos:]
        else:
            ch = ch[:pos] + ch[pos+l:]

    for _ in range(args.invs):
        l = genl()
        syn -= l
        pos = random.randrange(lch-l)
        ch = ch[:pos] + ch[pos+l-1:pos-1:-1] + ch[pos+l:] # first is always inclusive

    for _ in range(args.dels):
        l = genl()
        syn -= l
        pos = random.randrange(lch-l)
        ch = ch[:pos] + ch[pos+l:]

    for _ in range(args.inss):
        l = genl()
        pos = random.randrange(lch)
        ch = ch[:pos] + rand_seq(l) + ch[pos:]

    for _ in range(args.hdrs):
        l = genl()
        syn -= l
        pos = random.randrange(lch-l)
        ch = ch[:pos] + rand_seq(l) + ch[pos+l:]

    for _ in range(args.dups):
        l = genl()
        pos_from = random.randrange(lch-l)
        pos_to = random.randrange(lch-l)
        ch = ch[:pos_to] + ch[pos_from:pos_from+l] + ch[pos_from+l:]

    for _ in range(args.tands):
        l = genl()
        pos = random.randrange(lch-l)
        ch = ch[:pos+l] + ch[pos:pos+l] + ch[pos+l:]

    for _ in range(args.tposs):
        l = genl()
        syn -= 2*l
        pos_from = random.randrange(lch-l)
        pos_to = random.randrange(lch-2*l) # insertion happens after deletion
        seq = ch[pos_from:pos_from+l]
        ch = ch[:pos_from] + ch[pos_from+l:]
        ch = ch[:pos_to] + seq + ch[pos_to:]

    return ch, syn



def sim_genome(n):
    ret = copy.deepcopy(ref)
    chs = list(ref.keys())
    syns_ch = {}
    for id, ch in ref.items():
        seq, syn = sim_chr(ch)
        ret[id] = seq
        syns_ch[id] = syn

    # do translocations
    for _ in range(args.tlocs):
        ch_fr, ch_to = random.sample(chs, k=2)

        l = genl()
        syns_ch[ch_fr] -= l
        ch_fr = ret[ch_fr]
        ch_to = ret[ch_to]

        pos_fr = random.randrange(len(ch_fr) - l)
        pos_to = random.randrange(len(ch_to))

        ch_to = ch_to[:pos_to] + ch_fr[pos_fr:pos_fr+l] + ch_to[pos_to:]
        ch_fr = ch_fr[:pos_fr] + ch_fr[pos_fr+l:]

    # save to files
    with open(str(n)+'.fna', 'wt') as f:
        for ch in chs:
            # write chr header
            f.write(f">{ch} {syns_ch[ch]}\n")
            # write seq with width 80
            seq = ret[ch]
            cov = 0
            while cov < len(seq)-80:
                f.write(seq[cov:cov+80] + '\n')
                cov += 80
            f.write(seq[cov:] + '\n')
            

#for x in range(args.genomes):
#    sim_genome(x)
pool = multiprocessing.Pool(args.cores)
genomes = pool.map(sim_genome, range(args.genomes))
pool.close()
