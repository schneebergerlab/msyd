import pandas as pd
"""
The output of msyd call --realign should return all haplotypes blocks in the population.
By comparing non-core-syn hapblocks, we can identify structurally rearranged blocks.

Traditionally, in a pan-genome-graph, annotating sequence segments as (for example) a genomic insertion is not very meaningful. This is because a variation can only be defined within the scope of a reference. As pan-genomes do not have a defined reference sequence, calling variation does not work.

Here, we take a step back, instead of trying to define the sequence space of a population, we try to define the haplotype/variation space of a population within the context of a user-defined reference sequence. The idea is, if we can define variation in a population then we provide a more meaningful and intuitive way to understand the genomic properties of a population.

"""
from collections import deque
import pandas as pd
def rng_to_tuple(r):
    """
    required `r` format Chr:Start-End
    """
    c, r = r.split(':')
    s, e = tuple(map(int, r.split('-')))
    return(c, s, e)
# END

def read_pff2(pff):
    """
    Reads the newer variant of PFF file. A reader for the older variant is available here: msyd.pyxfiles.io.read_pff
    """
    out = deque()
    with open(pff) as fin:
        line = fin.readline().strip().split()
        samples = line[4:]
        for line in fin:
            # if line == '': continue
            # if line is None: continue
            line = line.strip().split()
            if line == []: continue
            try:
                anno = line[3]
            except IndexError:
                print(line)
            rreg = (line[0], int(line[1]), int(line[2]))
            qryreg = {samples[i]: rng_to_tuple(q.split(',')[0]) for i, q in enumerate(line[4:])  if q != '-'}
            refsample = set([q.split(',')[1] for q in line[4:] if q != '-'])
            # Test that all ranges have same reference
            try:
                assert len(refsample) == 1
            except AssertionError:
                logger.error(f'Single region has multiple references {line}')
                return
            out.append((anno, rreg, qryreg, list(refsample)[0]))
    return out
# END

def find_invhaps(f):
    """
    Read realigned PFF and call inversions. Later, this function would be called directly within the main.py and the multisyn blocks would be parsed as arguments.
    """
    from msyd.util import CustomFormatter
    from msyd.io import read_pff, readsyriout

    logger = CustomFormatter.getlogger('find_invhaps')

    def concatsyriout(syrifins, qrynames):
        """
        Given a list of syri output files, read them, select syn and SR regions, and return a concatenated dataframe
        """
        try:
            assert len(syrifins) == len(qrynames)
        except AssertionError:
            logger.error(f'Number of syri output files {len(syrifins)} does not match number of qrynames {len(qrynames)}')
        outdf = pd.DataFrame()
        for syrifin, qryname in zip(syrifins, qrynames):
            syridf, mapid = readsyriout(syrifin)
            syridf['qgen'] = qryname
            outdf = pd.concat([outdf, syridf])
        return outdf
    # END

    pfffin = '/home/ra98jam/projects/tests/BKP_colcc_10_thal.pff'
    synblocks = read_pff2(pfffin)

    # Inversion haplotypes would be called as inversion by Syri. Read the syri output and then match hap-blocks that are
    # inverted to each other.
    # Read syri annotations
    # TMP CODE
    indir = '/home/ra98jam/projects/msyd/results/thaliana/'
    inpdf = pd.read_table(f'{indir}genomes.csv', header=None)
    syrifins = [f'{indir}{s}' for s in inpdf.loc[:, 2].values.tolist()]
    qrynames = inpdf.loc[:, 0].values.tolist()
    # END TMP CODE
    # TODO : define syrifins and qrynames
    allsyri = concatsyriout(syrifins, qrynames)
    allsyri.sort_values(by='achr astart aend'.split(), inplace=True)
    # Remove "SYN" annotations
    allsyri = allsyri.loc[~(allsyri['type'] == 'SYN')]
    for row in synblocks:
        print(row)
        break




    return
# END
CP116280.1,784681,893258,OX291513.1,791184,899980,INV,IP-Evs-12
