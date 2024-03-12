import pandas as pd
"""
The output of msyd call --realign should return all haplotypes blocks in the population.
By comparing non-core-syn hapblocks, we can identify structurally rearranged blocks.

Traditionally, in a pan-genome-graph, annotating sequence segments as (for example) a genomic insertion is not very meaningful. This is because a variation can only be defined within the scope of a reference. As pan-genomes do not have a defined reference sequence, calling variation does not work.

Here, we take a step back, instead of trying to define the sequence space of a population, we try to define the haplotype/variation space of a population within the context of a user-defined reference sequence. The idea is, if we can define variation in a population then we provide a more meaningful and intuitive way to understand the genomic properties of a population.

"""


def find_invhaps(f):
    """
    Read realigned PFF and call inversions. Later, this function would be called directly within the main.py and the pansyn blocks would be parsed as arguments.
    """
    from msyd.scripts.util import CustomFormatter
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
    synblocks = read_pff(pfffin)

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

    for row in synblocks.itertuples(index=False):
        row = row[0]
        print(row)
        break




    return
# END
