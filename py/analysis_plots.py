# <editor-fold desc="Define Imports">
import numpy as np
from matplotlib import pyplot as plt
from msyd.coords import Pansyn
from collections import deque
from hometools.hometools import view
from hometools.plot import cleanax
# </editor-fold>

# <editor-fold desc="Define constants">

pfffin = 'chr02_all_hap.pff'
chrsize = 46102915
outdir = 'C:\\Users\\ra98jam\\Dropbox\\projects\\msyd\\figures\\potato\\'

# </editor-fold>

# <editor-fold desc="Define functions and classes">
# </editor-fold>

# Read PFF (from a single chr) analysis and save the coordinates in a deque
rec = deque()
with open(pfffin, 'r') as fin:
    snames = fin.readline().strip().split()
    for line in fin:
        line = line.strip().split()[1:]
        rec.append([None if c == '.' else list(map(int, c.split(',')[0].rsplit(':', maxsplit=1)[1].split('-'))) for c in line])
snames.pop(0)
snames[0] = 'DM'

total_ref = [i[0][1] - i[0][0] for i in rec]
fig, ax = plt.subplots(figsize=(4, 4))
ax = loghist(total_ref, ax=ax, bins=100)
ax = cleanax(ax)
ax.set_xlabel('Length of region on reference genome')
ax.set_ylabel('Number of regions')
plt.tight_layout()
plt.savefig(f'{outdir}ref_DM_coresyn_length.pdf')

# generate syntenic region length distribution and position barplots
maxx = max([r[1] - r[0] for row in rec for r in row if r is not None])
logbins = np.logspace(np.log10(1), np.log10(maxx), 101)
fig = plt.figure(figsize=(12, 12))
coresyn_index = deque()
for i, sname in enumerate(snames[1:]):
    i = i+1
    print(sname)
    total_ref = [r[i][1] - r[i][0] for r in rec if r[i] is not None]
    ax = fig.add_subplot(6, 6, i)
    ax.hist(total_ref, bins=logbins)
    ax.set_xscale('log')
    ax = cleanax(ax)
    ax.set_xlabel(sname)
    ax.set_ylabel('Count')
    plt.hist(([1, 2, 3, 5, 6], [1, 2, 3, 5, 6]), stacked=True, density=True)


    hist, bins = np.histogram(x, bins=bins)
    logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))
    ax.hist(x, bins=logbins)




plt.tight_layout()
plt.savefig(f'{outdir}ref_DM_coresyn_length.pdf')

    print(sname)


