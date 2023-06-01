#%%
from config import *

import pandas as pd
#%% Define list of libraries
libs = [f"R{i}" for i in range(5,9)]

#%% Load vireo assignments and filtered barcodes
vs = dict()
vireodir = Path("/ssd1/rng/data/incov_bcell_multiome/")
for lib in libs:
    df = pd.read_csv(vireodir / f"{lib}/vireo/donor_ids.tsv", sep='\t')
    df = df.set_index('cell')
    inner = pd.read_csv(
        f"../output/cell_calling/accepted_barcodes_inner_{lib}.csv", 
        header=None)
    # Annotate barcodes based on filtering method
    df['inner'] = 0
    df.loc[inner[0], 'inner'] = 1
    df['library'] = lib
    vs[lib] = df
df = pd.concat(vs.values())

# %% Load expected cell ratio
sample_sheet = pd.read_excel("../data/sample_sheet.xlsx")
sample_dict = dict()
for lib in libs:
    inds = sample_sheet['library'] == lib
    temp = dict()
    temp[sample_sheet.loc[inds, 'donor1'].values[0]] = 'donor1'
    temp[sample_sheet.loc[inds, 'donor2'].values[0]] = 'donor2'
    sample_dict[lib] = temp

# %%
for lib in libs:
    inds = df['library'] == lib
    df.loc[inds, 'donor_num'] = df.loc[inds, 'donor_id'].map(sample_dict[lib])
#%%
# %% Outer ratio
counts = df.groupby(['library'])['donor_num'].value_counts().unstack()
ratios = counts.div(counts.min(axis=1), axis=0)
ratios
# %% Inner ratio
counts = df.loc[df['inner']==1,].groupby(['library'])['donor_num'].value_counts().unstack()
ratios = counts.div(counts.min(axis=1), axis=0)
ratios
#%%
# %% Outer - Inner ratio
counts = df.loc[df['inner']==0,].groupby(['library'])['donor_num'].value_counts().unstack()
ratios = counts.div(counts.min(axis=1), axis=0)
ratios
# %%
counts = df.groupby(['library'])['donor_id'].value_counts().to_frame()
counts
# %%
counts = df.groupby(['library', 'inner'])['donor_id'].value_counts().to_frame()
counts
# %%
freq = df.groupby(['library', 'inner'])['donor_id'].value_counts(normalize=True).to_frame()

# %%
outer = df.groupby(['library'])['donor_id'].value_counts(
    normalize=True).to_frame()
inner = df.loc[df['inner']==1,].groupby(['library'])['donor_id'].value_counts(
    normalize=True).to_frame()
diff = df.loc[df['inner']==0,].groupby(['library'])['donor_id'].value_counts(
    normalize=True).to_frame()
merged = pd.concat([outer, inner, diff], axis=1, ignore_index=False)
merged.columns = ['outer', 'inner', 'diff']
merged
# %%
