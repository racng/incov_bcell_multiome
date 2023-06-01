#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from config import *

outdir = outputdir / "cell_calling"
outdir.mkdir(exist_ok=True)

#%% Funcitons for joint calling scatter plots
def plot_cellranger_jointcall(bcmetrics, library, 
    x='atac_peak_region_fragments', y='gex_umis_count', cell_col='is_cell',
    xlabel="ATAC peak region fragments", ylabel="mRNA UMI counts",
    ax=None, **kwargs):
    ax = bcmetrics.loc[bcmetrics[cell_col] == 0,:].plot.scatter(
        x=x, y=y, c='lightgrey', label="Non-cells", loglog=True, ax=ax, 
        **kwargs)
    ax = bcmetrics.loc[bcmetrics[cell_col] == 1,:].plot.scatter(
        x=x, y=y, c='firebrick', label="Cells", loglog=True, ax=ax, 
        **kwargs)
    ax.set(xlabel=xlabel, ylabel=ylabel, 
        xlim=[1e1,1e6], ylim=[1e1,1e6])
    n = bcs[library][cell_col].sum()
    ax.set_title(f"Joint Cell Calling ({library}: {n} cells)")
    ax.grid(True)
    return

# %% Functions for rank plots
def rank(bcmetrics, cell_col='is_cell'):
    ## Sort cells from most to fewest reads. 
    reads = bcmetrics[[cell_col, 'gex_raw_reads']].sort_values(
        'gex_raw_reads', ascending=False).reset_index()
    reads['rank'] = reads.index
    ## Sort cells from most to fewest UMIs. 
    umis = bcmetrics[[cell_col, 'gex_umis_count']].sort_values(
        'gex_umis_count', ascending=False).reset_index()
    umis['rank'] = umis.index
    ## Sort cells from most to fewest ATAC fragments. 
    atac = bcmetrics[[cell_col, 'atac_peak_region_fragments']].sort_values(
        'atac_peak_region_fragments', ascending=False).reset_index()
    atac['rank'] = atac.index
    
    return reads, umis, atac

def plot_rank(bcmetrics, library, cellcount=20000, cell_col='is_cell', 
    ax=None, ylim=(1e1, 1e5), plot_thr=True):
    reads, umis, atac = rank(bcmetrics, cell_col=cell_col)
    # for i, label, c in [(0, 'Non-cells', 'grey'), (1, 'Cells', None)]:
    palette = {}
    for data, y, l, c in zip(
        [reads, umis, atac], 
        ['gex_raw_reads', 'gex_umis_count', 'atac_peak_region_fragments'],
        ["Reads", "UMIs", "ATAC fragments"],
        ['tab:'+color for color in ['blue', 'orange', 'green']]):
        ## Reads rank plot
        ax = data.loc[data[cell_col] == 0,].plot.scatter(x='rank', y=y, 
            loglog=True, grid=True, ax=ax, label='Non-cells', c='lightgrey', alpha=0.5, zorder=2)
        ax = data.loc[data[cell_col] == 1,].plot.scatter(x='rank', y=y, 
            loglog=True, grid=True, ax=ax, label='Cells', c=c, alpha=0.5, zorder=2)
        palette[l] = c
    palette['Non-cell'] = 'lightgrey'
    if plot_thr:
        umithr = umis.loc[cellcount, 'gex_umis_count']
        atacthr = atac.loc[cellcount, 'atac_peak_region_fragments']

        ax.axhline(umithr, c='tab:orange')
        ax.axhline(atacthr, c='tab:green')

    ax.set(title=(f"{library} Reads vs Barcodes"), xlabel="Cell Barcodes", ylabel="Read/UMI/ATAC")
    handles = [plt.scatter([], [], color=palette[l], label=l) for l in palette]
    ax.legend(handles=handles, loc='lower left')
    ax.set_xlim(1e3, 1e5)
    ax.set_ylim(*ylim)
    ax.grid(True)
    ## Highlight region or likely real cells.  Should always be under 40k. 
    ax.add_patch( patches.Rectangle((0,0), 40000, 1e6, linewidth=4, 
        edgecolor='none', facecolor='powderblue') )
    ax.add_patch( patches.Rectangle((0,0), cellcount, 1e6, linewidth=4, 
        edgecolor='none', facecolor='lavender') )
    return

def rank_filter(bcmetrics, cellcount, how='outer'):
    reads, umis, atac = rank(bcmetrics)
    readthr = reads.loc[cellcount, 'gex_raw_reads']
    umithr = umis.loc[cellcount, 'gex_umis_count']
    atacthr = atac.loc[cellcount, 'atac_peak_region_fragments']
    print("Droplets containing more than", readthr, "reads")
    print("Droplets containing more than", umithr, "UMIs")
    print("Droplets containing more than", atacthr, "ATAC peak region fragments")
    keepumi = np.array(bcmetrics.gex_umis_count >= umithr)
    keepatac = np.array(bcmetrics.atac_peak_region_fragments >= atacthr)
    if how=='outer':
        bcmetrics['accept'] = np.array((keepumi | keepatac), dtype=int)
    if how=='inner':
        bcmetrics['accept'] = np.array((keepumi & keepatac), dtype=int)
    print(bcmetrics['accept'].sum(), "cells")
    return umithr, atacthr

def custom_jointcall(bcmetrics, library, cellcount, how):
    fig, ax = plt.subplots(figsize=(4,4))
    fig.patch.set_alpha(1.0)
    umithr, atacthr = rank_filter(bcmetrics, cellcount, how=how)
    plot_rank(bcmetrics, library, cellcount, ax=ax, cell_col='accept')
    plt.savefig(outdir / f"rankplot_custom_{how}_{library}.png")
    fig, ax = plt.subplots(figsize=(4,4))
    plot_cellranger_jointcall(bcmetrics, library, ax=ax, cell_col='accept')
    plt.savefig(outdir / f"jointcall_custom_{how}_{library}.png")
    keepbc = bcmetrics.loc[bcmetrics.accept == 1, :]
    keepbc.to_csv(outdir / f"accepted_barcodes_{how}_{library}.csv", 
        columns=["gex_barcode"], index=False, header=False)
    with open(outdir / f"threshold_{how}_{library}.tsv", "w") as f:
        f.write(f"gex_umis_count\t{umithr}\n")
        f.write(f"atac_peak_region_fragments\t{atacthr}\n")
    return

#%%
libs = [f"R{i}" for i in range(1,9)]
bcs = dict()
for library in libs:
    path = datadir / f"{library}/outs/per_barcode_metrics.csv"
    # Read barcode metrics file
    bcmetrics = pd.read_csv(path)
    bcs[library] = bcmetrics

#%% Read Metrics
for library in libs:
    print(library, "barcode count:", len(bcs[library]), '\t', 
        "cell count:", bcs[library]['is_cell'].sum())

# %% Visualize original cellranger-based joint cell calling
for library in bcs:
    fig, ax = plt.subplots(figsize=(4,4))
    plot_cellranger_jointcall(bcs[library], library, ax=ax)
    plt.savefig(outdir / f"jointcall_cellranger_{library}.png")

#%% Visualize original cellranger-based joint cell calling
for library in bcs:
    fig, ax = plt.subplots(figsize=(4,4))
    plot_cellranger_jointcall(bcs[library], library, ax=ax, 
        x="atac_peak_region_cutsites", xlabel="ATAC transposition events in peaks")
    plt.savefig(outdir / f"crossSens_cellranger_{library}.png")

#%% Plot rank plots colored by cell ranger based joint cell calling 
for library in bcs:
    fig, ax = plt.subplots(figsize=(4,4))
    fig.patch.set_alpha(1.0)
    plot_rank(bcs[library], library, cellcount=bcs[library]['is_cell'].sum(), ax=ax)
    plt.savefig(outdir / f"rankplot_cellranger_{library}.png")

#%% Custom joint-calling using umi and fragment thresholds
for how in ['inner', 'outer']:
    library = "R2"
    custom_jointcall(bcs[library], library, cellcount=24000, how=how)
    library = "R4"
    custom_jointcall(bcs[library], library, cellcount=25000, how=how)
    library = "R8"
    custom_jointcall(bcs[library], library, cellcount=33000, how=how)
    for library in ["R1", "R3", "R5", "R6", "R7"]:
        custom_jointcall(bcs[library], library, 
            cellcount=bcs[library]['is_cell'].sum(), how=how)
# %%
#%% Check barcodes
# for lib in libs:
#     df1 = pd.read_csv(
#         f"../output/cell_calling_v1/accepted_barcodes_outer_{lib}.csv", 
#         header=None)
#     df2 = pd.read_csv(
#         f"../output/cell_calling/accepted_barcodes_outer_{lib}.csv", 
#         header=None)
#     assert((df1[0]==df2[0]).all())