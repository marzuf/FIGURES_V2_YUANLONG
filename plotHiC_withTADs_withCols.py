
# python plotHiC_withTADs_withCols.py

import os
import sys
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import rotate
from scipy.sparse import coo_matrix, triu

def to_matrix(data, chrom1, chrom2, 
              resolution, assembly="hg19",
              start1=None, end1=None, 
              start2=None, end2=None):
    
    if all([x is not None for x in [start1, end1, start2, end2]]):
        start_bin1 = start1//resolution
        end_bin1 = end1//resolution
        start_bin2 = start2//resolution
        end_bin2 = end2//resolution
    else:
        chromsize = load_chromsizes(assembly)
        chromsize1 = chromsize[chrom1]
        chromsize2 = chromsize[chrom2]
        start_bin1 = 0
        end_bin1 = chromsize1//resolution
        start_bin2 = 0
        end_bin2 = chromsize2//resolution
    
    n_bins_1 = end_bin1 - start_bin1 + 1#int(chromsize1 / resolution) + 1
    n_bins_2 = end_bin2 - start_bin2 + 1#int(chromsize2 / resolution) + 1
    
    bin1s = data.bin1//resolution - start_bin1#(data.bin1 / resolution).astype(int).values
    bin2s = data.bin2//resolution - start_bin2#(data.bin2 / resolution).astype(int).values
    values = data.value.values
    m = coo_matrix( (values, (bin1s, bin2s)), shape=(n_bins_1, n_bins_2) )
    if chrom1 == chrom2:
        m = triu(m, k=1).T + m
    return m

expr_ds = None
# python plotHiC_withTADs_withCols.py ENCSR489OCU_NCI-H460_40kb chr11_TAD390 mega_ENCSR489OCU_NCI-H460_mat_chr11_40kb_ob.txt
myargs = sys.argv
if(len(myargs) == 5):
    hic_ds = myargs[1]
    expr_ds = myargs[2]
    select_tad = myargs[3]
    matrix_file = myargs[4]
elif(len(myargs) == 4):
    hic_ds = myargs[1]
    select_tad = myargs[2]
    matrix_file = myargs[3]
else:    
    hic_ds = "ENCSR489OCU_NCI-H460_40kb"
    expr_ds = "TCGAlusc_norm_lusc"
    matrix_file = "mega_ENCSR489OCU_NCI-H460_mat_chr10_40kb_ob.txt"
    select_tad = "chr10_TAD268"
    #matrix_file = "mega_ENCSR489OCU_NCI-H460_mat_chr11_40kb_ob.txt"
    #select_tad = "chr11_TAD390"


# SOME PARAMETERS FOR PLOTTING
other_col_tad = "black"
select_col_tad = "green"
#other_col_lab = "black"
#select_col_lab = "green"
nAround_toplot = 0
nSurroundBins = 2
shrink = 1
tad_lwd = 1.2
lab_offset = -0.8
labSize = 10    
labBox = True
addGenes = True
genesOffset = 0.5
genesSpacing = 0.5
symbolOffset = 0.1
gene_col = "blue"
gene_lwd = 1
plot_labs = True
plotTit = hic_ds
if expr_ds:
    plotTit += " - " + expr_ds

if addGenes:
    assert expr_ds

# PARAMETERS TO SET 
setDir = "/media/electron"
setDir = ""
resolution = 40000
tad_file  = setDir + "/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/" + hic_ds + "/genes2tad/all_assigned_regions.txt"
assert os.path.exists(tad_file)

if addGenes:
    gene2tad_file  = setDir + "/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/" + hic_ds + "/genes2tad/all_genes_positions.txt"
    assert os.path.exists(gene2tad_file)
    entrezDT_file = setDir + "/mnt/ed4/marie/entrez2synonym/entrez/ENTREZ_POS/gff_entrez_position_GRCh37p13_nodup.txt"
    assert os.path.exists(entrezDT_file)
    pipGene_file = setDir + "/mnt/etemp/marie/v2_Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER/" + hic_ds + "/" + expr_ds + "/0_prepGeneData" + "/pipeline_geneList.txt"
    assert os.path.exists(pipGene_file)

out_dir = "PLOTHIC_WITHTADS_WITHCOLS"
output_file =  out_dir + "/" + hic_ds + "_" + select_tad + ".pdf"


if not os.path.exists(out_dir):
    os.makedirs(out_dir)

if addGenes:
    gene2tad_dt = pd.read_csv(gene2tad_file, sep="\t", header=None, names = ['entrezID', 'chromo', 'start', 'end', 'region'] )
    tad_genes_dt = gene2tad_dt[gene2tad_dt['region'] == select_tad]
    gff_dt = pd.read_csv(entrezDT_file, sep="\t")
    gff_dt = gff_dt[['entrezID', 'symbol']]
    tad_genes_symbols_dt = tad_genes_dt.merge(gff_dt, on='entrezID', how='left')        
    pipGenes_dt = pd.read_csv(pipGene_file, sep="\t")
    tad_genes_symbols_dt = tad_genes_symbols_dt[tad_genes_symbols_dt['entrezID'].isin(list(pipGenes_dt['value']))].reset_index(drop=True)
    assert tad_genes_symbols_dt.shape[0] > 0
    
    
tad_dt = pd.read_csv(tad_file, sep="\t", header=None, names = ['chromo', 'region', 'start', 'end'] )

select_idxs=  tad_dt[tad_dt['region'] == select_tad].index
assert len(select_idxs) == 1
select_idx =  int(tad_dt[tad_dt['region'] == select_tad].index[0]) 

chrom = tad_dt['chromo'][select_idx]
assert re.match(chrom, select_tad)

#start_pos_list = list(tad_dt['start'][(select_idx-nAround_toplot):(select_idx+nAround_toplot+1)])
start_pos_list = list(tad_dt['start'][(select_idx-nAround_toplot):(select_idx+nAround_toplot+1)] - 1)
end_pos_list = list(tad_dt['end'][(select_idx-nAround_toplot):(select_idx+nAround_toplot+1)])
#end_pos_list = list(tad_dt['end'][(select_idx-nAround_toplot):(select_idx+nAround_toplot+1)] + 1)
lab_list = ([""] * nAround_toplot) + [select_tad] + ([""] * nAround_toplot)
col_list_tad = ([other_col_tad] * nAround_toplot) + [select_col_tad] + ([other_col_tad] * nAround_toplot)
#col_list_lab = ([other_col_lab] * nAround_toplot) + [select_col_lab] + ([other_col_lab] * nAround_toplot)



map_start = start_pos_list[0] - nSurroundBins*resolution
map_end = end_pos_list[len(end_pos_list)-1] + nSurroundBins*resolution

assert len(start_pos_list) == len(end_pos_list) == len(lab_list) == len(col_list_tad)

tad_toplot = [
    start_pos_list,
    end_pos_list,
    lab_list,
    col_list_tad
]
    
# LOAD AND SUBSET Hi-C DATA
data = pd.read_csv(matrix_file, sep="\t", header=None, names=['bin1', 'bin2', 'value'])
data = data[(data.bin1>=map_start) & (data.bin2<=map_end)]
m = to_matrix(data, chrom, chrom, resolution, start1=map_start, end1=map_end, start2=map_start, end2=map_end)

# retrieve the height of the plot (half-diagonal)
diagonal_height = np.sqrt(m.shape[0]**2 / 2)

my_map = plt.get_cmap("Reds")
my_map.set_under('white')
#fig, ax = plt.subplots(1, 1, figsize=(40, 40))
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
X = np.log10(np.triu(m.todense()*shrink) + 1)
ax.matshow( rotate(X, 45), cmap=my_map, vmin=0.01 )

if addGenes:
    plt.ylim(diagonal_height )
    #plt.ylim(diagonal_height + genesOffset + genesSpacing*tad_genes_symbols_dt.shape[0])
else:
    plt.ylim(diagonal_height)

plt.xticks([])
plt.yticks([])
plt.axis('on')


for tad_start, tad_end , tad_lab, tad_col in zip(tad_toplot[0], tad_toplot[1], tad_toplot[2], tad_toplot[3]):

    tad_start_bin = (tad_start - map_start)//resolution  
    tad_end_bin = (tad_end - map_start)//resolution


    conv_start = np.sqrt(2) * tad_start_bin  # distance to bin start: on the side of the square, it is the # of bins but in triangle plot, it is the diagonal
    conv_end = np.sqrt(2) * tad_end_bin

    tad_height = diagonal_height-(conv_end-conv_start)*0.5  # needs the "diagonal_heigh" because the way it is plotted and the orientation of the y-axis
    
    tad_triangle = np.array([
        [conv_start,diagonal_height],
        [0.5*(conv_start+conv_end), tad_height],
        [conv_end, diagonal_height]
    ])

    tadTri = plt.Polygon(tad_triangle, facecolor='none', edgecolor=tad_col, linewidth=tad_lwd)
    ax.add_patch(tadTri)

    if labBox:
        plt.text(0.5*(conv_start+conv_end),tad_height + lab_offset,tad_lab, fontsize=labSize, horizontalalignment='center', verticalalignment='bottom',
             bbox=dict(facecolor=tad_col, alpha=0.5), fontweight='bold')
    else:
        plt.text(0.5*(conv_start+conv_end),tad_height + lab_offset,tad_lab, fontsize=labSize, horizontalalignment='center', verticalalignment='bottom', fontweight='bold')
    
    
if addGenes:
    from matplotlib.patches import Rectangle
    from matplotlib import collections  as mc
    #lines = [[(0, 50), (20, 50)], [(50, 50), (70, 50)]]
    #c = np.array([(1, 0, 0, 1), (0, 1, 0, 1), (0, 0, 1, 1)])
    #c = ['black', 'green']
    #lc = mc.LineCollection(lines, colors=c, linewidths=2)
    #fig, ax = pl.subplots()
    #ax.add_collection(lc)
    genePos = diagonal_height.item() + genesOffset
    
    
#     first_start = start_pos_list[0]
#     first_end =  0.5*( start_pos_list[0] + end_pos_list[len(end_pos_list)-1])
#     first_lab = "FOO_FIRST_GENE"
#     conv_gene_start = np.sqrt(2) * (first_start - map_start)/resolution
#     conv_gene_end = np.sqrt(2) * (first_end - map_start)/resolution
#     print("FIRST");print(conv_gene_start);print(genePos);print(conv_gene_end-conv_gene_start)
#     ax.add_patch(Rectangle((conv_gene_start, genePos), (conv_gene_end-conv_gene_start), height=0.2,
#                            facecolor=gene_col, clip_on=False))
#     plt.text(0.5*(conv_gene_start+conv_gene_end), genePos-symbolOffset, first_lab, fontsize=8, horizontalalignment='center', verticalalignment='center')
#     genePos += genesSpacing
        

    
    
    for i in range(tad_genes_symbols_dt.shape[0]):
        conv_gene_start = np.sqrt(2) * (tad_genes_symbols_dt['start'][i] - map_start)/resolution
        conv_gene_end = np.sqrt(2) * (tad_genes_symbols_dt['end'][i] - map_start)/resolution
        #plt.plot([conv_gene_start, genePos], [conv_gene_end, genePos], color=gene_col, lw=gene_lwd)
        #lines = [[(conv_gene_start, genePos), (conv_gene_end, genePos)]]
        #c = [gene_col]
        #lc = mc.LineCollection(lines, colors=c, linewidths=2)
        #ax.add_collection(lc)
        #plt.add_collection(lc)
        #plt.text(conv_gene_start, genePos,tad_genes_symbols_dt['symbol'][i], fontsize=10, horizontalalignment='right', verticalalignment='center', fontweight='bold')
        print("i="+str(i));print(conv_gene_start);print(genePos);print(conv_gene_end-conv_gene_start)
        ax.add_patch(Rectangle((conv_gene_start, genePos), (conv_gene_end-conv_gene_start), height=0.2,
                               facecolor=gene_col, clip_on=False))
        plt.text(0.5*(conv_gene_start+conv_gene_end), genePos-symbolOffset,tad_genes_symbols_dt['symbol'][i], fontsize=8, horizontalalignment='center', verticalalignment='center')

        genePos += genesSpacing
        
        
        
        
#     last_start = 0.5*( start_pos_list[0] + end_pos_list[len(end_pos_list)-1])
#     last_end =  end_pos_list[len(end_pos_list)-1]
#     last_lab = "FOO_LAST_GENE"
#     conv_gene_start = np.sqrt(2) * (last_start - map_start)/resolution
#     conv_gene_end = np.sqrt(2) * (last_end - map_start)/resolution
#     print("LAST");print(conv_gene_start);print(genePos);print(conv_gene_end-conv_gene_start)
#     ax.add_patch(Rectangle((conv_gene_start, genePos), (conv_gene_end-conv_gene_start), height=0.2,
#                            facecolor=gene_col, clip_on=False))
#     plt.text(0.5*(conv_gene_start+conv_gene_end), genePos-symbolOffset, last_lab, fontsize=8, horizontalalignment='center', verticalalignment='center')
#     genePos += genesSpacing
    
        
        

   
plt.title(plotTit, fontweight="bold")    

if output_file:
	plt.savefig(output_file, bbox_inches='tight',transparent=True)
    
	print("... saved: " + output_file + "\n")
    
