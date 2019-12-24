
# python plotHiC_withTADs_withCols.py

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



#def load_tad_list(tad_file, with_labs, **args):
#    tad_dt = pd.read_csv(tad_file, **args)
#    
#    if tad_dt.shape[1] - int(with_labs) == 3:
#        #return([np.array(tad_dt.iloc[:,1]), np.array(tad_dt.iloc[:,2])])
#        tad_pos = [list(tad_dt.iloc[:,1]), list(tad_dt.iloc[:,2])]
#    elif tad_dt.shape[1] - int(with_labs) == 2:
#        #return([np.array(tad_dt.iloc[:,0]), np.array(tad_dt.iloc[:,1])])
#        tad_pos = [list(tad_dt.iloc[:,0]), list(tad_dt.iloc[:,1])]
#    else:
#            return None
#        
#    if with_labs:
#        tad_pos.append(list(tad_dt.iloc[:,tad_dt.shape[1]-1]))
#    return tad_pos

# all expected to have chromo-start-end-label-color

def load_tad_list(tad_file, with_labs, **args):
    tad_dt = pd.read_csv(tad_file, **args)
    return [list(tad_dt.iloc[:,1]), list(tad_dt.iloc[:,2]), list(tad_dt.iloc[:,3]), list(tad_dt.iloc[:,4])]




# SETTINGS FOR HI-C MATRIX
#matrix_file = "/media/electron//mnt/ptemp/luca/hic_enhance/data/hic/Rao2014-GM12878-full/Rao2014_observed_NONE_1_1_50000.txt"
matrix_file = "mat_chr11_40kb_ob.txt"
resolution = 40000

# SETTINGS FOR TAD DATA
tad_file = "ENCSR489OCU_NCI-H460_40kb_coord11.txt"  # can be 2-col or 3-col
plot_labs = True

# SET WHAT TO PLOT:
chrom = 'chr11'
map_start=102200001
map_end=103320000


# OUTFILE 
output_file = "./somewhere.pdf"


# SOME PARAMETERS FOR PLOTTING
shrink = 1
lwd = 5

# LOAD AND SUBSET Hi-C DATA
data = pd.read_csv(matrix_file, sep="\t", header=None, names=['bin1', 'bin2', 'value'])
data = data[(data.bin1>=map_start) & (data.bin2<=map_end)]
m = to_matrix(data, chrom, chrom, resolution, start1=map_start, end1=map_end, start2=map_start, end2=map_end)

# LOAD THE LIST OF TADS TO PLOT
tad_toplot = load_tad_list(tad_file, with_labs=True, sep="\t", header=None)


# retrieve the height of the plot (half-diagonal)
diagonal_height = np.sqrt(m.shape[0]**2 / 2)


my_map = plt.get_cmap("Reds")
my_map.set_under('white')
#fig, ax = plt.subplots(1, 1, figsize=(40, 40))
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
X = np.log10(np.triu(m.todense()*shrink) + 1)
ax.matshow( rotate(X, 45), cmap=my_map, vmin=0.01 )

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

    tadTri = plt.Polygon(tad_triangle, facecolor='none', edgecolor=tad_col, linewidth=lwd)
    ax.add_patch(tadTri)

    plt.text(0.5*(conv_start+conv_end),tad_height,tad_lab, fontsize=10, horizontalalignment='center', verticalalignment='bottom')

    
#if plot_labs:
#    for tad_start, tad_end, tad_lab in zip(tad_toplot[0], tad_toplot[1], tad_toplot[2]):
#        tad_start_bin = (tad_start - map_start)//resolution  
#        tad_end_bin = (tad_end - map_start)//resolution

#        conv_start = np.sqrt(2) * tad_start_bin  # distance to bin start: on the side of the square, it is the # of bins but in triangle plot, it is the diagonal
#        conv_end = np.sqrt(2) * tad_end_bin
#        

#        tad_height = diagonal_height-(conv_end-conv_start)*0.5  # needs the "diagonal_heigh" because the way it is plotted and the orientation of the y-axis
#        
#        plt.text(0.5*(conv_start+conv_end),tad_height,tad_lab, fontsize=10, horizontalalignment='center', verticalalignment='bottom')
#        #
#        # plt.text(0.5*(conv_start+conv_end),tad_height,tad_lab, fontsize=10)
    

if output_file:
	plt.savefig(output_file)
	print("... saved: " + output_file + "\n")
    

    


