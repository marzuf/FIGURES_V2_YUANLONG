
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

path = "/home/luca/projects/GenomicInteractions/data/hic/Rao2014-GM12878-full/Rao2014_observed_NONE_1_1_50000.txt"
output_path = "./somewhere.pdf"
chrom = 'chr1'
start=1000000
end=70000000

tad_start = 1100000
tad_end = 12500000

resolution = 50000
shrink=1

tad_pos_1 = (tad_start - start)//resolution
tad_pos_2 = (tad_end - start)//resolution

data = pd.read_csv(path, sep="\t", header=None, names=['bin1', 'bin2', 'value'])
data = data[(data.bin1>=start) & (data.bin2<=end)]


m = to_matrix(data, chrom, chrom, resolution, start1=start, end1=end, start2=start, end2=end)
diagonal_height = np.sqrt(m.shape[0]**2 / 2)

triangle1 = np.array([
        [0,diagonal_height],
        [np.sqrt(2)*tad_pos_1/2, diagonal_height - np.sqrt(2)*tad_pos_1/2],
        [np.sqrt(2)*tad_pos_1, diagonal_height]
    ])

my_map = plt.get_cmap("Reds")
my_map.set_under('white')
fig, ax = plt.subplots(1, 1, figsize=(40, 40))
X = np.log10(np.triu(m.todense()*shrink) + 1)
ax.matshow( rotate(X, 45), cmap=my_map, vmin=0.01 )
t1 = plt.Polygon(triangle1, facecolor='none', edgecolor='black')
ax.add_patch(t1)
plt.ylim(diagonal_height)
plt.xticks([])
plt.yticks([])
plt.axis('off')