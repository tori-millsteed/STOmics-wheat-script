
"""
Created January-March 2025

Author: Tori Millsteed
t.millsteed@uq.edu.au
"""

#STOmics wheat analysis

#Bin50, filtered min 30 MID, Normalised by log1P method, spatial Leiden clusters

#First the 'st' conda environment was created to work from

#python3.8.5 is the best version for stereopy
module load python3/3.8.5
cd /data/path/where/st/environment/is/stored
source st/bin/activate

#Once inside the environmnet, install ipython and stereopy
pip install ipython
pip install stereopy

#Open python terminal
python3

#From here to 'save as h5ad' the code is taken from the Stereopy library basic workflow for square bin sample
import os
from natsort import natsorted
import stereo as st
from stereo.core.ms_data import MSData
from stereo.core.ms_pipeline import slice_generator

#Read in the data
data_path = '/data/path/to/file.gef'
st.io.read_gef_info(data_path)

#Select bin level from Bin10-Bin200
data = st.io.read_gef(file_path=data_path, bin_size=50)

#Change data to anndata object 
data=data.to_ann_based()

#Save the raw data
data.tl.raw_checkpoint()
data.tl.raw

#Filter by minimum MID count per cell (bin) to reduce background noise
data.tl.filter_cells(min_counts=30)

#Filter by coordinates to focus in on the part of the chip you are interested in
data.tl.filter_coordinates(min_x=7176, max_x=16425, min_y=5300, max_y=12200)

#Normalise the data
data.tl.normalize_total(target_sum=None)
data.tl.log1p()

#Calculate highly variable genes
data.tl.highly_variable_genes(
        min_mean=0.0125,
        max_mean=3,
        min_disp=0.5,
        n_top_genes=2000,

        res_key='highly_variable_genes'
        )
data.plt.highly_variable_genes(res_key='highly_variable_genes', out_path='/outpath/to/file.png')

#Scale data to improve processing efficiency
data.tl.scale(max_value=10, zero_center=False)

#Calculate Principle component analysis
data.tl.pca(
        use_highly_genes=False,
        n_pcs=30,
        res_key='pca'
        )
data.tl.key_record
data.plt.elbow(pca_res_key='pca', out_path='/outpath/to/file.png')

#Generate neighbours
data.tl.neighbors(
        pca_res_key='pca',
        n_pcs=30,
        res_key='neighbors'
        )

#Generate spatial neighbours
data.tl.spatial_neighbors(
        neighbors_res_key='neighbors',
        res_key='spatial_neighbors'
        )
data.tl.spatial_neighbors

#Generate umap
data.tl.umap(pca_res_key='pca', neighbors_res_key='spatial_neighbors', res_key='umap')

#Generate neighbours Leiden clusters
data.tl.leiden(neighbors_res_key='neighbors', res_key='leiden')

#Plot neighbours Leiden clusters
data.plt.cluster_scatter(res_key='leiden', out_path='/outpath/to/file.png')

#Plot umap
data.plt.umap(res_key='umap', cluster_key='spatial_leiden', out_path='/outpath/to/file.png') 

#Generate spatial neighbours Leiden clusters
data.tl.leiden(neighbors_res_key='spatial_neighbors', res_key='spatial_leiden')

#Plot spatial neighbours Leiden clusters
data.plt.cluster_scatter(res_key='spatial_leiden', out_path='/outpath/to/file.png')

#Generate marker genes for the 'spatial leiden' clusters
data.tl.find_marker_genes(
        cluster_res_key='spatial_leiden',
        method='t_test',
        use_highly_genes=False,
        use_raw=True
        )

#Plot marker genes
data.plt.marker_genes_scatter(res_key='marker_genes', markers_num=5, out_path='/outpath/to/file.png')

#Save processed data object, with normalisation, filtering and cluster data, as h5ad
adata = st.io.stereo_to_anndata(data, flavor='scanpy', output='/outpath/to/file.h5ad')



#CODE BELOW IS TO GENERATE RANDOMLY SAMPLED BINS FOR PAIRWISE PERMUTATION TESTS BETWEEN TARGET GENES AND CLUSTERS

#Set parameters to generate randomly sampled bins
#downsample size can be used to set the number of bins generated (we used a fraction instead in code below)
#set random seed so results can be replicated
cluster_IDs='spatial_leiden'
downsample_size=78
random_seed=222

#Randomely generate sample bins as 5% of bins froms each cluster
sample_bins=data._ann_data.obs.groupby([cluster_IDs]).apply(lambda x: x.sample(frac=0.05, random_state=random_seed)).index.tolist()

#Read the bins generated
sample_bins

#Code for saving the cluster and bin IDs generated. This is so that we can use the locations of the bins and clusters from the normalised data to pull out the read values from those locations for the unnormalised data

import pandas as pd
cluster_and_bin_df = pd.DataFrame(sample_bins, columns = ['cluster_id', 'bin_id']) #This tells it to put the list of tuples into a data frame, with the column names for clarity
cluster_and_bin_df.to_csv("cluster_and_bin_ids_normalise_gef.csv", index=False)
#There will now be a file called "cluster_and_bin_ids_normalise_gef.csv" in the same directory as everything else
#You only need to run this once per sample/gene set, and then it will be saved 

# Get the bin names for indexing 
sample_bins=[index[1] for index in sample_bins]


#THEN READ IN THE ORIGINAL, UNNORMALISED GEF FILE AGAIN

#Read in the data
data_path = '/data/path/to/file.gef'
st.io.read_gef_info(data_path)

#Select bin level from Bin10-Bin200
data = st.io.read_gef(file_path=data_path, bin_size=50)

#Change data to anndata object 
data=data.to_ann_based()

#Use this code to retrieve the bin and cluster list from the normalised data
import pandas as pd
sample_bins = list(pd.read_csv("cluster_and_bin_ids_normalise_gef.csv", dtype={'bin_id': str})['bin_id']) #Opens the data frame, sets the bin_id column values to strings (otherwise they are ints and it breaks other bits of code), then pulls out the bin_id column, and then turns it into list format, which gives us the same format as the code above


# Subset the data by bin names of downsample
data_subset_1=data.sub_by_name(sample_bins)


# Define target gene list for future subsetting
genes_of_interest = ('LOC123101925', 'LOC543308', 'LOC100125699', 'LOC123057832', 'LOC123064838', 'LOC123073994', 'LOC123187748', 'LOC123043995', 'LOC123051864', 'LOC123104584', 'LOC123179983', 'LOC543476', 'LOC543084', 'LOC123182837', 'LOC123185730', 'LOC123041664', 'LOC123049628', 'LOC123055206', 'LOC123132411', 'LOC123181969', 'LOC123188447', 'LOC123044698', 'LOC542962')

# Subset the data by the defined gene list
data_subset_2= data_subset_1.tl.filter_genes(gene_list=genes_of_interest)

# Covert the resultant data to pandas dataframe if needed
data_subset_df=data_subset_2.to_df()

#Save data frame as csv 
data_subset_df.to_csv('/outpath/to/file.csv', index=True)
