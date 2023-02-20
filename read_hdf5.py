# U01_data.read_hdf5 created by bathy at 2/19/2023
import h5py
import numpy as np
import pandas as pd
import time

start= time.time()

filename = r"C:\Users\bathy\Downloads\12be772d-49bd-4abc-b08e-76684c2ef1f2_KPMPAtlas_PREMIER_062921.h5Seurat"
fo = h5py.File(filename, 'r')

data_array = np.array(fo['assays']['RNA']['scale.data'], dtype=np.float64)
print(time.time()-start); start=time.time()

gene_list = [i.decode('utf8') for i in fo['assays']['RNA']['scaled.features']]
meta_keys = fo['meta.data'].keys()
umap1_key = 'UMAP_1' if 'UMAP_1' in meta_keys else 'umap1'
umap2_key = 'UMAP_2' if 'UMAP_2' in meta_keys else 'umap2'

umap1 = np.array(fo['meta.data'][umap1_key], dtype=np.float64)
umap2 = np.array(fo['meta.data'][umap2_key], dtype=np.float64)

df1 = pd.DataFrame(data=data_array, index=[i.decode('utf8') for i in fo['cell.names']], columns=gene_list)

df1['umap1'] = umap1
df1['umap2'] = umap2
df1.reset_index().to_feather('kpmp_atlas_data.feather')
