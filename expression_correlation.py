# U01_data.expression_correlation created by bathy at 2/19/2023
# video_detection.expression_correlation created by bathy at 11/11/2021
import pandas as pd
from scipy import spatial
from data_store import matrisome_dict, annotation_l1, exist_gene
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch


def cluster_corr(corr_array, inplace=False):
    """
    Rearranges the correlation matrix, corr_array, so that groups of highly
    correlated variables are next to eachother

    Parameters
    ----------
    corr_array : pandas.DataFrame or numpy.ndarray
        a NxN correlation matrix

    Returns
    -------
    pandas.DataFrame or numpy.ndarray
        a NxN correlation matrix with the columns and rows rearranged
    """
    pairwise_distances = sch.distance.pdist(corr_array)
    linkage = sch.linkage(pairwise_distances, method='complete')
    cluster_distance_threshold = pairwise_distances.max() / 2
    idx_to_cluster_array = sch.fcluster(linkage, cluster_distance_threshold,
                                        criterion='distance')
    idx = np.argsort(idx_to_cluster_array)

    if not inplace:
        corr_array = corr_array.copy()

    if isinstance(corr_array, pd.DataFrame):
        return corr_array.iloc[idx, :].T.iloc[idx, :]
    return corr_array[idx, :][:, idx]


def get_expression(gene_name_list, filename = 'kpmp_atlas_data.feather',):
    result={}
    gene_name_list = [i for i in gene_name_list if i in exist_gene]
    df_data = pd.read_feather(filename, columns=gene_name_list)
    for each in gene_name_list:
        result[each]=df_data[each].to_numpy()
    return result


if __name__ == '__main__':
    matrix_data={}
    annotation_data={}
    for each in matrisome_dict:
        matrix_data[each]=get_expression(matrisome_dict[each])
    for each in annotation_l1:
        annotation_data[each] = get_expression(annotation_l1[each])

    correlation_matrix_max = np.zeros([len(matrix_data), len(annotation_data)])
    correlation_matrix_mean = np.zeros([len(matrix_data), len(annotation_data)])
    correlation_matrix_median = np.zeros([len(matrix_data), len(annotation_data)])
    i=0
    for each in matrix_data:
        j = 0
        for each2 in annotation_data:
            score_array=[]
            for x in matrix_data[each]:
                for y in annotation_data[each2]:
                    score =1-spatial.distance.cosine(matrix_data[each][x], annotation_data[each2][y])
                    score_array.append(score)
            correlation_matrix_max[i,j]=max(score_array)
            correlation_matrix_mean[i, j] = np.mean(score_array)
            correlation_matrix_median[i, j] = np.median(score_array)
            j+=1
        i+=1

    print(correlation_matrix_max, correlation_matrix_mean, correlation_matrix_median)
    df1 = pd.DataFrame(correlation_matrix_max.T,  index= annotation_data.keys(), columns=matrix_data.keys())
    df1.to_csv('max.csv')
    df2 = pd.DataFrame(correlation_matrix_mean.T, index=annotation_data.keys(), columns=matrix_data.keys())
    df2.to_csv('mean.csv')
    df3 = pd.DataFrame(correlation_matrix_median.T, index=annotation_data.keys(), columns=matrix_data.keys())
    df3.to_csv('median.csv')
