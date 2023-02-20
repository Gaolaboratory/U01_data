# U01_data.cluster_map created by bathy at 2/19/2023
# video_detection.cluster_map created by bathy at 11/11/2021

import pandas as pd
import matplotlib.pyplot as plt
import time
import seaborn as sns
from data_store import exist_gene, matrisome_dict
import numpy as np


def plot_cluster(csv_name):
    sns.set(rc={'figure.figsize': (10, 15)})
    df1 = pd.read_csv(csv_name, index_col=0)
    plt.clf()
    ax = sns.clustermap(df1, cmap='viridis')
    for tick in ax.ax_heatmap.get_xticklabels():
        tick.set_rotation(20)
    ax.ax_heatmap.tick_params(labelsize=15)
    plt.show()


def seaborn_plot(gene_list, filename = 'kpmp_atlas_data.feather', out_name='test.png'):
    plt.clf()
    sns.set(rc={'figure.figsize': (12, 12)})
    sns.set_style("white")
    sns.set_context("poster")
    start = time.time()
    gene_list=[i for i in gene_list if i in exist_gene]
    gene_list.append('umap1')
    gene_list.append('umap2')
    df_data = pd.read_feather(filename, columns=gene_list)
    median_score = []
    for idx, row in df_data.iterrows():
        score = np.max(row[gene_list])
        median_score.append(score)
    df_data['score'] = median_score
    df_data = df_data.sort_values(by='score', ascending=True)
    print(time.time() - start); start = time.time()
    ax = sns.scatterplot(data=df_data, x='umap1', y='umap2', hue='score', s=3, legend=False, palette='viridis')
    ax.set(yticklabels=[], xticklabels=[], xlabel='', ylabel='')  # remove the tick labels
    ax.tick_params(left=False, bottom=False)# remove the ticks
    #plt.text(-15,11, "Secreted", fontdict={'size':70})
    figure = plt.gcf()
    figure.set_size_inches(12, 12)
    plt.savefig(out_name, dpi=150, bbox_inches='tight')
    #plt.show()


def range_plot(gene_list, filename = 'kpmp_atlas_data.feather', out_name='test.png', range=(-15,15,-15,15)):
    gene_list=[i for i in gene_list if i in exist_gene]
    plt.clf()
    count_dict={}
    for each in gene_list:
        count_dict[each]=0
    sns.set(rc={'figure.figsize': (12, 12)})
    sns.set_style("white")
    sns.set_context("poster")
    start = time.time()
    gene_list.append('umap1')
    gene_list.append('umap2')
    df_data = pd.read_feather(filename, columns=gene_list)
    df_data = df_data[df_data['umap1'].between(range[0],range[1])]
    df_data = df_data[df_data['umap2'].between(range[2], range[3])]
    median_score = []
    for idx, row in df_data.iterrows():
        score = np.sum(row[gene_list])
        median_score.append(score)
        for each in gene_list:
            if not each.startswith('umap'):
                count_dict[each]+=row[each]
    print(count_dict)
    df_data['score'] = median_score
    df_data = df_data.sort_values(by='score', ascending=True)
    print(time.time() - start); start = time.time()
    ax = sns.scatterplot(data=df_data, x='umap1', y='umap2', hue='score', s=5, legend=False, palette='viridis')
    ax.set(yticklabels=[], xticklabels=[], xlabel='', ylabel='')  # remove the tick labels
    ax.tick_params(left=False, bottom=False)# remove the ticks
    #plt.text(-15,11, "Secreted", fontdict={'size':70})
    figure = plt.gcf()
    figure.set_size_inches(12, 12)
    plt.savefig(out_name, dpi=150, bbox_inches='tight')


if __name__ == '__main__':
    range_plot(['COL4A1', 'COL4A2', 'COL4A4', 'COL6A3', 'COL15A1', 'COL5A1', 'COL12A1', 'COL6A1', 'COL4A6', 'COL6A2', 'COL4A3', 'COL1A2'],
               out_name='a.png',
               range=(5,12,-5,0))

    # for key in matrisome_dict:
    #     range_plot(matrisome_dict[key], out_name=key.replace(' ','_')+'.png', range=(5,12,-5,-0))

    #plot_cluster('median.csv')