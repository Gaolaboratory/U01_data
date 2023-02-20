# U01_data.scatter_plot created by bathy at 2/19/2023
# video_detection.scatter_plot created by bathy at 11/10/2021
import pandas as pd
from bokeh.plotting import figure, show
from bokeh.transform import linear_cmap, jitter
from bokeh.palettes import Spectral10
from bokeh.io import export_png
import time
from selenium import webdriver
import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm
from data_store import exist_gene
import matplotlib.colors as mcolors


def hex_to_rgb(value):
    '''
    Converts hex to rgb colours
    value: string of 6 characters representing a hex colour.
    Returns: list length 3 of RGB values'''
    value = value.strip("#")  # removes hash symbol if present
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_dec(value):
    '''
    Converts rgb to decimal colours (i.e. divides each value by 256)
    value: list (length 3) of RGB values
    Returns: list (length 3) of decimal values'''
    return [v / 256 for v in value]


def get_continuous_cmap(hex_list, float_list=None):
    ''' creates and returns a color map that can be used in heat map figures.
        If float_list is not provided, colour map graduates linearly between each color in hex_list.
        If float_list is provided, each color in hex_list is mapped to the respective location in float_list.

        Parameters
        ----------
        hex_list: list of hex code strings
        float_list: list of floats between 0 and 1, same length as hex_list. Must start with 0 and end with 1.

        Returns
        ----------
        colour map'''
    rgb_list = [rgb_to_dec(hex_to_rgb(i)) for i in hex_list]
    if float_list:
        pass
    else:
        float_list = list(np.linspace(0, 1, len(rgb_list)))

    cdict = dict()
    for num, col in enumerate(['red', 'green', 'blue']):
        col_list = [[float_list[i], rgb_list[i][num], rgb_list[i][num]] for i in range(len(float_list))]
        cdict[col] = col_list
    cmp = mcolors.LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=10016)
    return cmp


# hex_list = ['#ffffff', '#ffff00', '#ff00ff', '#00df00', '#ff0000', '#000000']
# #hex_list = ['#283F3B', '#556F44', '#659B5E', '#95BF74', '#99DDC8', '#000000']
# plt.register_cmap("mycolormap", get_continuous_cmap(hex_list))
# cpal=sns.color_palette("mycolormap", n_colors=10016)


def to_list(comma_string):
    return comma_string.split(',')


def bokeh_plot(gene_name='COL4A1', filename='kpmp_atlas_data.feather'):
    driver = webdriver.Chrome(r'C:\Users\bathy\Downloads\chromedriver.exe')
    start = time.time()
    df_data = pd.read_feather(filename, columns=[gene_name, 'umap1', 'umap2'])
    df_data = df_data.sort_values(by=gene_name, ascending=False)
    print(time.time() - start);
    start = time.time()
    mapper = linear_cmap(field_name=gene_name, palette=Spectral10, low=min(df_data[gene_name]), high=max(df_data[gene_name]))
    TOOLS = "hover,crosshair,pan,wheel_zoom,zoom_in,zoom_out,box_zoom,undo,redo,reset,tap,save,box_select,poly_select,lasso_select,"
    p = figure(tools=TOOLS, width=2000, height=1250)
    p.scatter(x=jitter('umap1', width=0.5), y=jitter('umap2', width=0.5), radius=0.03, fill_alpha=0.4,
              line_color=None, source=df_data, color=mapper)
    export_png(p, filename='test.png', webdriver=driver)
    show(p)


def seaborn_plot(gene_name='COL4A1', filename='kpmp_biopsy_data.feather', path='pics\\'):
    plt.clf()
    sns.set(rc={'figure.figsize': (12, 12)})
    sns.set_style("white")
    sns.set_context("poster")
    start = time.time()
    df_data = pd.read_feather(filename, columns=[gene_name, 'umap1', 'umap2'])
    df_data = df_data.sort_values(by=gene_name, ascending=True)
    print(time.time() - start, gene_name);
    start = time.time()
    ax = sns.scatterplot(data=df_data, x='umap1', y='umap2', hue=gene_name, s=3, legend=False, palette="magma", x_jitter=True, y_jitter=True, alpha=1)
    ax.set(yticklabels=[], xticklabels=[], xlabel='', ylabel='')  # remove the tick labels
    ax.tick_params(left=False, bottom=False)  # remove the ticks
    plt.text(-14, 11, gene_name, fontdict={'size': 70})
    figure = plt.gcf()
    figure.set_size_inches(12, 12)
    plt.savefig(path + gene_name + '.png', dpi=150, bbox_inches='tight')


def seaborn_plot_overlap(gene_list, filename='kpmp_biopsy_data.feather', path='pics\\', outname=''):
    plt.clf()
    sns.set(rc={'figure.figsize': (12, 12)})
    sns.set_style("white")
    sns.set_context("poster")
    start = time.time()
    gene_list=[i for i in gene_list if i in exist_gene]
    column_list = gene_list+['umap1','umap2']
    df_data = pd.read_feather(filename, columns=column_list)
    df_data['combined'] = df_data[gene_list].agg(sum, axis=1)
    df_data = df_data.sort_values(by='combined', ascending=True)
    print(time.time() - start, gene_list);
    start = time.time()
    ax = sns.scatterplot(data=df_data, x='umap1', y='umap2', hue='combined', s=3, legend=False, palette='magma', x_jitter=True, y_jitter=True, alpha=0.5)
    ax.set(yticklabels=[], xticklabels=[], xlabel='', ylabel='')  # remove the tick labels
    ax.tick_params(left=False, bottom=False)  # remove the ticks
    #plt.text(-16, 13, outname, fontdict={'size': 70})
    figure = plt.gcf()
    figure.set_size_inches(12, 12)
    if outname=='':
        outname = '-'.join(gene_list)
    plt.savefig(path + outname.replace(' ','_').replace('/','') + '.png', dpi=150, bbox_inches='tight')
    #plt.show()


def plot_annotation(annotation_dict, filename='kpmp_biopsy_data.feather', path='pics\\'):
    for each in annotation_dict:
        gene_list = annotation_dict[each]
        seaborn_plot_overlap(gene_list, outname = each, path='pic3\\new\\')


if __name__ == '__main__':
    import h5py
    from data_store import matrisome_dict, annotation_l1

    # filename = r"C:\Users\bathy\Downloads\a87273b3-ec0d-4419-9253-2dc1dcf4a099_WashU-UCSD_KPMP-Biopsy_10X-R_05142021.h5Seurat"
    # fo = h5py.File(filename, 'r')
    # matrisome_list = to_list("KRT13,EEF1A1,S100A6,EHF,AMBRA1,LCN2,SLPI,TACSTD2,RPL41,RPLP1,KCNT2,ALDH1A2,CFH,PAWR,FRMD4A,SLC4A11,LINC01435,PDE1A,RHEX,SYNE1")
    # for key in matrisome_dict:
    #     matrisome_list = matrisome_dict[key]
    #     gene_list = [i.decode('utf8') for i in fo['assays']['RNA']['scaled.features']]
    #     gene_list = [i for i in gene_list if i in matrisome_list]
    #     pout = 'C:\\Users\\bathy\\PycharmProjects\\video_detection\\pic3\\%s\\' % key
    #     if not os.path.isdir(pout):
    #         os.mkdir(pout)
    #     for each in gene_list:
    #         fig_name = pout + each + '.png'
    #         if not os.path.exists(fig_name):
    #             seaborn_plot(gene_name=each, path=pout)
        # seaborn_plot('WISP2')
    # integrin=['ITGA1', 'ITGA10', 'ITGA11', 'ITGA2', 'ITGA2B', 'ITGA3', 'ITGA4', 'ITGA5', 'ITGA6', 'ITGA7', 'ITGA8', 'ITGA9', 'ITGAD', 'ITGAE', 'ITGAL', 'ITGAM', 'ITGAV', 'ITGAX', 'ITGB1', 'ITGB1-DT', 'ITGB1BP1', 'ITGB1BP2', 'ITGB2', 'ITGB2-AS1', 'ITGB3', 'ITGB3BP', 'ITGB4', 'ITGB5', 'ITGB5-AS1', 'ITGB6', 'ITGB7', 'ITGB8', 'ITGBL1']
    # actin=['ACTA1', 'ACTA2', 'ACTB', 'ACTBL2', 'ACTC1', 'ACTG1', 'ACTG2', 'ACTL10', 'ACTL6A', 'ACTL6B', 'ACTL7A', 'ACTL7B', 'ACTL8', 'ACTN1', 'ACTN1-AS1', 'ACTN2', 'ACTN3', 'ACTN4', 'ACTR10', 'ACTR1A', 'ACTR1B', 'ACTR2', 'ACTR3', 'ACTR3B', 'ACTR3C', 'ACTR5', 'ACTR6', 'ACTR8', 'ACTRT3']
    # collagen = matrisome_dict['Collagen']
    # for each in integrin+actin:
    #     for each2 in collagen:
    #         seaborn_plot_overlap([each, each2], path='pic3\\overlap\\')

    # plot_annotation(annotation_l1)

    seaborn_plot_overlap(['COL4A5','COL4A6'],path='pic3\\new\\')