#!/usr/bin/env python


# ----------------------------------------------------------------------------
# Copyright (c) 2016--, tax-credit development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import seaborn as sns
from seaborn import violinplot, heatmap
from pylab import scatter, xlabel, ylabel, xlim, ylim
import matplotlib.pyplot as plt
from scipy.stats import kruskal, linregress


def lmplot_from_data_frame(df, x, y, group_by, style_theme="whitegrid",
                           regress=False):
    '''Make seaborn lmplot from pandas dataframe.
    df: pandas dataframe
    x: str
        x axis variable
    y: str
        y axis variable
    group_by: str
        df variable to use for separating plot panels with FacetGrid
    style_theme: str
        seaborn plot style theme
    '''
    sns.set_style(style_theme)
    sns.lmplot(x, y, col=group_by, data=df, ci=None, size=5,
               scatter_kws={"s": 50, "alpha": 1}, sharey=True)
    sns.plt.show()

    if regress is True:
        return calculate_linear_regress(df, x, y, group_by)


def pointplot_from_data_frame(df, x_axis, y_vars, group_by, color_by,
                              color_pallette, style_theme="whitegrid",
                              plot_type=sns.pointplot):
    '''Generate seaborn pointplot from pandas dataframe.
    df = pandas dataframe
    x_axis = x axis variable
    y_vars = LIST of variables to use for plotting y axis
    group_by = df variable to use for separating plot panels with FacetGrid
    color_by = df variable on which to plot and color subgroups within data
    color_pallette = color palette to use for plotting. Either a dict mapping
                     color_by groups to colors, or a named seaborn palette.
    style_theme = seaborn plot style theme
    plot_type = allows switching to other plot types, but this is untested
    '''
    sns.set_style(style_theme)
    for y_var in y_vars:
        grid = sns.FacetGrid(df, col=group_by, hue=color_by,
                             palette=color_pallette)
        grid = grid.map(sns.pointplot, x_axis, y_var, marker="o", ms=4)
    sns.plt.show()


def heatmap_from_data_frame(df, metric, rows=["Method", "Parameters"],
                            cols=["Dataset"], vmin=0, vmax=1, cmap='Reds'):
    """Generate heatmap of specified metric by (method, parameter) x dataset

    df: pandas.DataFrame
    rows: list
        df column names to use for categorizing heatmap rows
    cols: list
        df column names to use for categorizing heatmap rows
    metric: str
        metric to plot in the heatmap

    """
    df = df.pivot_table(index=rows, columns=cols, values=metric)
    df.sort_index()

    height = len(df.index) * 0.35
    width = len(df.columns) * 1

    # Based on SO answer: http://stackoverflow.com/a/18238680
    ax = plt.figure(figsize=(width, height))
    ax = heatmap(df, cmap=cmap, linewidths=0, square=True, vmin=vmin,
                 vmax=vmax)

    ax.set_title(metric, fontsize=20)

    plt.show()


def boxplot_from_data_frame(df,
                            group_by="Method",
                            metric="Precision",
                            hue=None,
                            y_min=0.0,
                            y_max=1.0,
                            plotf=violinplot,
                            color='grey',
                            x_tick_label_rotation=45):
    """Generate boxplot or violinplot of metric by group

    To generate boxplots instead of violin plots, pass plotf=seaborn.boxplot

    hue, color variables all pass directly to equivalently named
        variables in seaborn.violinplot().

    group_by = "x"
    metric = "y"
    """

    x_tick_labels = df[group_by].unique()
    x_tick_labels.sort()

    ax = violinplot(x=group_by, y=metric, hue=hue, data=df, color=color)
    ax.set_ylim(bottom=y_min, top=y_max)
    ax.set_ylabel(metric)
    ax.set_xlabel(group_by)
    ax.set_xticklabels(x_tick_labels, rotation=x_tick_label_rotation)
    ax


def generate_pr_scatter_plots(query_prf,
                              subject_prf,
                              query_color="b",
                              subject_color="r",
                              x_label="Precision",
                              y_label="Recall"):
    """ Generate scatter plot of precision versus recall for query and subject
    results

        query_prf : pandas.DataFrame
         Precision, recall, and f-measure values as returned from compute_prfs
         for query data
        subject_prf : pandas.DataFrame
         Precision, recall, and f-measure values as returned from compute_prfs
         for subject data
        query_color : str, optional
         The color of the query points
        subject_color : str, optional
         The color of the subject points
        x_label : str, optional
         x axis label for the plot
        y_label : str, optional
         y axis label for the plot

    """
    # Extract the query precisions and recalls and
    # generate a scatter plot
    query_precisions = query_prf['Precision']
    query_recalls = query_prf['Recall']
    scatter(query_precisions,
            query_recalls,
            c=query_color)

    # Extract the subject precisions and recalls and
    # generate a scatter plot
    subject_precisions = subject_prf['Precision']
    subject_recalls = subject_prf['Recall']
    scatter(subject_precisions,
            subject_recalls,
            c=subject_color)

    xlim(0, 1)
    ylim(0, 1)
    xlabel(x_label)
    ylabel(y_label)


def calculate_linear_regress(df, x, y, group_by):
    '''Calculate slope, intercept from series of lines
    df: pandas dataframe
    x: str
        x axis variable
    y: str
        y axis variable
    group_by: str
        df variable to use for separating data subsets
    '''
    results = []
    for group in df[group_by].unique():
        df_mod = df[df[group_by] == group]
        slope, intercept, r_value, p_value, std_err = linregress(df_mod[x],
                                                                 df_mod[y])
        results.append((group, slope, intercept, r_value, p_value, std_err))
    result = pd.DataFrame(results, columns=[group_by, "Slope", "Intercept",
                                            "R", "P-val", "Std Error"])
    return result
