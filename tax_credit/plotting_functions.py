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
import numpy as np
from seaborn import violinplot, heatmap
import matplotlib.pyplot as plt
from scipy.stats import (kruskal,
                         linregress,
                         mannwhitneyu,
                         wilcoxon,
                         ttest_ind,
                         ttest_rel)
from statsmodels.sandbox.stats.multicomp import multipletests
from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa
from skbio.stats.distance import anosim
from biom import load_table
from glob import glob
from os.path import join, split
from itertools import combinations
from IPython.display import display, Markdown
from bokeh.plotting import figure, show, output_file
from bokeh.models import (HoverTool,
                          WheelZoomTool,
                          PanTool,
                          ResetTool,
                          SaveTool,
                          ColumnDataSource)
from bokeh.io import output_notebook


def lmplot_from_data_frame(df, x, y, group_by=None, style_theme="whitegrid",
                           regress=False, hue=None, color_palette=None):
    '''Make seaborn lmplot from pandas dataframe.
    df: pandas.DataFrame
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
    lm = sns.lmplot(x, y, col=group_by, data=df, ci=None, size=5,
                    scatter_kws={"s": 50, "alpha": 1}, sharey=True, hue=hue,
                    palette=color_palette)
    sns.plt.show()

    if regress is True:
        try:
            reg = calculate_linear_regress(df, x, y, group_by)
        except ValueError:
            reg = calculate_linear_regress(df, x, y, hue)
    else:
        reg = None

    return lm, reg


def pointplot_from_data_frame(df, x_axis, y_vars, group_by, color_by,
                              color_palette, style_theme="whitegrid",
                              plot_type=sns.pointplot):
    '''Generate seaborn pointplot from pandas dataframe.
    df = pandas.DataFrame
    x_axis = x axis variable
    y_vars = LIST of variables to use for plotting y axis
    group_by = df variable to use for separating plot panels with FacetGrid
    color_by = df variable on which to plot and color subgroups within data
    color_palette = color palette to use for plotting. Either a dict mapping
                     color_by groups to colors, or a named seaborn palette.
    style_theme = seaborn plot style theme
    plot_type = allows switching to other plot types, but this is untested
    '''
    grid = dict()
    sns.set_style(style_theme)
    for y_var in y_vars:
        grid[y_var] = sns.FacetGrid(df, col=group_by, hue=color_by,
                                    palette=color_palette)
        grid[y_var] = grid[y_var].map(
            sns.pointplot, x_axis, y_var, marker="o", ms=4)
    sns.plt.show()
    return grid


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

    ax = plt.figure(figsize=(width, height))
    ax = heatmap(df, cmap=cmap, linewidths=0, square=True, vmin=vmin,
                 vmax=vmax)

    ax.set_title(metric, fontsize=20)

    plt.show()

    return ax


def boxplot_from_data_frame(df,
                            group_by="Method",
                            metric="Precision",
                            hue=None,
                            y_min=0.0,
                            y_max=1.0,
                            plotf=violinplot,
                            color='grey',
                            color_palette=None,
                            label_rotation=45):
    """Generate boxplot or violinplot of metric by group

    To generate boxplots instead of violin plots, pass plotf=seaborn.boxplot

    hue, color variables all pass directly to equivalently named
        variables in seaborn.violinplot().

    group_by = "x"
    metric = "y"
    """

    sns.set_style("whitegrid")
    ax = violinplot(x=group_by, y=metric, hue=hue, data=df, color=color,
                    palette=color_palette, order=sorted(df[group_by].unique()))
    ax.set_ylim(bottom=y_min, top=y_max)
    ax.set_ylabel(metric)
    ax.set_xlabel(group_by)
    for lab in ax.get_xticklabels():
        lab.set_rotation(label_rotation)

    plt.show()

    return ax


def calculate_linear_regress(df, x, y, group_by):
    '''Calculate slope, intercept from series of lines
    df: pandas.DataFrame
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


def per_level_kruskal_wallis(df,
                             y_vars,
                             group_by,
                             dataset_col='Dataset',
                             level_name="level",
                             levelrange=range(1, 7),
                             alpha=0.05,
                             pval_correction='fdr_bh'):

    '''Test whether 2+ population medians are different.

    Due to the assumption that H has a chi square distribution, the number of
    samples in each group must not be too small. A typical rule is that each
    sample must have at least 5 measurements.

    df = pandas.DataFrame
    y_vars = LIST of variables (df column names) to test
    group_by = df variable to use for separating subgroups to compare
    dataset_col = df variable to use for separating individual datasets to test
    level_name = df variable name that specifies taxonomic level
    levelrange = range of taxonomic levels to test.
    alpha = level of alpha significance for test
    pval_correction = type of p-value correction to use
    '''
    dataset_list = []
    p_list = []
    for dataset in df[dataset_col].unique():
        df1 = df[df[dataset_col] == dataset]
        for var in y_vars:
            dataset_list.append((dataset, var))
            for level in levelrange:
                level_subset = df1[level_name] == level

                # group data by groups
                group_list = []
                for group in df1[group_by].unique():
                    group_data = df1[group_by] == group
                    group_results = df1[level_subset & group_data][var]
                    group_list.append(group_results)

                # kruskal-wallis tests
                try:
                    h_stat, p_val = kruskal(*group_list, nan_policy='omit')
                # default to p=1.0 if all values = 0
                # this is not technically correct, from the standpoint of p-val
                # correction below makes p-vals very slightly less significant
                # than they should be
                except ValueError:
                    h_stat, p_val = ('na', 1)  # noqa

                p_list.append(p_val)

    # correct p-values
    rej, pval_corr, alphas, alphab = multipletests(np.array(p_list),
                                                   alpha=alpha,
                                                   method=pval_correction)

    range_len = len([i for i in levelrange])
    results = [(dataset_list[i][0], dataset_list[i][1],
                *[pval_corr[i * range_len + n] for n in range(0, range_len)])
               for i in range(0, len(dataset_list))]
    result = pd.DataFrame(results, columns=[dataset_col, "Variable",
                                            *[n for n in levelrange]])
    return result


def seek_tables(expected_results_dir, table_fn='merged_table.biom'):
    '''Find and deliver merged biom tables'''
    table_fps = glob(join(expected_results_dir, '*', '*', table_fn))
    for table in table_fps:
        reference_dir, _ = split(table)
        dataset_dir, reference_id = split(reference_dir)
        _, dataset_id = split(dataset_dir)
        yield table, dataset_id, reference_id


def batch_beta_diversity(expected_results_dir, method="braycurtis",
                         permutations=99, col='method', dim=2,
                         colormap={'expected': 'red', 'rdp': 'seagreen',
                                   'sortmerna': 'gray', 'uclust': 'blue',
                                   'blast': 'purple'}):

    '''Find merged biom tables and run beta_diversity_through_plots'''
    for table, dataset_id, reference_id in seek_tables(expected_results_dir):
        display(Markdown('## {0} {1}'.format(dataset_id, reference_id)))
        s, r, pc, dm = beta_diversity_pcoa(table, method=method, col=col,
                                           permutations=permutations, dim=dim,
                                           colormap=colormap)
        sns.plt.show()
        sns.plt.clf()


def make_distance_matrix(biom_fp, method="braycurtis"):
    '''biom.Table --> skbio.DistanceMatrix'''
    table = load_table(biom_fp)

    # extract sample metadata from table, put in df
    table_md = {s_id: dict(table.metadata(s_id)) for s_id in table.ids()}
    s_md = pd.DataFrame.from_dict(table_md, orient='index')

    # extract data from table and multiply, assuming that table contains
    # relative abundances (which cause beta_diversity to fail)
    table_data = [[int(num * 100000) for num in table.data(s_id)]
                  for s_id in table.ids()]

    # beta diversity
    dm = beta_diversity(method, table_data, table.ids())

    return dm, s_md


def beta_diversity_pcoa(biom_fp, method="braycurtis", permutations=99, dim=2,
                        col='method', colormap={'expected': 'red',
                                                'rdp': 'seagreen',
                                                'sortmerna': 'gray',
                                                'uclust': 'blue',
                                                'blast': 'purple'}):

    '''From biom table, compute Bray-Curtis distance; generate PCoA plot;
    and calculate adonis differences.

    biom_fp: path
        Path to biom.Table containing sample metadata.
    method: str
        skbio.Diversity method to use for ordination.
    permutations: int
        Number of permutations to perform for anosim tests.
    dim: int
        Number of dimensions to plot. Currently supports only 2-3 dimensions.
    col: str
        metadata name to use for distinguishing groups for anosim tests and
        pcoa plots.
    colormap: dict
        map groups names (must be group names in col) to colors used for plots.
    '''

    dm, s_md = make_distance_matrix(biom_fp, method=method)

    # pcoa
    pc = pcoa(dm)

    # anosim tests
    results = anosim(dm, s_md, column=col, permutations=permutations)
    print('R = ', results['test statistic'], '; P = ', results['p-value'])

    if dim == 2:
        # bokeh pcoa plots
        pc123 = pc.samples.ix[:, ["PC1", "PC2", "PC3"]]
        smd_merge = s_md.merge(pc123, left_index=True, right_index=True)
        smd_merge['Color'] = [colormap[x] for x in smd_merge['method']]
        title = smd_merge['reference'][0]
        labels = ['PC {0} ({1:.2f})'.format(d + 1, pc.proportion_explained[d])
                  for d in range(0, 2)]
        circle_plot_from_dataframe(smd_merge, "PC1", "PC2", title,
                                   columns=["method", "sample_id", "params"],
                                   color="Color", labels=labels)

    else:
        # skbio pcoa plots
        pcoa_plot_skbio(pc, s_md, col='method')

    return s_md, results, pc, dm


def circle_plot_from_dataframe(df, x, y, title=None, color="Color",
                               columns=["method", "sample_id", "params"],
                               labels=None, plot_width=400, plot_height=400,
                               fill_alpha=0.2, size=10, output_fn=None):
    '''Make bokeh circle plot from dataframe, use df columns for hover tool.
    df: pandas.DataFrame
        Containing all sample data, including color categories.
    x: str
        df category to use for x-axis coordinates.
    y: str
        df category to use for y-axis coordinates.
    title: str
        Title to print above plot.
    color: str
        df category to use for coloring data points.
    columns: list
        df categories to add as hovertool metadata.
    labels: list
        Axis labels for x and y axes. If none, default to column names.
    output_fn: path
        Filepath for output file. Defaults to None.

    Other parameters feed directly to bokeh.plotting.
    '''
    if labels is None:
        labels = [x, y]

    source = ColumnDataSource(df)
    hover = HoverTool(tooltips=[(c, '@' + c) for c in columns])

    TOOLS = [hover, WheelZoomTool(), PanTool(), ResetTool(), SaveTool()]

    fig = figure(title=title, tools=TOOLS, plot_width=plot_width,
                 plot_height=plot_height)

    # Set asix labels
    fig.xaxis.axis_label = labels[0]
    fig.yaxis.axis_label = labels[1]

    # Plot x and y axes
    fig.circle(x, y, source=source, color=color, fill_alpha=fill_alpha,
               size=size)

    if output_fn is not None:
        output_file(output_fn)

    output_notebook()
    show(fig)


def pcoa_plot_skbio(pc, s_md, col='method'):
    '''Input principal coordinates, display figure.

    pc: skbio.OrdinationResults
        Sample coordinates.
    s_md: pandas.DataFrame
        Sample metadata.
    col: str
        Category in s_md to use for coloring groups.
    '''

    # make labels for PCoA plot
    pcl = ['PC {0} ({1:.2f})'.format(d + 1, pc.proportion_explained[d])
           for d in range(0, 3)]
    fig = pc.plot(s_md, col, axis_labels=(pcl[0], pcl[1], pcl[2]),
                  cmap='jet', s=50)
    fig


def average_distance_boxplots(expected_results_dir, group_by="method",
                              standard='expected', metric="distance",
                              params='params', beta="braycurtis",
                              reference_filter=True, reference_col='reference',
                              references=['gg_13_8_otus',
                                          'unite_20.11.2016_clean_fullITS'],
                              paired=True, use_best=True, parametric=True,
                              plotf=violinplot, label_rotation=45,
                              color_palette=None, y_min=0.0, y_max=1.0,
                              color=None, hue=None):

    '''Distance boxplots that aggregate and average results across multiple
    mock community datasets.

    reference_filter: bool
        Filter by reference dataset to only include specific references in
        results?
    reference_col: str
        df column header containing reference set information.
    references: list
        List of strings containing names of reference datasets to include.
    paired: bool
        Perform paired or unpaired comparisons?
    parametric: bool
        Perform parametric or non-parametric statistical tests?
    use_best: bool
        Compare average distance distributions across all methods (False) or
        only the best parameter configuration for each method? (True)
    '''
    box = dict()
    best = dict()
    # Aggregate all distance matrix data
    archive = pd.DataFrame()
    for table, dataset_id, reference_id in seek_tables(expected_results_dir):
        dm, sample_md = make_distance_matrix(table, method=beta)
        per_method = per_method_distance(dm, sample_md, group_by=group_by,
                                         standard=standard, metric=metric)
        archive = pd.concat([archive, per_method])

    # filter out auxiliary reference database results
    if reference_filter is True:
        archive = archive[archive[reference_col].isin(references)]

    # plot results for each reference db separately
    for reference in archive[reference_col].unique():
        display(Markdown('## {0}'.format(reference)))
        archive_subset = archive[archive[reference_col] == reference]

        # for each method find best average method/parameter config
        if use_best:
            best[reference], param_report = isolate_top_params(
                archive_subset, group_by, params, metric)
            # display(pd.DataFrame(param_report, columns=[group_by, params]))

            method_rank = _show_method_rank(
                best[reference], group_by, params, metric,
                [group_by, params, metric], ascending=False)

        else:
            best[reference] = archive_subset

        results = per_method_pairwise_tests(best[reference], group_by=group_by,
                                            metric=metric, paired=paired,
                                            parametric=parametric)

        box[reference] = boxplot_from_data_frame(
            best[reference], group_by=group_by, color=color, hue=hue,
            y_min=None, y_max=None, plotf=plotf, label_rotation=label_rotation,
            metric=metric, color_palette=color_palette)

        if use_best:
            box[reference] = _add_significance_to_boxplots(
                results, method_rank, box[reference], method='method')

        sns.plt.show()
        sns.plt.clf()

        display(results)

    return box, best


def _add_significance_to_boxplots(pairwise, rankings, ax, method='Method'):

    x_labels = [a.get_text() for a in ax.get_xticklabels()]
    methods = [m for m in rankings[method]]

    ranks = []
    # iterate range instead of methods, so that we can use pop for comparisons
    # against shrinking list of methods
    methods_copy = methods.copy()
    for n in range(len(methods_copy)):
        method = methods_copy.pop(0)
        inner_rank = {method}
        for other_method in methods_copy:
            if method in pairwise.index.levels[0] and \
                    other_method in pairwise.loc[method].index:
                if pairwise.loc[method].loc[other_method]['FDR P'] > 0.05:
                    inner_rank.add(other_method)
            elif other_method in pairwise.index.levels[0] and \
                    method in pairwise.loc[other_method].index:
                if pairwise.loc[other_method].loc[method]['FDR P'] > 0.05:
                    inner_rank.add(other_method)
        # only add new set of equalities if it contains unique items
        if len(ranks) == 0 or not inner_rank.issubset(ranks[-1]):
            ranks.append(inner_rank)

    # provide unique letters for each significance group
    letters = 'abcdefghijklmnopqrstuvwxyz'

    sig_groups = {}
    for method in methods:
        sig_groups[method] = []
        for rank, letter in zip(ranks, letters):
            if method in rank:
                sig_groups[method].append(letter)
        sig_groups[method] = ''.join(sig_groups[method])

    # add significance labels above plot
    pos = range(len(x_labels))
    for tick, label in zip(pos, x_labels):
        ax.text(tick, ax.get_ybound()[1], sig_groups[label], size='medium',
                horizontalalignment='center', color='k', weight='semibold')

    return ax


def _show_method_rank(best, group_by, params, metric, display_fields,
                      ascending=False):
    '''Find the best param configuration for each method and show those
    configs, along with the parameters and metric scores.
    '''
    avg_best = best.groupby([group_by, params]).mean().reset_index()
    avg_best_sorted = avg_best.sort_values(by=metric, ascending=ascending)
    method_rank = avg_best_sorted.ix[:, display_fields]
    display(method_rank)
    return method_rank


def fastlane_boxplots(expected_results_dir, group_by="method",
                      standard='expected', metric="distance", hue=None,
                      plotf=violinplot, label_rotation=45,
                      y_min=0.0, y_max=1.0, color=None, beta="braycurtis"):

    '''per_method_boxplots for those who don't have time to wait.'''

    for table, dataset_id, reference_id in seek_tables(expected_results_dir):
        display(Markdown('## {0} {1}'.format(dataset_id, reference_id)))

        dm, sample_md = make_distance_matrix(table, method=beta)

        per_method_boxplots(dm, sample_md, group_by=group_by, metric=metric,
                            standard=standard, hue=hue, y_min=y_min,
                            y_max=y_max, plotf=plotf, color=color,
                            label_rotation=label_rotation)


def per_method_boxplots(dm, sample_md, group_by="method", standard='expected',
                        metric="distance", hue=None, y_min=0.0, y_max=1.0,
                        plotf=violinplot, label_rotation=45, color=None,
                        color_palette=None):
    '''Generate distance boxplots and Mann-Whitney U tests on distance matrix.

    dm: skbio.DistanceMatrix
    sample_md: pandas.DataFrame
        containing sample metadata
    group_by: str
        df category to use for grouping samples
    standard: str
        group name in group_by category to which all other groups are compared.
    metric: str
        name of distance column in output.

    To generate boxplots instead of violin plots, pass plotf=seaborn.boxplot

    hue, color variables all pass directly to equivalently named variables in
        seaborn.violinplot().
    '''
    box = dict()
    within_between = within_between_category_distance(dm, sample_md, 'method')

    per_method = per_method_distance(dm, sample_md, group_by=group_by,
                                     standard=standard, metric=metric)

    for d, g, s in [(within_between, 'Comparison', '1: Within- vs. Between-'),
                    (per_method, group_by, '2: Pairwise ')]:

        display(Markdown('## Comparison {0} Distance'.format(s + group_by)))
        box[g] = boxplot_from_data_frame(
            d, group_by=g, color=color, metric=metric, y_min=None, y_max=None,
            hue=hue, plotf=plotf, label_rotation=label_rotation,
            color_palette=color_palette)

        results = per_method_pairwise_tests(d, group_by=g, metric=metric)

        sns.plt.show()
        sns.plt.clf()
        display(results)

    return box


def per_method_distance(dm, md, group_by='method', standard='expected',
                        metric='distance', sample='sample_id'):
    '''Compile list of distances between groups of samples in distance matrix.
    returns dataframe of distances and group metadata.

    dm: skbio.DistanceMatrix
    md: pandas.DataFrame
        containing sample metadata
    group_by: str
        df category to use for grouping samples
    standard: str
        group name in group_by category to which all other groups are compared.
    metric: str
        name of distance column in output.
    sample: str
        df category containing sample_id names.
    '''
    results = []
    expected = md[md[group_by] == standard]
    observed = md[md[group_by] != standard]
    for group in observed[group_by].unique():
        group_md = observed[observed[group_by] == group]
        for i in list(expected.index.values):
            for j in list(group_md.index.values):
                if group_md.loc[j][sample] == expected.loc[i][sample]:
                    results.append((*[n for n in group_md.loc[j]], dm[i, j]))
    return pd.DataFrame(results, columns=[*[n for n in md.columns.values],
                                          metric])


def within_between_category_distance(dm, md, md_category, distance='distance'):
    '''Compile list of distances between groups of samples and within groups
    of samples.

    dm: skbio.DistanceMatrix
    md: pandas.DataFrame
        containing sample metadata
    md_category: str
        df category to use for grouping samples
    '''
    distances = []
    for i, sample_id1 in enumerate(dm.ids):
        sample_md1 = md[md_category][sample_id1]
        for sample_id2 in dm.ids[:i]:
            sample_md2 = md[md_category][sample_id2]
            if sample_md1 == sample_md2:
                comp = 'within'
                group = sample_md1
            else:
                comp = 'between'
                group = sample_md1 + '_' + sample_md2
            distances.append((comp, group, dm[sample_id1, sample_id2]))
    return pd.DataFrame(distances, columns=["Comparison", md_category,
                                            distance])


def per_method_pairwise_tests(df, group_by='method', metric='distance',
                              paired=False, parametric=True):
    '''Perform mann whitney U tests between group distance distributions,
    followed by FDR correction. Returns pandas dataframe of p-values.
    df: pandas.DataFrame
        results from per_method_distance()
    group_by: str
        df category to use for grouping samples
    metric: str
        df category to use as variable for comparison.
    paired: bool
        Perform Wilcoxon signed rank test instead of Mann Whitney U. df must be
        ordered such that paired samples will appear in same order in subset
        dataframes when df is subset by term f[df[group_by] == a[0]][metric].
    '''
    pvals = []
    groups = [group for group in df[group_by].unique()]
    combos = [a for a in combinations(groups, 2)]
    for a in combos:
        try:
            if paired is False and parametric is False:
                u, p = mannwhitneyu(df[df[group_by] == a[0]][metric],
                                    df[df[group_by] == a[1]][metric],
                                    alternative='two-sided')
            elif paired is False and parametric is True:
                u, p = ttest_ind(df[df[group_by] == a[0]][metric],
                                 df[df[group_by] == a[1]][metric],
                                 nan_policy='raise')
            elif paired is True and parametric is False:
                u, p = wilcoxon(df[df[group_by] == a[0]][metric],
                                df[df[group_by] == a[1]][metric])
            else:
                u, p = ttest_rel(df[df[group_by] == a[0]][metric],
                                 df[df[group_by] == a[1]][metric],
                                 nan_policy='raise')
        except ValueError:
            # default to p=1.0 if all values = 0
            # this is not technically correct, from the standpoint of p-val
            # correction below makes p-vals very slightly less significant
            # than they should be
            u, p = 0.0, 1.0
        pvals.append((a[0], a[1], u, p))

    result = pd.DataFrame(pvals, columns=["Method A", "Method B", "stat", "P"])
    result.set_index(['Method A', 'Method B'], inplace=True)
    try:
        result['FDR P'] = multipletests(result['P'], method='fdr_bh')[1]
    except ZeroDivisionError:
        pass

    return result


def isolate_top_params(df, group_by="Method", params="Parameters",
                       metric="F-measure", ascending=True):
    '''For each method in df, find top params for each method and filter df to
    contain only those parameters.

    df: pandas df
    group_by: str
        df category name to use for segregating groups from which top param is
        chosen.
    params: str
        df category name indicating parameters column.
    '''
    best = pd.DataFrame()
    param_report = []
    for group in df[group_by].unique():
        subset = df[df[group_by] == group]
        avg = subset.groupby(params).mean().reset_index()
        sorted_avg = avg.sort_values(by=metric, ascending=ascending)
        top_param = sorted_avg.reset_index()[params][0]
        param_report.append((group, top_param))
        best = pd.concat([best, subset[subset[params] == top_param]])
    return best, param_report


def rank_optimized_method_performance_by_dataset(df,
                                                 dataset="Dataset",
                                                 method="Method",
                                                 params="Parameters",
                                                 metric="F-measure",
                                                 level="Level",
                                                 level_range=range(5, 7),
                                                 display_fields=["Method",
                                                                 "Parameters",
                                                                 "Precision",
                                                                 "Recall",
                                                                 "F-measure"],
                                                 ascending=False,
                                                 paired=True,
                                                 parametric=True,
                                                 hue=None,
                                                 y_min=0.0,
                                                 y_max=1.0,
                                                 plotf=violinplot,
                                                 label_rotation=45,
                                                 color=None,
                                                 color_palette=None):

    '''Rank the performance of methods using optimized parameter configuration
    within each dataset in dataframe. Optimal methods are computed from the
    mean performance of each method/param configuration across all datasets
    in df.

    df: pandas df
    dataset: str
        df category to use for grouping samples (by dataset)
    method: str
        df category to use for grouping samples (by method); these groups are
        compared in plots and pairwise statistical testing.
    params: str
        df category containing parameter configurations for each method. Best
        method configurations are computed by grouping method groups on this
        category value, then finding the best average metric value.
    metric: str
        df category containing metric to use for ranking and statistical
        comparisons between method groups.
    level: str
        df category containing taxonomic level information.
    level_range: range
        Perform plotting and testing at each level in range.
    display_fields: list
        List of columns in df to display in results.
    ascending: bool
        Rank methods my metric score in ascending or descending order?
    paired: bool
        Perform paired statistical test? See per_method_pairwise_tests()
    parametric: bool
        Perform parametric statistical test? See per_method_pairwise_tests()

    To generate boxplots instead of violin plots, pass plotf=seaborn.boxplot

    hue, color variables all pass directly to equivalently named variables in
        seaborn.violinplot(). See boxplot_from_data_frame() for more
        information.
    '''
    box = dict()
    for d in df[dataset].unique():
        for lv in level_range:
            display(Markdown("## {0} level {1}".format(d, lv)))
            df_l = df[df[level] == lv]
            best, param_report = isolate_top_params(df_l[df_l[dataset] == d],
                                                    method, params, metric,
                                                    ascending=ascending)

            method_rank = _show_method_rank(
                best, method, params, metric, display_fields,
                ascending=ascending)

            results = per_method_pairwise_tests(best, group_by=method,
                                                metric=metric, paired=paired,
                                                parametric=parametric)
            display(results)

            box[d] = boxplot_from_data_frame(
                best, group_by=method, color=color, metric=metric, y_min=y_min,
                y_max=y_max, label_rotation=label_rotation, hue=hue,
                plotf=plotf, color_palette=color_palette)

            box[d] = _add_significance_to_boxplots(
                results, method_rank, box[d])

            sns.plt.show()

    return box
