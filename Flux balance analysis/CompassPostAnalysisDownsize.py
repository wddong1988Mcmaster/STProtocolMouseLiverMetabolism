
#
# #### Different cell types compairism


# def install_reqs():
#   !pip install pandas
#   !pip install "matplotlib>=3.4"
#   !pip install numpy
#   !pip install statsmodels
#   !pip install scipy
# install_reqs()

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


from compass_analysis import cohens_d, wilcoxon_test, get_reaction_consistencies, get_metareactions, labeled_reactions, amino_acid_metab

from matplotlib import __version__ as matplotlibversion
if matplotlibversion < "3.4":
    print("Matplotlib versions older than 3.4 may not be able to generate figure 2E, as they do not support alpha arrays")

import os
os.system("compass --data extdata1small/LivMeta/linear_gene_expression_matrix.tsv \
          --model RECON2_mat --species mus_musculus --media default-media --lambda 0.25 \
          --and-function mean --output-dir extdata1small/LivMeta --penalty-diffusion knn --num-neighbors 10 \
          --isoform-summing legacy --num-processes 50")

reaction_penalties = pd.read_csv("extdata1small/LivMeta/reactions.tsv", sep="\t", index_col = 0)
cell_metadata = pd.read_csv("extdata1small/LivMeta/cell_metadata.csv", index_col = 0)

# Compare different cell types
LivT_cells = cell_metadata.index[cell_metadata['celltype'] == '0: Hepatocytes1']
LivC_cells = cell_metadata.index[cell_metadata['celltype'] == '10: B cells']

# The reaction metadata for RECON2
reaction_metadata = pd.read_csv("extdata1small/RECON2/reaction_metadata.csv", index_col = 0)

# one example of a row for the reaction metadata:
reaction_metadata.loc[['r0281']]

#This function is repeated here for clarity
def get_reaction_consistencies(compass_reaction_penalties, min_range=1e-3):
    """
        Converts the raw penalties outputs of compass into scores per reactions where higher numbers indicate more activity
    """
    df = -np.log(compass_reaction_penalties + 1)
    df = df[df.max(axis=1) - df.min(axis=1) >= min_range]
    df = df - df.min().min()
    return df

reaction_consistencies = get_reaction_consistencies(reaction_penalties)

# We use the unpaired Wilcoxon rank-sum test (equivlanet to the Mann–Whitney U test) to analyze
# pathogenic Th17p cells compared to the non-pathogenic Th17n
wilcox_results = wilcoxon_test(reaction_consistencies, LivT_cells, LivC_cells)
wilcox_results['metadata_r_id'] = ""
for r in wilcox_results.index:
    if r in reaction_metadata.index:
        wilcox_results.loc[r, 'metadata_r_id'] = r
    elif r[:-4] in reaction_metadata.index:
        wilcox_results.loc[r, 'metadata_r_id'] = r[:-4]
    else:
        print("Should not occur")

# we join the metadata to the reactions in a new dataframe W
W = wilcox_results.merge(reaction_metadata, how='left',
                         left_on='metadata_r_id', right_index=True, validate='m:1')
W = W[W['confidence'].isin([0, 4])]
W = W[~W['EC_number'].isna()]
W.loc[(W['formula'].map(lambda x: '[m]' not in x)) & (W['subsystem'] == "Citric acid cycle"), 'subsystem'] = 'Other'

# one example of a row of the reuslting dataframe for a reaction
wilcox_results.loc[['r0281_pos']]
print(wilcox_results.loc[['r0281_pos']])

reaction_metadata.loc['r0281']['formula']
print(reaction_metadata.loc['r0281']['formula'])

def plot_differential_scores(data, title, c):
    plt.figure(figsize=(10,10))
    axs = plt.gca()
    axs.scatter(data['cohens_d'], -np.log10(data['adjusted_pval']), c=c)
    axs.set_xlabel("Cohen's d", fontsize=16)
    axs.set_ylabel("-log10 (Wilcoxon-adjusted p)", fontsize=16)
    #Everything after this could be tweaked depending on the application
    axs.set_xlim(-2.2, 2.2)
    axs.axvline(0, dashes=(3,3), c='black')
    axs.axhline(1, dashes=(3,3), c='black')
    axs.set_title(title, fontdict={'fontsize':20})
    axs.annotate('', xy=(0.5, -0.08), xycoords='axes fraction', xytext=(0, -0.08),
            arrowprops=dict(arrowstyle="<-", color='#348C73', linewidth=4))
    axs.annotate('Hepatocytes', xy=(0.55, -0.12), xycoords='axes fraction', fontsize=16)
    axs.annotate('', xy=(0.5, -0.08), xycoords='axes fraction', xytext=(1, -0.08),
            arrowprops=dict(arrowstyle="<-", color='#E92E87', linewidth=4))
    axs.annotate('B cells', xy=(0.05, -0.12), xycoords='axes fraction', fontsize=16)
    for r in data.index:
        if r in labeled_reactions:
            x = data.loc[r, 'cohens_d']
            y = -np.log10(data.loc[r, 'adjusted_pval'])
            offset = (20, 0)
            if x < 0:
                offset = (-100, -40)
            axs.annotate(labeled_reactions[r], (x,y), xytext = offset,
                         textcoords='offset pixels', arrowprops={'arrowstyle':"-"})


# print(W['subsystem'].unique()) (optional)

filtered_data = pd.concat([W[W['subsystem'] == "Glycolysis/gluconeogenesis"],
             W[W['subsystem'] == "Citric acid cycle"],
            W[W['subsystem'].isin(amino_acid_metab)],
           W[W['subsystem'] == "Fatty acid oxidation"],
         W[W['subsystem'] == "Cholesterol metabolism"],
        W[W['subsystem'] == "Fatty acid synthesis"]])

data = W[W['subsystem'] == "Glycolysis/gluconeogenesis"]
plot_differential_scores(data, title='Glycolysis', c="#695D73")

# save as PNG (or png, jpg, jpeg, pdf, svg)
plt.savefig('Glycolysis.png', dpi=300)

# show figure (optional)
# plt.show()

data = W[W['subsystem'] == "Citric acid cycle"]
plot_differential_scores(data, title="TCA Cycle", c="#D3A991")
plt.savefig('TCA Cycle.png', dpi=300)

data = W[W['subsystem'].isin(amino_acid_metab)].copy()
data['adjusted_pval'] = data['adjusted_pval'].clip(1e-12)
plot_differential_scores(data, "Amino Acid Metabolism", c="#BF1E2E")

data = W[W['subsystem'] == "Fatty acid oxidation"]
plot_differential_scores(data, "Fatty Acid Oxidation", c="#040772")
plt.savefig('FAO.png', dpi=300)

data = W[W['subsystem'] == "Cholesterol metabolism"]
plot_differential_scores(data, "Cholesterol metabolism", c="#040772")
plt.savefig('CholMeta.png', dpi=300)

data = W[W['subsystem'] == "Fatty acid synthesis"]
plot_differential_scores(data, "DNL (Fatty acid synthesis)", c="#040772")
plt.savefig('DNL.png', dpi=300)


# Another FIGURE: show overall changes in all pathways
data = W[~W['subsystem'].isin(["Miscellaneous", "Unassigned"])]
data = data[~data['subsystem'].map(lambda x: "Transport" in x or "Exchange" in x or x == "Other")]
items, counts = np.unique(data['subsystem'], return_counts=True)
items = [items[i] for i in range(len(items)) if counts[i] > 5] #filter(n() > 5) %>%
data = data[data['subsystem'].isin(items)]

plt.figure(figsize=(12, 12))
axs = plt.gca()
#Sorts the reactions for plotting
d = data[data['adjusted_pval'] < 0.1].groupby('subsystem')['cohens_d'].median().abs()
axs.scatter(d[d.argsort], d[d.argsort].index, alpha=0)
color = data['cohens_d'].map(lambda x: 'r' if x >= 0 else 'b')
alpha = data['adjusted_pval'].map(lambda x: 1.0 if x < 0.1 else 0.25)
axs.scatter(data['cohens_d'], data['subsystem'], c=color, alpha=alpha)
axs.set_xlabel("Cohen's d")
axs.set_title("Differential Metabolic Activity in hepatocytes vs B cells", fontsize=14)
plt.tight_layout()
plt.subplots_adjust(bottom=0.2)
plt.show()
plt.savefig('FigA.png', dpi=300, bbox_inches='tight')


