import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load data
# headgroup
headgroup_datafiles = ['s1_contacts_headgroup_frq.tsv', 's2_contacts_headgroup_frq.tsv', 's3_contacts_headgroup_frq.tsv']
glycerol_datafiles = ['s1_contacts_glycerolBB_frq.tsv', 's2_contacts_glycerolBB_frq.tsv', 's3_contacts_glycerolBB_frq.tsv']
tail_datafiles = ['s1_contacts_tail_frq.tsv', 's2_contacts_tail_frq.tsv', 's3_contacts_tail_frq.tsv']

# Read to dataframes
def sort_dfs(datafiles):
    dfs = []
    for i in datafiles:
        df = pd.read_csv(i, comment="#", header=None, names=['Res1', 'Res2', 'Freq'], delim_whitespace=True)
        df['new_res1'] = np.array([j[2:] for j in df['Res1']])
        dfs.append(df)
    
    # Find overlapping records in the first column
    overlapping_records = set(dfs[0]['new_res1']).intersection(*[set(df['new_res1']) for df in dfs[1:]])

    # Filter DataFrames to include only overlapping records
    dfs_filtered = [df[df['new_res1'].isin(overlapping_records)] for df in dfs]
    
    final_data = {'new_res1': [], 'Mean': [], 'Std': []}
    # Iterate over overlapping records and calculate mean and std
    for record in overlapping_records:
        record_values = []
        for df in dfs_filtered:
            record_values.extend(df.loc[df['new_res1'] == record, 'Freq'].values)
        # print(len(record_values))
        if len(record_values) == 12:
            record_mean = sum(record_values) / len(record_values)
            record_std = (sum((x - record_mean) ** 2 for x in record_values) / len(record_values)) ** 0.5 / 12**0.5
            final_data['new_res1'].append(record)
            final_data['Mean'].append(record_mean)
            final_data['Std'].append(record_std)
    final_df = pd.DataFrame(final_data).sort_values(by='Mean', ascending=False).reset_index(drop=True)



    # merged_df = pd.concat(dfs, ignore_index=True)
    # grouped = merged_df.groupby('new_res1')['Freq'].mean().reset_index()
    # grouped_sorted = grouped.sort_values(by='Freq', ascending=False).reset_index(drop=True)
    # grouped_sorted = grouped_sorted[grouped_sorted['Freq']>=0.0]
    # return grouped_sorted
    return final_df

amino_acid_mapping = {
    'ALA': 'A',  # Alanine
    'ARG': 'R',  # Arginine
    'ASN': 'N',  # Asparagine
    'ASP': 'D',  # Aspartic Acid
    'CYS': 'C',  # Cysteine
    'GLN': 'Q',  # Glutamine
    'GLU': 'E',  # Glutamic Acid
    'GLY': 'G',  # Glycine
    'HIS': 'H',  # Histidine
    'HSD': 'H',  # Histidine
    'ILE': 'I',  # Isoleucine
    'LEU': 'L',  # Leucine
    'LYS': 'K',  # Lysine
    'MET': 'M',  # Methionine
    'PHE': 'F',  # Phenylalanine
    'PRO': 'P',  # Proline
    'SER': 'S',  # Serine
    'THR': 'T',  # Threonine
    'TRP': 'W',  # Tryptophan
    'TYR': 'Y',  # Tyrosine
    'VAL': 'V'   # Valine
}

# Set color thresholds
def get_colors(df):
    cutoffs = [0.6, 0.3]
    cmap = plt.get_cmap('RdYlGn')
    colors = [cmap((v - min(df['Mean'])) / (max(df['Mean']) - min(df['Mean']))) for v in df['Mean']]
    for i, v in enumerate(df['Mean']):
        if v >= cutoffs[0]:
            colors[i] = '#D22027' # '#fb8d62'
        elif v >= cutoffs[1]:
            colors[i] = '#385989' # '#8da0cb'
        else:
            colors[i] = '#7FA5B7' # '#66c2a5'
    return colors

# Calculate and plot
headgroup_df = sort_dfs(headgroup_datafiles)
glycerol_df = sort_dfs(glycerol_datafiles)
tail_df = sort_dfs(tail_datafiles)

fig, axs = plt.subplots(1, 3,figsize=(16, 5.33), sharey=True)
axs = axs.ravel()

cutoff_frq_list = [0.2, 0.2, 0.2]

for df,cutoff_frq,ax in zip([headgroup_df, glycerol_df, tail_df], cutoff_frq_list, axs):
    df_screened = df[df['Mean']>cutoff_frq]
    bar_width = 0.075 * df_screened.shape[0]
    df_screened_copy = df_screened.copy()
    print(df_screened_copy)
    colors = get_colors(df_screened_copy)
    
    res1 = [amino_acid_mapping[i[:3]]+i[4:] for i in df_screened_copy['new_res1']]
    df_screened_copy['res1'] = res1
    
    df_screened_copy.plot(kind='bar', x='res1', y='Mean', yerr='Std', ax=ax, color=colors, legend=False, rot=0, width=bar_width, capsize=5)
    ax.set_xlabel('Residue Index', fontsize=12, fontweight='bold')
    ax.set_ylabel('Contact frequency', fontsize=12, fontweight='bold')
    
    ax.tick_params(axis="both", labelsize=12)
    
    # for p in ax.patches:
        # ax.annotate('{:.2f}'.format(p.get_height()), (p.get_x() + p.get_width() / 2., p.get_height()+0.005), \
        #             ha='center', va='center', xytext=(0, 5), textcoords='offset points', fontsize=10) #, fontweight='bold')
    ax.grid(axis='y', ls='--', alpha=0.7)
    ax.set_ylim([0, 1])
plt.tight_layout()
plt.savefig("pose1_contacts.svg")
plt.show()



