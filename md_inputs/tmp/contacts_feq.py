# Author: Jian Huang
# Date: Aug 11 2023
# E-mail: jianhuang@umass.edu

# TODO:
# base on given contact file, plot each PIP2 contacts separately
#    <-- see s1_all_frq.tsv

import pandas as pd
import matplotlib.pyplot as plt
import sys

contact_datafile = sys.argv[1]
# print(len(sys.argv))
if len(sys.argv) == 3:
    savefig = sys.argv[2]

# load 
df = pd.read_csv(contact_datafile, comment="#", header=None, names=['Res1', 'Res2', 'Freq'], delim_whitespace=True)
df_I = df[df['Res2'].str.startswith('I')]
df_I_sorted = df_I.sort_values(by='Freq', ascending=False, ignore_index=True)
df_J = df[df['Res2'].str.startswith('J')]
df_J_sorted = df_J.sort_values(by='Freq', ascending=False, ignore_index=True)
df_K = df[df['Res2'].str.startswith('K')]
df_K_sorted = df_K.sort_values(by='Freq', ascending=False, ignore_index=True)
df_L = df[df['Res2'].str.startswith('L')]
df_L_sorted = df_L.sort_values(by='Freq', ascending=False, ignore_index=True)


# plot
fig, ax = plt.subplots(2,2,figsize=(12, 8))
# Set color thresholds
def get_colors(df):
    cutoffs = [0.6, 0.3]
    cmap = plt.get_cmap('RdYlGn')
    colors = [cmap((v - min(df['Freq'])) / (max(df['Freq']) - min(df['Freq']))) for v in df['Freq']]
    for i, v in enumerate(df['Freq']):
        if v >= cutoffs[0]:
            colors[i] = 'g'
        elif v >= cutoffs[1]:
            colors[i] = 'orange'
        else:
            colors[i] = 'r'
    return colors

# only plot Freq > 0.1
cutoff_frq = 0.1

df_I_screened = df_I_sorted[df_I_sorted['Freq']>cutoff_frq]
I_colors = get_colors(df_I_screened)
df_I_screened.plot(kind='bar', x='Res1', y='Freq', ax=ax[0][0], color=I_colors, legend=False, rot=45)
ax[0][0].set_xlabel('PIP2(I)', fontsize=10, fontweight='bold')
ax[0][0].set_ylabel('Contact frq', fontsize=14, fontweight='bold')


df_J_screened = df_J_sorted[df_J_sorted['Freq']>cutoff_frq]
J_colors = get_colors(df_J_screened)
df_J_screened.plot(kind='bar', x='Res1', y='Freq', ax=ax[0][1], color=J_colors, legend=False, rot=45)
ax[0][1].set_xlabel('PIP2(J)', fontsize=10, fontweight='bold')
ax[0][1].set_ylabel('Contact frq', fontsize=14, fontweight='bold')

df_K_screened = df_K_sorted[df_K_sorted['Freq']>cutoff_frq]
K_colors = get_colors(df_K_screened)
df_K_screened.plot(kind='bar', x='Res1', y='Freq', ax=ax[1][0], color=K_colors, legend=False, rot=45)
ax[1][0].set_xlabel('PIP2(K)', fontsize=10, fontweight='bold')
ax[1][0].set_ylabel('Contact frq', fontsize=14, fontweight='bold')

df_L_screened = df_L_sorted[df_L_sorted['Freq']>cutoff_frq]
L_colors = get_colors(df_L_screened)
df_L_screened.plot(kind='bar', x='Res1', y='Freq', ax=ax[1][1], color=L_colors, legend=False, rot=45)
ax[1][1].set_xlabel('PIP2(L)', fontsize=10, fontweight='bold')
ax[1][1].set_ylabel('Contact frq', fontsize=14, fontweight='bold')

# annotate each bar
for i in ax.flatten():
    for p in i.patches:
        i.annotate('{:.2f}'.format(p.get_height()), (p.get_x() + p.get_width() / 2., p.get_height()), \
                    ha='center', va='center', xytext=(0, 5), textcoords='offset points', fontsize=5, fontweight='bold')
        i.tick_params(axis='both', which='major', labelsize=14)



plt.tight_layout()
if len(sys.argv) == 3:
    plt.savefig(savefig)
plt.show()

