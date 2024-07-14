import pandas as pd
import matplotlib.pyplot as plt
import sys
import numpy as np

# raw contacts data from get_dynamic_contacts.py calculation
contacts_data = sys.argv[1]
# output_prefix = sys.argv[2]

headgroup_contact_output = contacts_data.split('.')[0] + '_headgroup.tsv'
glycerolBB_contact_output = contacts_data.split('.')[0] + '_glycerolBB.tsv'
tail_contact_output = contacts_data.split('.')[0] + '_tail.tsv'

# Read data
contacts_df = pd.read_csv(contacts_data, comment="#", header=None,\
        names=['frame', 'type', 'protein_atom', 'pip2_atom'], delim_whitespace=True)
contacts_df_copy = contacts_df.copy()

print("\n### Raw contacts data \n")
print(contacts_df_copy)

pip2_headgroup_atoms = ['C11', 'C12', 'C13', 'C14', 'C15', 'C16', 'O12', 'O2', 'O3', 'O4', 'O5', 'O6', 'P', 'P4', 'P5', 'O13', 'O14',\
                       'OP42', 'OP43', 'OP44', 'OP52', 'OP53', 'OP54']

glycerol_backbone_atoms = ['C1', 'C2', 'C3', 'O21', 'O31', 'O11', 'O22', 'O32', 'C31', 'C21']

tails = ['C32', 'C33', 'C34', 'C35', 'C36', 'C37', 'C38', 'C39', 'C310', 'C311', 'C312', 'C313', 'C314', 'C315', 'C316',\
        'C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C210', 'C211', 'C212', 'C213', 'C214', 'C215', 'C216', 'C217', 'C218']

# ref: https://github.com/getcontacts/getcontacts
contact_type = ['hblb', 'hbls', 'vdw']

contacts_df_copy['pip2_atom_new'] = contacts_df_copy['pip2_atom'].apply(lambda x: x.split(':')[-1])
head_contact_df = contacts_df_copy[ (contacts_df_copy['pip2_atom_new'].isin(pip2_headgroup_atoms)) & (contacts_df_copy['type'].isin(contact_type))]
glycerol_contact_df = contacts_df_copy[(contacts_df_copy['pip2_atom_new'].isin(glycerol_backbone_atoms)) & (contacts_df_copy['type'].isin(contact_type)) ]
tail_contact_df = contacts_df_copy[(contacts_df_copy['pip2_atom_new'].isin(tails)) & (contacts_df_copy['type'].isin(contact_type))]

head_contact_final = head_contact_df.drop('pip2_atom_new', axis=1)
glycerol_contact_final = glycerol_contact_df.drop('pip2_atom_new', axis=1)
tail_contact_final = tail_contact_df.drop('pip2_atom_new', axis=1)


print("saving output...")
with open(headgroup_contact_output, 'w') as f:
    header1 = "# total_frames:12001 beg:0 end:12000 stride:1 interaction_types:sb,vdw,hb\n"
    header2 = "# Columns: frame, interaction_type, atom_1, atom_2[, atom_3[, atom_4]]\n"
    f.write(header1)
    f.write(header2)
    head_contact_final.to_csv(f, sep='\t', index=False, mode='a', header = False)

with open(glycerolBB_contact_output, 'w') as f:
    header1 = "# total_frames:12001 beg:0 end:12000 stride:1 interaction_types:sb,vdw,hb\n"
    header2 = "# Columns: frame, interaction_type, atom_1, atom_2[, atom_3[, atom_4]]\n"
    f.write(header1)
    f.write(header2)
    glycerol_contact_final.to_csv(f, sep='\t', index=False, mode='a', header = False)

with open(tail_contact_output, 'w') as f:
    header1 = "# total_frames:12001 beg:0 end:12000 stride:1 interaction_types:sb,vdw,hb\n"
    header2 = "# Columns: frame, interaction_type, atom_1, atom_2[, atom_3[, atom_4]]\n"
    f.write(header1)
    f.write(header2)
    tail_contact_final.to_csv(f, sep='\t', index=False, mode='a', header = False)

