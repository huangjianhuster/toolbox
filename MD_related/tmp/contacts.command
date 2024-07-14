get_dynamic_contacts.py --topology 5ns.pdb --trajectory s1_1ns.xtc --itype sb hb --sele "resname POPI and chain I J K L" --sele2 protein --output sb_hb.tsv
get_contact_frequencies.py --input_files sb_hb.tsv --output_file sb_hb_frq.tsv
grep -v "TIP3" sb_hb.tsv > sb_hb_trimWater.tsv
python contacts_seq.py sb_hb_frq.tsv


# for the raw trajectory
get_dynamic_contacts.py --topology 5ns.pdb --trajectory s1.xtc --itypes sb hb --sele "resname POPI and chain I J K L" --sele2 protein --output sb_hb_update.tsv
# remove water mediated interaction <-- Do not want to consider
grep -v "TIP3" sb_hb_update.tsv > sb_hb_trimWater.tsv
get_contact_frequencies.py --input_files sb_hb_trimWater.tsv --output_file sb_hb_trimWater_frq.tsv
python contacts_seq.py sb_hb_trimWater_frq.tsv


# P4 atom
get_dynamic_contacts.py --topology 5ns.pdb --trajectory s1.xtc --itypes sb hb --sele "resname POPI and chain I J K L and name P4" --sele2 protein --output P4_sb_hb.tsv &
grep -v "TIP3" P4_sb_hb.tsv > P4_sb_hb_trimWater.tsv


# P5 atom
get_dynamic_contacts.py --topology 5ns.pdb --trajectory s1.xtc --itypes sb hb --sele "resname POPI and chain I J K L and name P5" --sele2 protein --output P5_sb_hb.tsv &
grep -v "TIP3" P5_sb_hb.tsv > P5_sb_hb_trimWater.tsv


# P atom
get_dynamic_contacts.py --topology 5ns.pdb --trajectory s1.xtc --itypes sb hb --sele "resname POPI and chain I J K L and name P" --sele2 protein --output P_sb_hb.tsv &
grep -v "TIP3" P_sb_hb.tsv > P_sb_hb_trimWater.tsv
python contacts_feq.py sb_hb_frq.tsv


(for manuscript)
1. calculate PIP2 contacts with protein
get_dynamic_contacts.py --topology 5ns.pdb --trajectory s1.xtc --itypes sb vdw hb --sele "resname POPI and chain I J K L and not hydrogen" --sele2 "protein" --output s1_contacts.tsv &> log &

python group_pip2_contacts.py s1_contacts.tsv
get_contact_frequencies.py --input_files s1_contacts_headgroup.tsv --output_file s1_contacts_headgroup_frq.tsv
get_contact_frequencies.py --input_files s1_contacts_glycerolBB.tsv --output_file s1_contacts_glycerolBB_frq.tsv
get_contact_frequencies.py --input_files s1_contacts_tail.tsv --output_file s1_contacts_tail_frq.tsv

bash contacts.sh

2. plot
python plot_grouped_contacts.py
