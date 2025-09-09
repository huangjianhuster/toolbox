# Modeller script for simple tasks

## 1. get sequence from a PDB file

```python
python get_sequence_from_pdb.py [input_pdb_path]
```
example: `python get_sequence_from_pdb.py 1ubq.pdb`.

We get the following:

```text
{'A': 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'}`
```

## 2. fix sidechains and optimize side chain conformations

```python
python fix_opt_sidchains.py [input_pdb_path]
```

Of note, the given input pdb is better to have all backbone atoms. During optimization, backbone atoms will be by default fixed while side chains will move around to minimize the conformational energy.

## 3. extract backbone of frames from a xtc file
dependencies: MDAnalysis
scenario: extract hyres backbone and then use `fix_opt_sidechains.py` to fix side chains

```python
# randomly select frames from xtc
python extract_frames_from_xtc.py [PSF_PATH] [XTC_PATH] --nframes 3 --seed 42 --prefix random
# extract the 1st, 10th and 20th frames
python extract_frames_from_xtc.py [PSF_PATH] [XTC_PATH] --frames 1 10 20 --prefix output
```

### 3.1 use `convpdb.pl` to post-processing the extracted backbone frames
```bash
# renumber chain ID from the last letter of segment ID
convpdb.pl -chainfromseg [EXTRACTED_PDB] > ouput.pdb
# convert HSD back to HIS
sed -i s/HSD/HIS/g output.pdb
```

