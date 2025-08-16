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
