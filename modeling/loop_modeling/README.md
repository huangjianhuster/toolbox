# loop modeling
1. generate the ali file (a format in modeller)
    given the pdb file is example.pdb:
    ```python
    python align_input.py pdb_file
    ```
    This will generate a file named ```example.ali```, which contains basic information about the PDB, the model and the sequence (in modeller ali file format). check: [modeller aln file](https://salilab.org/modeller/tutorial/basic.html).

2. prepare your sequence and perform alignment
    The sequence of final model you want should also be stored in ali format (sequence.ali). The header should be something like below:
    ```
    >P1;code
    sequence:sequence_name:::::::0.00: 0.00
    RHEM...CAMVIFT*
    ```
    For proteins having multiple domains/chains, use "/" as the separator in between chain sequences and the final chain should still end with a "*".

3. perform sequence alignment
    ```python 
    python align_generate.py example.pdb sequence.ali
    ```
    This will generate ```seq_alignment.ali``` and ```seq_alignment.pap``` files.

    A careful check of the ```seq_alignment.ali``` is required. I usually will use [T-coffee](https://tcoffee.crg.eu/) to quickly align the sequences as a third party reference. Modeller use "-" to represent mismatch of gaps. Make sure all the gap alignments are good, and manually edit the ```seq_alignment.ali``` if required.
    The final seq_alignment.ali should have something like the following:

    ```
    >P1;code
    structureX:example.pdb:::::::0.00: 0.00
    MSEAAHVLITGAAGQIGYILSHWI-----YGDRQVYLHLLDIPPAMNRLTALTMELEDCAFPHLAGFVATTDPKA
    AFKDIDCAFLVASMPLKPGQVRADLISSNSVIFKNTGEYLSKWAKPSVKVLVIGNPDNTNCEIAMLHAKNLKPEN
    FSSLSMLDQNRAYYEVASKLGVDVKDVHDIIVWGNHGESMVADLTQATFTKEGKTQKVVDVLDHDYVFDTFFKKI
    GHRAWDILEHRGFTSAASPTKAAIQHMKAWLFGTAPGEVLSMGIPVPEGNPYGIKPGVVF-----VDKEGKIHVV
    EGFKVNDWLREKLDFTEKDLFHEKEIALNHLAQGG*

    >P1;sequence_code
    sequence:TvLDH:::::::0.00: 0.00
    MSEAAHVLITGAAGQIGYILSHWIASGELYGDRQVYLHLLDIPPAMNRLTALTMELEDCAFPHLAGFVATTDPKA
    AFKDIDCAFLVASMPLKPGQVRADLISSNSVIFKNTGEYLSKWAKPSVKVLVIGNPDNTNCEIAMLHAKNLKPEN
    FSSLSMLDQNRAYYEVASKLGVDVKDVHDIIVWGNHGESMVADLTQATFTKEGKTQKVVDVLDHDYVFDTFFKKI
    GHRAWDILEHRGFTSAASPTKAAIQHMKAWLFGTAPGEVLSMGIPVPEGNPYGIKPGVVFSFPCNVDKEGKIHVV
    EGFKVNDWLREKLDFTEKDLFHEKEIALNHLAQGG*
    ```
    Here, the ```code``` part was the sequence extracted from the PDB file which has two missing loops. the ```sequence_code``` part is what we want the final model to be, which is a completed structure without any missing residues.

4. build the final models
    ```python
    python align_build.py seq_alignment.ali code sequence_code
    ```
    Remember, ```code``` and ```sequence_code``` are from the ```seq_alignment.ali``` file mentioned above.
    ```

# Homology modeling
Homology modeling is similar to the above loop modeling code. Except in the ```align_build.py``` script, directly use ```automodel``` without any atom selection.
