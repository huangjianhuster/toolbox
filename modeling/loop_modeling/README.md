# loop modeling

1. generate the ali file (a format in modeller)
   given the pdb file is example.pdb:

   ```python
   python align_input.py pdb_file
   ```

   This will generate a file named ``example.ali``, which contains basic information about the PDB, the model and the sequence (in modeller ali file format). check: [modeller aln file](https://salilab.org/modeller/tutorial/basic.html).
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

   This will generate ``seq_alignment.ali`` and ``seq_alignment.pap`` files.

   A careful check of the ``seq_alignment.ali`` is required. I usually will use [T-coffee](https://tcoffee.crg.eu/) to quickly align the sequences as a third party reference. Modeller use "-" to represent mismatch of gaps. Make sure all the gap alignments are good, and manually edit the ``seq_alignment.ali`` if required.
   If there is one or more gaps you do not want to model, in the "sequence_code" part, you need to make certain gaps using "-".

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

   Here, the ``code`` part was the sequence extracted from the PDB file which has two missing loops. the ``sequence_code`` part is what we want the final model to be, which is a completed structure without any missing residues.

4. build the final models
   We can use ``AutoModel`` from MODELLER (see: https://salilab.org/modeller/manual/node23.html) to build the loop region:

   ```python
   python align_build.py seq_alignment.ali code sequence_code
   ```

   Remember, ``code`` and ``sequence_code`` are from the ``seq_alignment.ali`` file mentioned above.

   We can also use ``LoopModel`` from MODELLER (see: https://salilab.org/modeller/manual/node35.html)

   ```python
   python align_build2.py seq_alignment.ali code sequence_code

   ```


    **Tailor the script to your system**: you need to change the residue selection to match the missing loop residue index inside the two scripts:

    ```python
    def select_atoms(self):
        # Select residue 30 ~ 40 in chain A
        return Selection(self.residue_range('30:A', '40:A'))
    ```

    Here, the script by default selects the loop between chain A residue index 30 ~ chain A residue index 40. Running it will build the loop between this range specifically. You can also define multiple ranges (check: `align_build2.py`)


    The difference between`AutoModel` and `LoopModel`: the `AutoModel` builds homology models using the given alignment and templates, though the modling region can be selected. The loop region will be automatically dealt with as part of the global model-building process, meaning the loop refinement is limited. This is usually used for quick and dirty loop building. The `LoopModel` is a derivative class of ` AutoModel` and it has all its methods. Moreover, `LoopModel` has built-in protocols for further loop refinement (more thorough sampling, scoring and energy optimization)). However, it takes longer time.

# Homology modeling

Homology modeling is similar to the above loop modeling code. Except in the ``align_build.py`` script, directly use ``automodel`` without any atom selection.
