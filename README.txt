## Modifications

See original Readme below. This version is modified to support running on arbitrary MSAs,
assuming they conform to a standard format.

This format requires the first sequence to be the seed, the last sequence to be a 'dummy'
(e.g. a duplicate of the second last sequence), and the labels to have the form
<LabelA>_<LabelB>_<PairId>_<SpeciesId>
Where SpeciesId is a numeric species identifier.

To run the code on an MSA, the command is:

```bash
matlab -nodisplay -r "MI_IPA_main(8,'Data/GSIC_GSID.fasta','Res')"
```

## Original Readme

This code contains a Matlab implementation of the MI-IPA as described in "Inferring interaction partners from protein sequences using mutual information", by Anne-Florence Bitbol, PLoS Comput Biol 14(11):e1006401 (2018) - https://doi.org/10.1371/journal.pcbi.1006401
Archived version: https://doi.org/10.5281/zenodo.1421781
Briefly, this is an iterative pairing algorithm (IPA) based on mutual information (MI) that aims to predict interaction partners among paralogs from two protein families, just from their sequences. Here it is applied to our "standard dataset" of 5064 sequences of cognate histidine kinases and response regulators from the P2CS database (http://www.p2cs.org/).

This algorithm is a MI-based variant of the IPA introduced in "Inferring interaction partners from protein sequences", by Anne-Florence Bitbol, Robert S. Dwyer, Lucy J. Colwell, and Ned S. Wingreen, Proc. Natl. Acad. Sci. U.S.A 113 (43) 12180-12185 (2016), DOI: 10.1073/pnas.1606762113.

In order to use the code, please run "MI_IPA_main" under Matlab.

The source code is freely available under the GNU GPLv3 license (unless otherwise indicated in specific files).

If you find this code useful for your research, please cite the associated reference, "Inferring interaction partners from protein sequences using mutual information", by Anne-Florence Bitbol, PLoS Comput Biol 14(11):e1006401 (2018) - https://doi.org/10.1371/journal.pcbi.1006401
