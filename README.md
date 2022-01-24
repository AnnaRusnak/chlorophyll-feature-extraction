# chlorohyll-feature-extraction

This script iterates through a dataset containing all of the chlorophyll a and chlorophyll b binding proteins from the PDB databank and extracts features used for the random forest classifier. 

Using Bio.PDB, it extracts all of the atoms that lie within 5A radius from a particular chlorophyll. The atomics interactions between the Carbon, Nitrogen, Oxygen, and Sulfur atoms of a particular chlorophyll and neighboring C, N, O, and S atoms are quantified and saved. The resulting dataset contains 16 features in form of atomic interactions between chlorophyll atoms and their neighbors, and the type of chlorophyll as a label. 
