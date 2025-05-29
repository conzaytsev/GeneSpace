# Gene Space
Python module for analysis of gene sequences.

## Requirements
`Python >3.9`, `numpy`, `scikit-lean`

## Dataset
This model requires a gene dataset for the specific organism.
Dataset can be either presented as a fasta file with gene sequences, or as a .csv file containing three columns: protein ID, expression number, gene sequence. This data format is compatible with the CodonExpressionIndex module https://github.com/conzaytsev/CodonExpressionIndex.

## Installation
From PyPI:

    pip install genespace
    
For the latest development version:
    
    pip install git+https://github.com/conzaytsev/GeneSpace.git

## Quickstart
Importing the module:

    import genespace

Organism Subspace can be created based on the nucleotide / codon / codon pair frequencies:

    model = genespace.NucleotideSubspace(path_to_dataset, 0.1)

    model = genespace.CodonSubspace(path_to_dataset, 0.1)

    model = genespace.CodonPairSubspace(path_to_dataset, 0.1)
        
To use the Codon Productivity model:

    model = cei.CodonProductivity(path_to_dataset)
        
In order to use the default dataset for E. Coli ATCC 25922 leave the brackets empty.
        
    model = cei.CodonExpressionIndex()
    model = cei.CodonProductivity()

Codon scores for each of the models as a dict:

    model.scores
        
Predict number of protein copies per cell based on the nucleotide sequence:

    model.predict('ATG...')
