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
nu is a parameter for the OneClassSVM, which defines the maximum fraction of genes that would not fit inside the subspace.

    model = genespace.NucleotideSubspace(path_to_dataset, nu=0.1)
    model = genespace.CodonSubspace(path_to_dataset, nu=0.1)
    model = genespace.CodonPairSubspace(path_to_dataset, nu=0.1)
        
To show the equation for the subspace boundary:

    model.hyperplane
        
Calculate the estimate for the subspace volume:
        
    model.volume()

Calculate the fraction of the training genes that fit inside the subspace:

    model.probability()
        
Calculate the distance from a nucleotide sequence to the edge of the subspace:
Positive values indicate that the sequence is inside the subspace, negative - on the outside.

    model.distance('ATG...')

Calculate the distance from a nucleotide sequence to the center of the Gene Space:

    model.center_distance('ATG...')
