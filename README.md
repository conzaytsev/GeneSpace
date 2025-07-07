# Gene Space
Python module for analysis of gene sequences.

## Requirements
`Python >3.9`, `numpy`, `scikit-lean`

## Dataset
This model requires a gene dataset for the specific organism.
Dataset can be either presented as a fasta file with a set of gene sequences, or as a .csv file containing three columns: ID, expression number, gene sequence (compatible with the CodonExpressionIndex module https://github.com/conzaytsev/CodonExpressionIndex).

## Installation
From PyPI:

    pip install genespace
    
For the latest development version:
    
    pip install git+https://github.com/conzaytsev/GeneSpace.git

## Quickstart
Importing the module:

    import genespace

Organism specific Subspaces can be created based on the nucleotide / codon / codon pair frequencies:
nu is a parameter for the OneClassSVM function, which defines the maximum fraction of genes left outside the subspace.

    model = genespace.NucleotideSubspace(path_to_dataset, nu=0.1)

or

    model = genespace.CodonSubspace(path_to_dataset, nu=0.1)

or

    model = genespace.CodonPairSubspace(path_to_dataset, nu=0.1)
        
To get the hyperplane coefficients for the subspace boundary:

Data format: [[B], [A_1, A2, A3,..., An]] for the hyperplane equation A1X1 + A2X2 + A3X3 + ... + AnXn + B = 0


    model.hyperplane
        
Estimate for the subspace volume:
        
    model.volume()

Calculate the fraction of the training genes that fit inside the subspace:

    model.probability()
        
Calculate the distance from a nucleotide sequence to the edge of the subspace:
Positive values indicate that the sequence is inside the subspace, negative - on the outside.

    model.distance('ATG...')

Calculate the distance from a nucleotide sequence to the center of the Gene Space:

    model.center_distance('ATG...')
