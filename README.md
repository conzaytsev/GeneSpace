# Gene Space
Python module for analysis of gene sequences.

## Requirements
`Python >3.9`, `numpy`, `scikit-lean`

## Citation


## Dataset
This model requires a gene dataset for the specific organism.
Dataset can be either presented as a fasta file with a set of gene sequences, or as a .csv file containing three columns: ID, expression number, gene sequence (compatible with the CodonExpressionIndex module https://github.com/conzaytsev/CodonExpressionIndex).

## Installation
From PyPI:

    pip install genespace
    
For the latest development version:
    
    pip install git+https://github.com/conzaytsev/GeneSpace.git

## Quickstart
### Importing the module

    import genespace

### Subspace creation
Organism specific Subspaces can be created based on the nucleotide / codon / codon pair frequencies:
nu is a parameter for the OneClassSVM function, which defines the maximum fraction of genes left outside the subspace.

    model = genespace.NucleotideSubspace(path_to_dataset, nu=0.1)

or

    model = genespace.CodonSubspace(path_to_dataset, nu=0.1)

or

    model = genespace.CodonPairSubspace(path_to_dataset, nu=0.1)

### Subspace boundary
Hyperplane coefficients for the subspace boundary:

    model.hyperplane

Data format: [[B], [A₁, A₂, A₃,..., Aₙ]] for the hyperplane equation A₁X₁ + A₂X₂ + A₃X₃ + ... + AₙXₙ + B = 0

### Subspace volume
Estimate for the subspace volume from the number of the Gene Space vertices inside of it
        
    model.volume()

### Fraction of training genes inside the subspace

    model.probability()
        
### Distance from a nucleotide sequence to the edge of the subspace
Positive values indicate that the sequence is inside the subspace, negative - on the outside.

    model.distance('ATG...')

### Distance from a nucleotide sequence to the center of the Gene Space

    model.center_distance('ATG...')
