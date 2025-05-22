# Running the script  find_interface

find_interface is a script that using the output provided by Alphafold3, provides different plots and files useful for the identification of interacting interfaces

## Installation

Packages required to run the find_interface script are listed in the conda environment file environment.yml or requirements.txt 

Conda environment can be activated with the instruction:

```console
conda env create -f environment.yml --name alphabridge

conda activate alphabridge
```

## How to run the script

1. Download and unzip the folder from the AlphaFold server.

To run the script from the find_interface working directory, use

```console
python3 define_interfaces.py -i path-to-AlphaFold3-folder
```

## DEMO EXAMPLE

```console
python3 define_interfaces.py -i /data/AF3/fold_1a02
```

The output will be saved inside a folder called Alphabrige, placed inside in the same path provided as an input.

## OUTPUT FILES

The following files will be saved inside the input folder:

1. alphabridge_data.json
2. Confidance-contact_plot.png
3. Confidence matrix.png
4. pae.png
5. contact_matrix.png
6. ribbon_plot.png

alphabridge_data.json

File in json format with all information regarding the predicted structure used as an input, and all interacting interfaces identified by a specific cut-off contact threshold.

Confidance-contact_plot.png

PNG image with a combination of the Predicted Merged Confidence (PMC) matrix (upper right half) and the Predicted Distance Error (PDE) matrix (lower left half). This two-dimensional graph where each token is represented along both the x and y axes. The upper section presents data on the overall confidence in the local structure and the relative positions of two tokens combining the PAE and pLDDT information. The lower section of the plot illustrates the closeness of pairwise contacts between two tokens in the predicted structure using the PDE matrix.

Confidence matrix.png

PNG image with a combination of the PAE and pLDDT information, in a new matrix we call the Predicted Merged Confidence (PMC)

pae.png

PNG image with Predicted Aligned Error (PAE) matrix

contact_matrix.png

PNG image with Predicted Distance Error (PDE) matrix 

ribbon_plot.png

PNG image displaying a circular layout showing all assigned interfaces and contact links under a specific threshold
