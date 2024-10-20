# Running the script  find_interface

find_interface is a script that using the output provided by Alphafold3, provides  different plots and files useful for the identification of interacting interfaces

## Installation

Packages required to run the find_interface script are listed in the conda environment file find_interface_env.yml. Conda environment can be activated with the instruction:

```console
conda env create -n find_interface --file find_interface_env.yml

conda activate find_interface
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

The output will be saved in the same folder provided as an input.

## OUTPUT FILES

The following files will be saved inside the input folder:

1. Confidance-contact_plot.png
2. Binding_interfaces.csv
3. interface_df_per_token.csv
4. matrix_info.json

Confidance-contact_plot
PNG image with The Confidence-Contact Integration (CCI) plot (name subject to change). A two-dimensional graph where each token is represented along both the x and y axes. The plot is bisected diagonally into two sections. The upper section presents data on the overall confidence in the local structure and the relative positions of two tokens. The lower section of the Confidence-Contact Integration (CCI) plot illustrates the closeness of pairwise contacts between two tokens in the predicted structure.

Binding_interfaces

CSV file where all interfaces are listed, along with the involved entities involved in the interaction

Interface_df_per_token

CSV file with per-token level information (Maarten this is the file you use as an input for creating the ribbon plot)

matrix_info

Json file with all information used to run this script. (Maarten if you want to recreate the Confidance-contact_plot, you need to use the information stored under the keys 'confidance_matrix' and 'contact_matrix' inside the json file)
