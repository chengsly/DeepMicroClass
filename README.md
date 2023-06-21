# DeepMicroClass

This is the repository for DeepMicroClass, a deep learning based method that classifies metagenomic contigs into five sequence classes, e.g., viruses infecting prokaryotic or eukaryotic hosts, eukaryotic or prokaryotic chromosomes.

The paper corresponding to this repository is available at [bioRxiv](https://www.biorxiv.org/content/10.1101/2021.10.26.466018v1).

Please direct any questions to [Dr. Fengzhu Sun](mailto:fsun@usc.edu).

## Installation

The software package can be installed using the following command:

```sh
pip install DeepMicroClass
```

Consider using `virtualenv` or `conda` to create a virtual environment for clean environment.

## Usage

### Prediction

``` text
usage: DeepMicroClass predict [-h] --input INPUT --output_dir OUTPUT_DIR [--model MODEL] [--encoding {onehot,embedding}] [--mode {hybrid,single}] [--single-len SINGLE_LEN] [--device {cpu,cuda}]

options:
  -h, --help            show this help message and exit
  --input INPUT, -i INPUT
                        Path to the input fasta file
  --output_dir OUTPUT_DIR, -o OUTPUT_DIR
                        Path to the output directory
  --model MODEL, -m MODEL
                        Path to the trained model
  --encoding {onehot,embedding}, -e {onehot,embedding}
                        Encoding method
  --mode {hybrid,single}, -md {hybrid,single}
                        Prediction mode
  --single-len SINGLE_LEN, -sl SINGLE_LEN
                        Length to use in the single mode
  --device {cpu,cuda}, -d {cpu,cuda}
                        Device to use
```

The most straight forward way to run the prediction is to use the following command:

```sh
DeepMicroClass predict -i <input_fasta> -o <output_dir>
```

### Prediction File Format

Each row of the output prediction file contains the scores from the neural network for class label assignment. The higher the score, the more likely the class label for a given sequence.  

* Sequence Name: the name of the sequence in the input fasta file
* Eukaryote: the score the neural network model assigns to this sequence being Eukaryote.
* EukaryoteVirus: the score for this sequence being EukaryoteVirus
* Plasmid: the score for this sequence being Plasmid
* Prokaryote: the score for this sequence being prokaryote
* ProkaryoteVirus: the score for this sequence being ProkaryoteVirus  

As an example:

```tsv
Sequence Name   Eukaryote       EukaryoteVirus  Plasmid Prokaryote      ProkaryoteVirus
NC_008199#3k#12001#15000        0.0001718740027172454   1.2298997305038642e-05  0.1367506879128939      0.4601334227976679      5.402931783348322
```

To assign a label to a given sequence, we find the maximum of all the assigned scores. In the above case, ProkaryoteVirus has the highest score, so the sequence is assigned as belonging to the class ProkaryoteVirus.
