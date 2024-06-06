# DeepMicroClass

This is the repository for DeepMicroClass, a deep learning based method that classifies metagenomic contigs into five sequence classes, e.g., viruses infecting prokaryotic or eukaryotic hosts, eukaryotic or prokaryotic chromosomes.

The paper corresponding to this repository has been published at _NAR Genomics and Bioinformatics_. The exact code version used in the paper can be found at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10989619.svg)](https://doi.org/10.5281/zenodo.10989619). If you find DeepMicroClass useful, please consider to cite: 
- Shengwei Hou, Tianqi Tang, Siliangyu Cheng, Yuanhao Liu, Tian Xia, Ting Chen, Jed A Fuhrman, Fengzhu Sun, DeepMicroClass sorts metagenomic contigs into prokaryotes, eukaryotes and viruses, NAR Genomics and Bioinformatics, Volume 6, Issue 2, June 2024, lqae044, https://doi.org/10.1093/nargab/lqae044

Please direct any questions to [Dr. Fengzhu Sun](mailto:fsun@usc.edu).

## Installation
Consider using `virtualenv` or `conda` to create a virtual environment for clean environment. The following workflow works on our tested server with GPU Nvidia Quadro RTX 6000 and CUDA version 10.1. Since DeepMicroClass is using Pytorch, it is possible that our Pytorch version is not compiled with the version of your CUDA drive. Therefore, we suggest you install Pytorch compiled with your own device.

You can find the Pytorch that has been compiled with your version of the CUDA driver from [Previous PyTorch Versions | PyTorch](https://pytorch.org/get-started/previous-versions/). Most Pytorch versions are available only for specific CUDA versions. For example, pytorch=1.0.1 is not available for CUDA=9.2

Or, you can update your CUDA driver. To query the CUDA version that matches your driver, just open a terminal and run:
```sh
nvidia-smi
```
To check the CUDA version in your environment variables, which is the CUDA version you are probably using, just run:
```sh
nvcc -V
```


### To use GPU (CPU is also optional)

Create a new environment:
```sh
conda create -n DeepMicroClass_python=3.10 python=3.10
conda activate DeepMicroClass_python=3.10
```

Check if `pip` is under the new environment:
```sh
which pip
```
Your `pip` should be under the new environment:
```text
/home/xiatian/anaconda3/envs/DeepMicroClass_python=3.10/bin/pip
```
If not, install `pip` under the new environment:
```sh
conda install pip
```
Then check it again.


After that, install DeepMicroClass by running:
```sh
pip install DeepMicroClass
```


Then, install Pytorch packages suitable for CUDA=11.6 and python=3.10.13. 
```sh
conda install pytorch==1.13.1 torchvision==0.14.1 torchaudio==0.13.1 pytorch-cuda=11.6 -c pytorch -c nvidia
```
These versions are compatible with some Python versions like 3.10.13 and 3.9.18. If you were going to install Pytorch packages suitable for your CUDA, you should try other versions of Python in the first step.

If you encountered errors like:
```text
LibMambaUnsatisfiableError: Encountered problems while solving:
  - package torchvision-0.14.1-py310_cpu requires python >=3.10,<3.11.0a0, but none of the providers can be installed
  - package torchaudio-0.13.1-py310_cpu requires python >=3.10,<3.11.0a0, but none of the providers can be installed

Could not solve for environment specs
The following packages are incompatible
...
```
You may need to try other versions of Python.

### To use CPU ONLY
CPU mode works for `python=3.9-3.12`. The workflow is likely the same.
```sh
conda create -n DeepMicroClass_python=3.10 python=3.10
conda activate DeepMicroClass_python=3.10
```


Check if `pip` is under the new environment:
```sh
which pip
```
Your `pip` should be under the new environment:
```text
/home/your_username/anaconda3/envs/DeepMicroClass_python=3.10/bin/pip
```
If not, install `pip` under the new environment:
```sh
conda install pip
```
Then check it again.


After that, install DeepMicroClass by:
```sh
pip install DeepMicroClass
```

Install module `requests`:
```sh
pip install requests
```

### Testing
To test DeepMicroClass, You can simply run:
```sh
DeepMicroClass test
```
It will use `test.fa` as input, and use all default settings to classify the input sequences, then write the test result in the current working directory. 

To find the expected output, You can run:
```sh
DeepMicroClass -h
```

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

The most straightforward way to run the prediction is to use the following command:

```sh
DeepMicroClass predict -i <input_fasta> -o <output_dir>
```

### Prediction File Format

Each row of the output prediction file contains the scores from the neural network for class label assignment. The higher the score, the more likely the class label for a given sequence.  

* `Sequence Name`: the name of the sequence in the input fasta file
* `Eukaryote`: the score of the current sequence being Eukaryote.
* `EukaryoteVirus`: the score of the current sequence being EukaryoteVirus
* `Plasmid`: the score of the current sequence being Plasmid
* `Prokaryote`: the score of the current sequence being prokaryote
* `ProkaryoteVirus`: the score of the current sequence being ProkaryoteVirus  

As an example:

```tsv
Sequence Name   Eukaryote       EukaryoteVirus  Plasmid Prokaryote      ProkaryoteVirus
NC_008199#3k#12001#15000        0.0001718740027172454   1.2298997305038642e-05  0.1367506879128939      0.4601334227976679      5.402931783348322
```

To assign a label to a given sequence, we find the maximum of all the assigned scores. In the above case, `ProkaryoteVirus` has the highest score, so the sequence is assigned as belonging to the class `ProkaryoteVirus`.



