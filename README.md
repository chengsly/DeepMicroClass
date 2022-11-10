### Installation

dependencies: install using anaconda, using the tentative name, "DeepMicrobeFinder".

```sh
conda create --name def python=3.6
conda activate def
pip install tensorflow==1.15 keras==2.2.4 numpy scipy pandas sklearn biopython torch pgzip
pip install 'h5py==2.10.0' --force-reinstall
pip install 'scipy==1.4.1' --force-reinstall
```

### Running Prediction

Several command line options are available for prediction:
* `-i`: input fasta file for prediction
* `-e`: encoding mode, currently support "one-hot" 
* `-d`: directory for models, currently use "models/one-hot-models"
* `-m`: prediction mode, "hybrid" or "single"
* `-l`: (optional), if choosing single mode, specify a length for model to use

A set of pretrained model with different lengths locate in the folder `model`.

Example usage: 

```sh
python predict.py -i test.fasta -e one-hot -d models/one-hot-models/ -m single -l 500 
python predict.py -i test.fasta -e one-hot -d models/one-hot-models/ -m single -l 1000
python predict.py -i test.fasta -e one-hot -d models/one-hot-models/ -m single -l 2000
python predict.py -i test.fasta -e one-hot -d models/one-hot-models/ -m single -l 3000
python predict.py -i test.fasta -e one-hot -d models/one-hot-models/ -m single -l 5000
python predict.py -i test.fasta -e one-hot -d models/one-hot-models/ -m single
python predict.py -i test.fasta -e one-hot -d models/one-hot-models/ -m hybrid
```

### Train your own model

to train your own model, first encode 

specify the encoding scheme and give directory for training and validation set. Commandline arguments:
```
Options:
  -h, --help            show this help message and exit
  -i INPUTDIR, --input=INPUTDIR
                        Specify input folder for data
  -e ENCODINGSCHEME, --encode=ENCODINGSCHEME
                        Specify one-hot, codon, or embedding
  -l CONTIGLENGTH, --length=CONTIGLENGTH
                        contigLength
  -t SEQUENCETYPE, --type=SEQUENCETYPE
                        sequence type, Protein or DNA
  -o OUTPUTDIR, --output=OUTPUTDIR
                        Output Directory
```

Example usage: 

__Encoding step__: use the batch script 
- first parameter as the input directory for fasta files
- second paramter as the encoding length
- third parameter as encoding scheme, one-hot or codon.

```sh
bash batch_encode.sh <input_dir> 3000 one-hot
```

encoded files will be located in the folder `encode` created under the same directory as input directory. Folder name will also contain the encoding scheme.

The input fasta file should have the name: 
* `EukTr.fasta` for eukaryote training fasta file, `EukVal.fasta` for validation.
* `ProkTr.fasta` for prokaryote training fasta file, `ProkVal.fasta` for valiation.
* `ProkVirusTr.fasta` for prokaryote virus training file, `ProkVirusVal.fasta`.
* `EukVirusTr.fasta` for eukaryote virus training file, `EukVirusVal.fasta` for eukaryote virus validation file.
* `PlasmidTr.fasta` for plasmid training file, `PlasmidVal.fasta` for plasmid validation file.

__Training step__: for training, we assume that the input folder contains all encoded code of each category (Prokaryote, Eukaryote, ProkVirus, EukVirus, Plasmid), with the naming convention listed above.

```
Options:
  -h, --help            show this help message and exit
  -i INPUTDIR, --input=INPUTDIR
                        Specify input folder for data
  -e ENCODINGSCHEME, --encode=ENCODINGSCHEME
                        Specify one-hot, codon, or embedding
  -l CONTIGLENGTH, --length=CONTIGLENGTH
                        contigLength
  -t SEQUENCETYPE, --type=SEQUENCETYPE
                        sequence type, Protein or DNA
  -o OUTPUTDIR, --output=OUTPUTDIR
                        Output Directory
```

example usage:
```sh
python training.py -i <input_dir> -e one-hot -l 3000 -t DNA -o <output_dir>
```

### Prediction File Format 

Each row of the output prediction file contains the scores from the neural network for class label assignment. The higher the score, the more likely the class label for a given sequence.  
- Sequence Name: the name of the sequence in the input fasta file 
- Eukaryote: the score the neural network model assigns to this sequence being Eukaryote. 
- EukaryoteVirus: the score for this sequence being EukaryoteVirus 
- Plasmid: the score for this sequence being Plasmid 
- Prokaryote: the score for this sequence being prokaryote
- ProkaryoteVirus: the score for this sequence being ProkaryoteVirus  

As an example:
```
Sequence Name   Eukaryote       EukaryoteVirus  Plasmid Prokaryote      ProkaryoteVirus
NC_008199#3k#12001#15000        0.0001718740027172454   1.2298997305038642e-05  0.1367506879128939      0.4601334227976679      5.402931783348322
```

To assign a label to a given sequence, we find the maximum of all the assigned scores. In the above case, ProkaryoteVirus has the highest score, so the sequence is assigned as belonging to the class ProkaryoteVirus. 



