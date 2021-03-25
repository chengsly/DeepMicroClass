### Installation

dependencies: install using anaconda, using the tentative name, "DeepEukFinder".

```sh
conda create --name def python=3.7
conda activate def
pip install tensorflow keras numpy scipy pandas sklearn biopython
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
python encode.py -i <input_dir> -e one-hot -l 3000 -t DNA -o <output_dir>
```
