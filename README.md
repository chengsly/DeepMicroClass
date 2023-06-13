### Installation

dependencies: install using anaconda, using the tentative name, "DeepMicrobeFinder".

```sh
conda env create -f DeepMicroClass.yml
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
python predict_pytorch.py -i test.fasta -e one-hot -d models/one-hot-models/ -m single -l 500 
python predict_pytorch.py -i test.fasta -e one-hot -d models/one-hot-models/ -m single -l 1000
python predict_pytorch.py -i test.fasta -e one-hot -d models/one-hot-models/ -m single -l 2000
python predict_pytorch.py -i test.fasta -e one-hot -d models/one-hot-models/ -m single -l 3000
python predict_pytorch.py -i test.fasta -e one-hot -d models/one-hot-models/ -m single -l 5000
python predict_pytorch.py -i test.fasta -e one-hot -d models/one-hot-models/ -m single
python predict_pytorch.py -i test.fasta -e one-hot -d models/one-hot-models/ -m hybrid
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



