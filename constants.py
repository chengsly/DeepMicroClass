####################################################################################################
#                           constants for class names                                              #
####################################################################################################
PROK_TR = "ProkTr"
EUK_TR = "EukTr"
PROK_VAL = "ProkVal"
EUK_VAL = "EukVal"

PLASMID_TR = "PlasmidTr"
PLASMID_VAL = "PlasmidVal"

PROK_VIR_TR = "ProkVirusTr"
PROK_VIR_VAL = "ProkVirusVal"
EUK_VIR_TR = "EukVirusTr"
EUK_VIR_VAL = "EukVirusVal"

TEST_INPUT = "TestInput"
PROK = "Prokaryote"
EUK = "Eukaryote"
PLASMID = "Plasmid"
PROKVIR = "ProkaryoteVirus"
EUKVIR = "EukaryoteVirus"
PRED = "Predicted Class"


ALLOWED_TYPES = [PROK_TR, PROK_VAL, EUK_TR, EUK_VAL, PLASMID_TR, PLASMID_VAL, 
                 PROK_VIR_TR, PROK_VIR_VAL, EUK_VIR_TR, EUK_VIR_VAL, TEST_INPUT]

NAME="Sequence Name"
FORWARD_NAME = "codefw.npy"
REVERSE_NAME = "codebw.npy"

####################################################################################################
#                           constants for encodings                                                #
####################################################################################################

ENCODING_SCHEME = ['embedding', 'one-hot', 'codon']

# embeddings
EMBDDING_DIM = 100
DUMMY_VEC_EMBEDDING = [0.01] * 100

# One hot
ALLOWED_CHARS_DNA = ['A', 'C', 'T', 'G']
ALLOWED_CHARS_PROTEIN = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
                         'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
import numpy as np
ENCODE_A = np.array([1,0,0,0], dtype=np.int8)
ENCODE_C = np.array([0,1,0,0], dtype=np.int8)
ENCODE_G = np.array([0,0,1,0], dtype=np.int8)
ENCODE_T = np.array([0,0,0,1], dtype=np.int8)
ENCODE_UNK = [1/4, 1/4, 1/4, 1/4]

# Codon
ENCODE_UNK_CODON = [1/64] * 64

# Protein
PROTEINUNK = [1/20] * 20

# Codon table 
CODON_ENCODED = ['F', 'L', 'I', 'M', 'V', 'S', 'P', 'T', 'A', 'Y',
		 'H', 'Q', 'N', 'K', 'D', 'E', 'C', 'W', 'R', 'S',
		 'R', 'G']
CODON_TABLE_UNK = [1/22] * 22


# dimensions
DNA_DIM = 4
CODON_DIM = 22
CODON_FULL_DIM = 64
PROTEIN_DIM = 20
EMBEDDING_DIM = 100














