from encoding_model import EncodingScheme
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import optparse
import constants
import os, sys

####################################################################################################
#                                       ARGS                                                       #
####################################################################################################
prog_base = os.path.split(sys.argv[0])[1]
parser = optparse.OptionParser()
parser.add_option("-i", "--filename", action = "store", type = "string", dest = "fileName",
									help = "fileName")
parser.add_option("-l", "--contigLength", action = "store", type = int, dest = "contigLength",
									help = "contigLength")
parser.add_option("-p", "--contigType", action = "store", type = "string", dest = "contigType",
									help = "contigType, must be one of the combination: (Prok|Euk|Plasmid|ProkVirus|EukVirus) + (Tr|Val)")
parser.add_option("-e", "--encode", action="store", type="string", dest = "encoding",
                                    help = "encoding scheme, one-hot, codon, embedding")
parser.add_option("-m", "--model", action="store", type="string", dest = "embeddingModel",
                                    help = "embeddingModel, input word2vec file" )
parser.add_option("-t", "--sequenceType", action="store", type="string", dest="sequenceType",
                                    help = "sequence type: protein or DNA")


(options, args) = parser.parse_args()
if (options.fileName is None or options.contigLength is None or options.encoding is None or options.sequenceType is None) :
	sys.stderr.write(prog_base + ": ERROR: missing required command-line argument")
	parser.print_help()
	sys.exit(0)

contigType = options.contigType
if (contigType not in constants.ALLOWED_TYPES):
    sys.stderr.write(prog_base + ": ERROR: unsupported type\n")
    parser.print_help()
    sys.exit(0)

contigLength = options.contigLength
contigLengthk = contigLength/1000
if contigLengthk.is_integer() :
    contigLengthk = int(contigLengthk)

if (options.encoding == "embedding" and options.embeddingModel == None): 
    sys.stderr.write(prog_base + " ERROR: must specify word2vec model file for embedding")
    parser.print_help()
    sys.exit(0)

fileName = options.fileName
fileDir = os.path.dirname(fileName)
encodeDir = "encode_" + options.encoding 
outDir = os.path.join(fileDir, encodeDir)
if not os.path.exists(outDir):
    os.makedirs(outDir)


encodingModel = EncodingScheme(options.encoding, options.sequenceType, options.embeddingModel)
encodingModel.encode_routine(fileName, contigLength, contigType, outDir)
print("Finished encoding.")
