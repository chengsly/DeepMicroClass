import os, sys, optparse
import keras
from tensorflow.keras.models import load_model
from encoding_model import EncodingScheme
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import constants
import multiprocessing


prog_base = os.path.split(sys.argv[0])[1]
parser = optparse.OptionParser()

parser.add_option("-i", "--in", action = "store", type="string", dest="inputFile", help="input fasta file")
parser.add_option("-d", "--modelDir", action="store", type="string", dest="modelDir", help="model directory for prediction")
parser.add_option("-e", "--encode", action="store", type="string", dest="encoding", help="encoding type, one-hot or codon")
parser.add_option("-m", "--mode", action="store", type="string", dest="predictionMode", help="prediction mode: single or hybrid")
parser.add_option("-l", "--length", action="store", type="int", dest="modelLength", help="in single mode, optionally choose one model with length")
parser.add_option("-o", "--outputDir", action="store", type="string", dest="outputDir", help="output directory for predicted results")

(options, args) = parser.parse_args()

if (options.inputFile is None or options.modelDir is None or 
    options.encoding is None or options.predictionMode is None):
    sys.stderr.write(prog_base + "ERROR: missing command line argument\n")
    parser.print_help()
    sys.exit(0)

if (options.encoding != "one-hot" and options.encoding != "codon"):
    sys.stderr.write(prog_base + "ERROR: must be codon-encoded or one-hot encoded\n")
    parser.print_help()
    sys.exit(0)

if (options.predictionMode != "single" and options.predictionMode != "hybrid"):
    sys.stderr.write(prog_base + "ERROR: must be single or hybrid mode\n")
    parser.print_help()
    sys.exit(0)

if (options.predictionMode == "single" and options.modelLength is None):
    sys.stderr.write(prog_base + "ERROR: specify a length from 500, 1000, 2000, 3000, 5000 when using single mode")
    parser.print_help()
    sys.exit(0)
    
if (options.outputDir is None):
    outputDir = os.path.join(os.getcwd(), "DeepMicrobeFinder_results")
else:
    outputDir = options.outputDir

inputFile = options.inputFile
if not os.path.exists(outputDir):
    os.makedirs(outputDir)

encoding = options.encoding
mode = options.predictionMode

####################################################################################################
#                                      STEP 1                                                      #
####################################################################################################


# Load Models: 500bp, 1000bp, 2000bp, 3000bp, 5000bp
modelDir = options.modelDir
print("Step 1/3: Loading models from {}".format(options.modelDir))
models = {}
model500Name = [x for x in os.listdir(modelDir) if "500.h5" in x][0]
model2000Name = [x for x in os.listdir(modelDir) if "2000.h5" in x][0]
model1000Name = [x for x in os.listdir(modelDir) if "1000.h5" in x][0]
model3000Name = [x for x in os.listdir(modelDir) if "3000.h5" in x][0]
model5000Name = [x for x in os.listdir(modelDir) if "5000.h5" in x][0]

if model500Name is None or model2000Name is None or model1000Name is None or model3000Name is None or model5000Name is None:
    print("ERROR: models for length 500, 1000, 2000, 3000, 5000 must be present")
    exit(0)

models["500"] = load_model(os.path.join(modelDir, model500Name))
models["1000"] = load_model(os.path.join(modelDir, model1000Name))
models["2000"] = load_model(os.path.join(modelDir, model2000Name))
models["3000"] = load_model(os.path.join(modelDir, model3000Name))
models["5000"] = load_model(os.path.join(modelDir, model5000Name))


####################################################################################################
#                                      STEP 2                                                      #
####################################################################################################
print("Step 2/3: Encode and predict...")
encodingModel = EncodingScheme(encoding, "DNA")

# helper function for single mode prediction
def predict_chunk_routine(fwdarr, bwdarr):
    sumEukScore = 0
    sumEukVirusScore = 0
    sumPlasmidScore = 0
    sumProkScore = 0
    sumProkVirScore = 0
    for i in range(int(len(seq)/5000)): # in this case we care only about 5000 chunks
        startIdx = i * 5000
        endIdx = (i+1) * 5000
        scores = models["5000"].predict([fwdarr[:,startIdx:endIdx,:], bwdarr[:,startIdx:endIdx,:]], batch_size=1)[0]
        sumEukScore += scores[0]
        sumEukVirusScore += scores[1]
        sumPlasmidScore += scores[2]
        sumProkScore += scores[3]
        sumProkVirScore += scores[4]
        return [sumEukScore, sumEukVirusScore, sumPlasmidScore, sumProkScore, sumProkVirScore]

# hybrid mode prediction
# for each hybrid-encoded sequence, output the total score
# greedy principle: always use the longest model if possible
# scores from strong predictors have more weight
def predict(encodedSeqfw, encodedSeqbw):
    # input is hybrid encoded: each item has different length
    sumEukScore = 0
    sumEukVirusScore = 0
    sumPlasmidScore = 0
    sumProkScore = 0
    sumProkVirScore = 0
    for i in range(len(encodedSeqfw)):
        scores=[]
        fwdarr = np.array([encodedSeqfw[i]])
        bwdarr = np.array([encodedSeqbw[i]])
        if (len(encodedSeqfw[i]) == 5000):
            scores = 5 * models["5000"].predict([fwdarr,bwdarr], batch_size=1)[0]
        elif len(encodedSeqfw[i]) == 3000:
            scores = 3 * models["3000"].predict([fwdarr,bwdarr], batch_size=1)[0]
        elif len(encodedSeqfw[i]) == 2000:
            scores = 2 * models["2000"].predict([fwdarr,bwdarr], batch_size=1)[0]
        elif len(encodedSeqfw[i]) == 1000:
            scores = 1 * models["1000"].predict([fwdarr,bwdarr], batch_size=1)[0]
        else:
            scores = 0.5 * models["500"].predict([fwdarr,bwdarr], batch_size=1)[0]  
        sumEukScore += scores[0]
        sumEukVirusScore += scores[1]
        sumPlasmidScore += scores[2]
        sumProkScore += scores[3]
        sumProkVirScore += scores[4]

    return [sumEukScore, sumEukVirusScore, sumPlasmidScore, sumProkScore, sumProkVirScore]


# single model preodict: predict the whole contig with one model, depending on the length
def predict_single(seq, length):
    # set up single mode encoding - normal encoding
    encodefw = encodingModel.encodeSeq(seq)
    seqR = Seq(seq).reverse_complement()
    encodebw = encodingModel.encodeSeq(seqR)
    seqLength = len(seq)
    # determine number of iterations for a fixed length chunk.
    count_iter = seqLength // length
    remainder_length = seqLength % length
    if remainder_length >= (2/3 * length):
        count_iter += 1
    sumEukScore = 0
    sumEukVirusScore = 0
    sumPlasmidScore = 0
    sumProkScore = 0
    sumProkVirScore = 0
    for i in range(count_iter):
        scores = []
        start_idx = i * length
        end_idx = (i+1) * length
        if end_idx >= seqLength:
            end_idx = seqLength
        fwdarr = np.array([encodefw[start_idx:end_idx]])
        bwdarr = np.array([encodebw[start_idx:end_idx]])
        scores = models[str(length)].predict([fwdarr, bwdarr], batch_size=1)[0]
        sumEukScore += scores[0]
        sumEukVirusScore += scores[1]
        sumPlasmidScore += scores[2]
        sumProkScore += scores[3]
        sumProkVirScore += scores[4]
    return [sumEukScore, sumEukVirusScore, sumPlasmidScore, sumProkScore, sumProkVirScore]


# steps for reading input file and start predictions
scores = []
names = []

with open(inputFile, 'r') as faLines :
    code = []
    codeR = []
    seqname = []
    head = ''
    lineNum = 0
    seq = ''
    flag = 0
    for line in faLines :
        lineNum += 1
        if flag == 0 and line[0] == '>' :
            head = line.strip()[1:]
            continue
        elif line[0] != '>' :
            seq = seq + line.strip()
            flag += 1
        elif flag > 0 and line[0] == '>' :
            countN = seq.count("N")
            if countN/len(seq) <= 0.3:
                # single mode prediction
                if mode == "single":
                    if len(seq) < 500:
                        print("contig {} length less than 500: ignore...".format(head))
                    else:
                        score = predict_single(seq, options.modelLength)
                        scores.append(score)
                        names.append(head)
                
                # hybrid mode prediction
                if mode == "hybrid":
                    if (len(seq) < 500):
                        print("contig {} length less than 500: ignore...".format(head))
                    else:
                        encodefw, encodebw = encodingModel.encodeContig(seq)
                        pred_score = predict(encodefw, encodebw)
                        scores.append(pred_score)
                        names.append(head)
            else :
                if countN/len(seq) > 0.3 :
                    print("   {} has >30% Ns, skipping it".format(head))
            
            flag = 0
            seq = ''
            head = line.strip()[1:]

    if flag > 0 :
        countN = seq.count("N")
        if countN/len(seq) <= 0.3: 
            # single mode prediction
            if mode == "single":
                if len(seq) < 500:
                    print("contig {} length less than 500: ignore...".format(head))
                else: 
                    score = predict_single(seq, options.modelLength)
                    scores.append(score)
                    names.append(head)
            
            # hybrid mode prediction
            if mode == "hybrid":
                if len(seq) < 500:
                    print("contig length {} less than 500: ignore...".format(head))
                else:
                    encodefw, encodebw = encodingModel.encodeContig(seq)
                    pred_score = predict(encodefw, encodebw)
                    scores.append(pred_score)
                    names.append(head)


####################################################################################################
#                                      STEP 3                                                      #
####################################################################################################
# summarize the result
print("Step 3/3: Summarize and Save to File...")

inputFileName = os.path.basename(inputFile)
modelLengthName = "all" if options.modelLength == None else options.modelLength
resultFileName_single = "{}_pred_{}_single_{}.txt".format(inputFileName, encoding, modelLengthName)
resultFileName_hybrid = "{}_pred_{}_hybrid.txt".format(inputFileName, encoding)

resultFileName = resultFileName_single if mode == "single" else resultFileName_hybrid
resultFilepath = os.path.join(outputDir, resultFileName)

with open(resultFilepath, "w") as f:
    f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(constants.NAME,constants.EUK, constants.EUKVIR, constants.PLASMID, constants.PROK, constants.PROKVIR))
    for i in range(len(scores)):
        f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(names[i],scores[i][0], scores[i][1], scores[i][2], scores[i][3], scores[i][4]))

print("Result saved in {}.".format(resultFilepath))

print("Cleaning up...")
for item in os.listdir(outputDir):
    if item.endswith(".npy"):
        os.remove(os.path.join(outputDir, item))

print("Prediction Finished.")
