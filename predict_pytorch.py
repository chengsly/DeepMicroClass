import os, sys, optparse
from encoding_model import EncodingScheme
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import constants
import multiprocessing
import argparse

import pytorch_lightning as pl
from model.DeepMicroClass import DMF, LightningDMF, DMFTransformer
import torch
import torch.nn as nn
import torch.nn.functional as F

import utils

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


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

parser = argparse.ArgumentParser(
    description='DeepMicrobeFinder: a deep learning tool for predicting contig origin',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument("-i", "--input", dest="input_file", type=str, required=True, help="Input fasta file")
parser.add_argument("-o", "--output", dest="output_dir", type=str, required=False, default="DeepMicrobeFinder_results", help="Output directory for predicted results")
parser.add_argument("--model", dest="model_path", type=str, required=False, default="model.ckpt", help="Path to the model")
parser.add_argument('--contig_length', dest="contig_length", type=int, required=False, default=5000, help="Contig length for prediction")

args = parser.parse_args()

input_file = args.input_file
output_dir = args.output_dir
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
model_path = args.model_path
contig_length = args.contig_length

####################################################################################################
#                                      STEP 1                                                      #
####################################################################################################


print("Step 1/3: Loading models from {}".format(options.modelDir))
models = {}

model = LightningDMF.load_from_checkpoint('data/pt_logs/checkpoint/epoch=2999-step=768000-val_f1=0.906-val_acc=0.907.ckpt', model=DMF())
# model = LightningDMF.load_from_checkpoint('data/pt_logs/checkpoint_new/epoch=729-step=186880-val_f1=0.948-val_acc=0.948.ckpt', model=DMF(), map_location=device)
model.to(device)
model.eval()

model_cpu = LightningDMF.load_from_checkpoint('data/pt_logs/checkpoint/epoch=2999-step=768000-val_f1=0.906-val_acc=0.907.ckpt', model=DMF())
model_cpu.eval()

models['500'] = model
models['1000'] = model
models['2000'] = model
models['3000'] = model
models['5000'] = model

####################################################################################################
#                                      STEP 2                                                      #
####################################################################################################
print("Step 2/3: Encode and predict...")
encodingModel = EncodingScheme(encoding, "DNA")

# hybrid mode prediction
# for each hybrid-encoded sequence, output the total score
# greedy principle: always use the longest model if possible
# scores from strong predictors have more weight
def predict(encodedSeqfw):
    # input is hybrid encoded: each item has different length
    sumEukScore = 0
    sumEukVirusScore = 0
    sumPlasmidScore = 0
    sumProkScore = 0
    sumProkVirScore = 0
    for i in range(len(encodedSeqfw)):
        scores=[]
        fwdarr = np.array([encodedSeqfw[i]])
        # print(fwdarr.shape)
        # fwdarr = utils.mutate_onehot(fwdarr, 0, 0.005)
        # fwdarr = fwdarr[None, :, :]
        with torch.no_grad():
            fwdarr = torch.from_numpy(fwdarr).float()[:,None,:,:].to(device)
            # fwdarr = fwdarr[:,:,:,np.newaxis]
            # bwdarr = bwdarr[:,:,:,np.newaxis]
            if (len(encodedSeqfw[i]) == 5000):
                scores = 5 * F.softmax(models["5000"](fwdarr), dim=1)[0]
            elif len(encodedSeqfw[i]) == 3000:
                scores = 3 * F.softmax(models["3000"](fwdarr), dim=1)[0]
            elif len(encodedSeqfw[i]) == 2000:
                scores = 2 * F.softmax(models["2000"](fwdarr), dim=1)[0]
            elif len(encodedSeqfw[i]) == 1000:
                scores = 1 * F.softmax(models["1000"](fwdarr), dim=1)[0]
            else:
                scores = 0.5 * F.softmax(models["500"](fwdarr), dim=1)[0]  
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
        scores = models[str(length)](fwdarr, batch_size=1)[0]
        sumEukScore += scores[0]
        sumEukVirusScore += scores[1]
        sumPlasmidScore += scores[2]
        sumProkScore += scores[3]
        sumProkVirScore += scores[4]
    return [sumEukScore, sumEukVirusScore, sumPlasmidScore, sumProkScore, sumProkVirScore]


# steps for reading input file and start predictions
scores = []
names = []

for record in SeqIO.parse(inputFile, "fasta"):
    seq = record.seq

    head = record.id

    seq = seq.upper()
    countN = seq.count("N")
    if countN / len(seq) > 0.3:
        print("{} has >30% Ns, skipping it".format(head))
        names.append(head)
        scores.append([0, 0, 0, 0, 0])
    elif len(seq) < 500:
        print("{} is too short(<500 bps), skipping it".format(head))
        names.append(head)
        scores.append([0, 0, 0, 0, 0])
    else:
        onehot = utils.seq2onehot(seq)
        if onehot.shape[0] < contig_length:
            continue
        else:
            onehot = onehot[:onehot.shape[0]//contig_length*contig_length, :4].reshape(-1, contig_length, 4)
        onehot = onehot[:, np.newaxis, :, :]
        onehot = torch.from_numpy(onehot).float()
        with torch.no_grad():
            try:
                pred = model(onehot.to(device))
                pred = pred.cpu().numpy()
            except:
                pred = model_cpu(onehot.to('cpu'))
                pred = pred.numpy()
        
        

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
            print(head)
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
                        scores.append([0, 0, 0, 0, 0])
                        names.append(head)
                    else:
                        score = predict_single(seq, options.modelLength)
                        scores.append(score)
                        names.append(head)
                
                # hybrid mode prediction
                if mode == "hybrid":
                    if (len(seq) < 500):
                        print("contig {} length less than 500: ignore...".format(head))
                        scores.append([0, 0, 0, 0, 0])
                        names.append(head)
                    else:
                        encodefw, encodebw = encodingModel.encodeContig(seq)
                        # encodefw = utils.seq2onehot(seq)
                        pred_score = predict(encodefw)
                        scores.append(pred_score)
                        names.append(head)
            else :
                if countN/len(seq) > 0.3 :
                    print("   {} has >30% Ns, skipping it".format(head))
                    scores.append([0, 0, 0, 0, 0])
                    names.append(head)
            
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
                    scores.append([0, 0, 0, 0, 0])
                    names.append(head)
                else: 
                    score = predict_single(seq, options.modelLength)
                    scores.append(score)
                    names.append(head)
            
            # hybrid mode prediction
            if mode == "hybrid":
                if len(seq) < 500:
                    print("contig length {} less than 500: ignore...".format(head))
                    scores.append([0, 0, 0, 0, 0])
                    names.append(head)
                else:
                    encodefw, encodebw = encodingModel.encodeContig(seq)
                    pred_score = predict(encodefw)
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
