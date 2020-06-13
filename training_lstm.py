import tensorflow as tf
from tensorflow.keras import layers, models, Model, Input
import numpy as np
import constants
from tensorflow.keras.callbacks import ModelCheckpoint, EarlyStopping
from tensorflow.keras.optimizers import Adam
import sklearn
import optparse
import utils
import os, sys
from sklearn import preprocessing


####################################################################################################################################
#                                                      Set up                                                                      #
####################################################################################################################################

POOL_FACTOR = 1
dropout_cnn = 0.1
dropout_pool = 0.1
dropout_dense = 0.1
learningrate = 0.001
channel_num = 100    # 100-dim dense vector
filters = [64, 128, 256]
kernels = [6, 3, 2]
ndense = 500
epochs = 50

#########################
#         ARGS          #
#########################
prog_base = os.path.split(sys.argv[0])[1]
parser = optparse.OptionParser()
parser.add_option("-i", "--input", action = "store", type = "string", dest="inputDir",
                                    help = "Specify input folder for data")
parser.add_option("-e", "--encode", action = "store", type = "string", dest = "encodingScheme",
									help = "Specify one-hot, codon, or embedding")
parser.add_option("-l", "--length", action = "store", type = int, dest = "contigLength",
									help = "contigLength")
parser.add_option("-t", "--type", action = "store", type = "string", dest = "sequenceType",
									help = "sequence type, Protein or DNA")
parser.add_option("-o", "--output", action = "store", type ="string", dest = "outputDir",
                                    help = "Output Directory")

(options, args) = parser.parse_args()
if (options.encodingScheme is None or options.contigLength is None 
    or options.sequenceType is None or options.outputDir is None or options.inputDir == None):
	sys.stderr.write(prog_base + ": ERROR: missing required command-line argument\n")
	parser.print_help()
	sys.exit(0)

if (options.sequenceType != 'Protein' and options.sequenceType != 'DNA'):
    sys.stderr.write(prog_base + ": ERROR: invalid sequence type")
    sys.exit(0)

if (options.encodingScheme not in constants.ENCODING_SCHEME):
    sys.stderr.write(prog_base + ": ERROR: invalid encoding type: {}".format(options.encodingScheme))
    parser.print_help()
    sys.exit(0)

if (options.encodingScheme == 'codon' and options.sequenceType == 'Protein'):
    sys.stderr.write(prog_base + ":ERROR: codon encode only applies to DNA sequences\n")
    sys.exit(0)

if (options.encodingScheme == 'embedding'):
    channel_num = constants.EMBEDDING_DIM
elif (options.encodingScheme == 'one-hot' and options.sequenceType == 'Protein'):
    channel_num = constants.PROTEIN_DIM
elif (options.encodingScheme == 'one-hot' and options.sequenceType == 'DNA'):
    channel_num = constants.DNA_DIM
elif (options.encodingScheme == 'codon' and options.sequenceType == 'DNA'):
    channel_num = constants.CODON_DIM

contigLength = options.contigLength
sequenceType = options.sequenceType
inputDir = options.inputDir
outputDir = options.outputDir

##########################################################################################################################################################
#                                                                   Models                                                                               #
##########################################################################################################################################################

outputModelName = "model_{}_{}.h5".format(options.sequenceType, options.encodingScheme)
checkpointer = ModelCheckpoint(os.path.join(options.outputDir, outputModelName), verbose=1,save_best_only=True)
earlystopper = EarlyStopping(monitor='val_accuracy', min_delta=0.0001, patience=10, verbose=1)

def get_output(input_layer, hidden_layers):
    output = input_layer
    for hidden_layer in hidden_layers:
        output = hidden_layer(output)
    return output

def buildModelBidirectionalPass():
    forward_input = Input(shape=(None, channel_num))
    reverse_input = Input(shape=(None, channel_num))
    hidden_layers = [
        layers.Conv1D(filters = filters[0], kernel_size = 6, activation='relu'),
        layers.MaxPooling1D(strides=2, pool_size=2),
        layers.BatchNormalization(),
        layers.Conv1D(filters = filters[1], kernel_size=3, activation="relu"),
        layers.MaxPooling1D(strides=1, pool_size=2),
        layers.BatchNormalization(),
        layers.Conv1D(filters = filters[2], kernel_size=2, activation='relu'),
        layers.GlobalMaxPooling1D(),
        layers.Dropout(0.1),
        layers.Flatten(),
        layers.Dense(500, activation='relu'),
        layers.Dropout(dropout_dense),
        layers.Dense(5, activation='softmax')
    ]
    forward_output = get_output(forward_input, hidden_layers) 
    reverse_output = get_output(reverse_input, hidden_layers)
    output = layers.average([forward_output, reverse_output])    
    model = Model(inputs=[forward_input, reverse_input], outputs=output)
    model.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])
    model.summary()
    return model


def buildModelLSTMBidrectional():
    forward_input = Input(shape=(None, channel_num))
    reverse_input = Input(shape=(None, channel_num))
    hidden_layers = [
        layers.Conv1D(filters=filters[0], kernel_size=6, activation='relu'),
        layers.MaxPooling1D(strides=2, pool_size=2),
        layers.Conv1D(filters=filters[1], kernel_size=3, activation='relu'),
        layers.MaxPooling1D(strides=2, pool_size=2),
        layers.Bidirectional(layers.LSTM(64)),
        layers.Dropout(0.5),
        layers.Dense(500, activation='relu'),
        layers.Dropout(0.1),
        layers.Dense(5, activation='softmax')
    ]
    forward_output = get_output(forward_input, hidden_layers)
    reverse_output = get_output(reverse_input, hidden_layers)
    output = layers.average([forward_output, reverse_output])
    model = Model(inputs=[forward_input, reverse_input], outputs=output)
    model.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])
    model.summary()
    return model



#######################################################################################################################################
#                                               Data                                                                                  #
#######################################################################################################################################

print("Processing Data...")
contigLengthk = contigLength/1000
if contigLengthk.is_integer() :
    contigLengthk = int(contigLengthk)
lengthK = "{}k".format(contigLengthk)

print("Processing sequences of length: {}k\n".format(lengthK))

# names
prokTrFileName = [x for x in os.listdir(inputDir) if constants.PROK_TR in x and lengthK in x and constants.FORWARD_NAME in x][0]
eukTrFileName = [x for x in os.listdir(inputDir) if constants.EUK_TR in x and lengthK in x and constants.FORWARD_NAME in x][0]
plasmidTrFileName = [x for x in os.listdir(inputDir) if constants.PLASMID_TR in x and lengthK in x and constants.FORWARD_NAME in x][0]
prokVirusTrFileName = [x for x in os.listdir(inputDir) if constants.PROK_VIR_TR in x and lengthK in x and constants.FORWARD_NAME in x][0]
eukVirusTrFileName = [x for x in os.listdir(inputDir) if constants.EUK_VIR_TR in x and lengthK in x and constants.FORWARD_NAME in x][0]

prokValFileName = [x for x in os.listdir(inputDir) if constants.PROK_VAL in x and lengthK in x and constants.FORWARD_NAME in x][0]
eukValFileName = [x for x in os.listdir(inputDir) if constants.EUK_VAL in x and lengthK in x and constants.FORWARD_NAME in x][0]
plasmidValFileName = [x for x in os.listdir(inputDir) if constants.PLASMID_VAL in x and lengthK in x and constants.FORWARD_NAME in x][0]
prokVirusValFileName = [x for x in os.listdir(inputDir) if constants.PROK_VIR_VAL in x and lengthK in x and constants.FORWARD_NAME in x][0]
eukVirusValFileName = [x for x in os.listdir(inputDir) if constants.EUK_VIR_VAL in x and lengthK in x and constants.FORWARD_NAME in x][0]

prokTrFileName_bw = utils.getBackwardName(prokTrFileName)
eukTrFileName_bw = utils.getBackwardName(eukTrFileName)
plasmidTrFileName_bw = utils.getBackwardName(plasmidTrFileName)
prokVirusTrFileName_bw = utils.getBackwardName(prokVirusTrFileName)
eukVirusTrFileName_bw = utils.getBackwardName(eukVirusTrFileName)

prokValFileName_bw = utils.getBackwardName(prokValFileName)
eukValFileName_bw = utils.getBackwardName(eukValFileName)
plasmidValFileName_bw = utils.getBackwardName(plasmidValFileName)
eukVirusValFileName_bw = utils.getBackwardName(eukVirusValFileName)
prokVirusValFileName_bw = utils.getBackwardName(prokVirusValFileName)



# Host - fw
print("Processing forward host data...")
print("Processing prokaryote training data from {}".format(prokTrFileName))
prokTr = np.load(os.path.join(inputDir, prokTrFileName), allow_pickle=True)
print("Processing eukaryote training data from {}".format(eukTrFileName))
eukTr = np.load(os.path.join(inputDir, eukTrFileName), allow_pickle=True)

# Host - bw
print("Processing reverse host data...")
print("Processing prokaryote training data from {}".format(prokTrFileName_bw))
prokTr_bw = np.load(os.path.join(inputDir, prokTrFileName_bw), allow_pickle=True)
print("Processing eukaryote training data from {}".format(eukTrFileName_bw))
eukTr_bw = np.load(os.path.join(inputDir, eukTrFileName_bw), allow_pickle=True)

# Plasmid 
print("Processing plasmid data")
plasmidTr = np.load(os.path.join(inputDir, plasmidTrFileName), allow_pickle=True)
print("Processing plasmid data - bw")
plasmidTr_bw = np.load(os.path.join(inputDir, plasmidTrFileName_bw), allow_pickle=True)

# Virus - fw
print("Processing virus data...")
print("Processing prokaryotevirus training data from {}".format(prokVirusTrFileName))
prokVirusTr = np.load(os.path.join(inputDir, prokVirusTrFileName), allow_pickle=True)
print("Processin eukaryotevirus training data from {}".format(eukVirusTrFileName))
eukVirusTr = np.load(os.path.join(inputDir, eukVirusTrFileName), allow_pickle=True)

# Virus - bw
print("Processing virus data...")
print("Processing prokaryotevirus training data from {}".format(prokVirusTrFileName_bw))
prokVirusTr_bw = np.load(os.path.join(inputDir, prokVirusTrFileName_bw), allow_pickle=True)
print("Processin eukaryotevirus training data from {}".format(eukVirusTrFileName_bw))
eukVirusTr_bw = np.load(os.path.join(inputDir, eukVirusTrFileName_bw), allow_pickle=True)

# VirusTr = np.concatenate((prokVirusTr, EukVirusTr), axis=0)
# print("Finished Processing Virus of shape {}".format(VirusTr.shape))
# del ProkaryoteVirusTr
# del EukaryoteVirusTr

# Validation
print("Processing Validation Data...")
prokVal = np.load(os.path.join(inputDir, prokValFileName), allow_pickle=True)
eukVal = np.load(os.path.join(inputDir, eukValFileName), allow_pickle=True)
plasmidVal = np.load(os.path.join(inputDir, plasmidValFileName), allow_pickle=True)
prokVirusVal = np.load(os.path.join(inputDir, prokVirusValFileName), allow_pickle=True)
eukVirusVal = np.load(os.path.join(inputDir, eukVirusValFileName), allow_pickle=True)

print("Processing Validation Data - bw...")
prokVal_bw = np.load(os.path.join(inputDir, prokValFileName_bw), allow_pickle=True)
eukVal_bw = np.load(os.path.join(inputDir, eukValFileName_bw), allow_pickle=True)
plasmidVal_bw = np.load(os.path.join(inputDir, plasmidValFileName_bw), allow_pickle=True)
prokVirusVal_bw = np.load(os.path.join(inputDir, prokVirusValFileName_bw), allow_pickle=True)
eukVirusVal_bw = np.load(os.path.join(inputDir, eukVirusValFileName_bw), allow_pickle=True)

# Training Data
labelsTr = ["Prokaryote"] * prokTr.shape[0] + ["Eukaryote"] * eukTr.shape[0] + \
    ["Plasmid"] * plasmidTr.shape[0] + ["ProkaryoteVirus"] * prokVirusTr.shape[0] + ["EukaryoteVirus"] * eukVirusTr.shape[0]

labelsVal = ["Prokaryote"] * prokVal.shape[0] + ["Eukaryote"] * eukVal.shape[0] + \
    ["Plasmid"] * plasmidVal.shape[0] + ["ProkaryoteVirus"] * prokVirusVal.shape[0] + ["EukaryoteVirus"] * eukVirusVal.shape[0]

X_tr = np.concatenate((prokTr, eukTr, plasmidTr, prokVirusTr, eukVirusTr), axis=0)
X_tr_bw = np.concatenate((prokTr_bw, eukTr_bw, plasmidTr_bw, prokVirusTr_bw, eukVirusTr_bw), axis=0)
X_val = np.concatenate((prokVal, eukVal, plasmidVal, prokVirusVal, eukVirusVal), axis=0)
X_val_bw = np.concatenate((prokVal_bw, eukVal_bw, plasmidVal_bw, prokVirusVal_bw, eukVirusVal_bw), axis=0)

del prokTr
del prokTr_bw
del eukTr
del eukTr_bw
del plasmidTr
del plasmidTr_bw
del prokVirusTr
del prokVirusTr_bw
del eukVirusTr
del eukVirusTr_bw

lb = preprocessing.LabelBinarizer()
Y_tr = lb.fit_transform(labelsTr)
Y_val = lb.fit_transform(labelsVal)

index_tr = list(range(0, X_tr.shape[0]))
np.random.shuffle(index_tr)
X_tr_shuf = X_tr[np.ix_(index_tr, range(X_tr.shape[1]), range(X_tr.shape[2]))]
X_tr_shuf_bw = X_tr_bw[np.ix_(index_tr, range(X_tr_bw.shape[1]), range(X_tr_bw.shape[2]))]
del X_tr
del X_tr_bw
Y_tr_shuf = Y_tr[index_tr]
del Y_tr
print(X_tr_shuf.shape)
print(X_tr_shuf_bw.shape)
print(Y_tr_shuf.shape)
print(X_val.shape)
print(X_val_bw.shape)
print(Y_val.shape)

############################################################################################################################
#                                             Training                                                                     #
############################################################################################################################

batch_size=100
pool_len1 = int((1000-500+1))

model = buildModelLSTMBidrectional()

model.fit(x = [X_tr_shuf, X_tr_shuf_bw], y = Y_tr_shuf, 
            batch_size=batch_size, epochs=60, verbose=2, 
            validation_data=([X_val, X_val_bw], Y_val), 
            callbacks=[checkpointer, earlystopper])
