# from dna2vec.multi_k_model import MultiKModel
from itertools import product
from numpy import array
from keras.utils import to_categorical
import sys
import constants
from Bio.Data import CodonTable
from Bio.Seq import Seq
import numpy as np
import os
import utils
import pgzip

class EncodingScheme:
    def __init__(self,encode, seqType,  modelurl=None):
        self.encode = encode
        self.seqType = seqType
        self.codons = [''.join(item) for item in list(product(constants.ALLOWED_CHARS_DNA, repeat = 3))]
        self.codonmap = CodonTable.unambiguous_dna_by_id[1].forward_table
        self.codonencode = dict(zip(constants.CODON_ENCODED, to_categorical([i for i in range(len(constants.CODON_ENCODED))])))
        self.proteinmap = dict(zip(constants.ALLOWED_CHARS_PROTEIN, to_categorical([i for i in range(len(constants.ALLOWED_CHARS_PROTEIN))])))

    def encodeCodon(self, seq):
        seq_code = list()
        for pos in range(len(seq)-3):
            threeMer = str(seq[pos:pos+3]).upper()
            if threeMer in self.codonmap:
                seq_code.append(self.codonencode[self.codonmap[threeMer]])
            else:
                seq_code.append(constants.CODON_TABLE_UNK)
        return seq_code

    def encodeOneHotDNA(self, seq):
        # seq_code = list()
        # for i in range(len(seq)):
        #     upperChar = seq[i].upper()
        #     if upperChar == 'A':
        #         seq_code.append(constants.ENCODE_A)
        #     elif upperChar == 'C':
        #         seq_code.append(constants.ENCODE_C)
        #     elif upperChar == 'G':
        #         seq_code.append(constants.ENCODE_G)
        #     elif upperChar == 'T':
        #         seq_code.append(constants.ENCODE_T)
        #     else:
        #         seq_code.append(constants.ENCODE_UNK)
        # # return seq_code
        # return np.vstack(seq_code)
        return utils.seq2onehot(seq)[:, :4]

 
    def encodeOneHotProtein(self, seq):
        seq_code = list()
        for i in range(len(seq)):
            if seq[i] in constants.ALLOWED_CHARS_PROTEIN:
                seq_code.append(self.proteinmap[seq[i]])
            else:
                seq_code.append(constants.PROTEINUNK)
        return seq_code
    

    # main encode function
    def encodeSeq(self, seq):
        if (self.seqType == 'Protein' and self.encode == 'one-hot'):
            return self.encodeOneHotProtein(seq)
        if (self.seqType == 'DNA' and self.encode == 'one-hot'):
            return self.encodeOneHotDNA(seq)
        if (self.seqType == 'DNA' and self.encode == 'codon'):
            return self.encodeCodon(seq)
        else:
            sys.stderr.write("ERROR: invalid encode type or sequence type: currently have {} and {}\n".format(self.seqType, self.encode))
            exit(0)


    # encode a fixed length of chunks
    def encodefixedchunks(self, seq, contigLength):
        if (len(seq) < contigLength):
            return [], []
        pos = 0
        posEnd = pos + contigLength 
        seqEncode = []
        seqEncodeR = []
        while posEnd <= len(seq):
            chunk = seq[pos: posEnd]
            chunk_code = self.encodeSeq(chunk)
            seqEncode.append(chunk_code)
            chunkR = Seq(chunk).reverse_complement()
            chunkR_code = self.encodeSeq(chunkR)
            seqEncodeR.append(chunkR_code)
            pos = posEnd
            posEnd += contigLength
        return seqEncode, seqEncodeR
         

    # batch encode a sequence when in prediction mode
    # input is a sequence of unknown length
    def encodeContig(self, seq):
        seqLen = len(seq)
        chunk5k = seqLen / 5000
        chunk5kres = seqLen % 5000 
        chunk3k = seqLen / 3000
        chunk3kres = seqLen % 3000
        chunk2k = seqLen / 2000
        chunk2kres = seqLen % 2000
        chunk1k = seqLen / 1000
        chunk1kres = seqLen % 1000
        chunk500 = seqLen / 500

        if chunk500 == 0:
            return [],[]
        if chunk1k == 0: # less than 1000
            return self.encodefixedchunks(seq, 500)
        if chunk2k == 0: # between 1000 and 2000
            chunk1klst, chunk1kRlst = self.encodefixedchunks(seq, 1000)
            chunk500lst, chunk500Rlst = self.encodefixedchunks(seq[-chunk1kres:], 500)
            return chunk1klst + chunk500lst, chunk1kRlst + chunk500Rlst
        if chunk3k == 0: # between 2000 and 3000
            chunk2klst, chunk2kRlst = self.encodefixedchunks(seq, 2000)
            chunk500lst, chunk500Rlst = self.encodefixedchunks(seq[-chunk2kres:], 500)
            return chunk2klst + chunk500lst, chunk2kRlst + chunk500Rlst
        if chunk5k == 0: # between 3000 and 5000
            chunk3klst, chunk3kRlst = self.encodefixedchunks(seq, 3000)
            chunk1klst, chunk1kRlst = self.encodefixedchunks(seq[-chunk3kres:], 1000)
            chunk500lst, chunk500Rlst = self.encodefixedchunks(seq[-chunk1kres:], 500)
            return chunk3klst + chunk1klst + chunk500lst, chunk3kRlst + chunk1kRlst + chunk500Rlst
        # above 5000 long contigs: N = 5000 * a + b
        chunk5klst, chunk5kRlst = self.encodefixedchunks(seq, 5000)
        # print(len(chunk5klst))
        # print(chunk5kres)
        if chunk5kres >= 3000: # remaining between 3000 and 5000
            chunk3klst, chunk3kRlst = self.encodefixedchunks(seq[-chunk5kres:], 3000)
            chunk1klst, chunk1kRlst = self.encodefixedchunks(seq[-chunk2kres:], 1000)
            chunk500lst, chunk500Rlst = self.encodefixedchunks(seq[-chunk1kres:], 500)
            if chunk5kres >= 4000:
                return chunk5klst + chunk3klst + chunk1klst + chunk500lst, chunk5kRlst + chunk3kRlst + chunk1kRlst + chunk500lst
            else:
                return chunk5klst + chunk3klst + chunk500lst, chunk5kRlst + chunk3kRlst + chunk500lst
        if chunk5kres >= 2000:  # remaining between 2000 and 3000 
            chunk2klst, chunk2kRlst = self.encodefixedchunks(seq[-chunk5kres:], 2000)
            chunk500lst, chunk500Rlst = self.encodefixedchunks(seq[-chunk1kres:], 500)
            return chunk5klst + chunk2klst + chunk500lst, chunk5kRlst + chunk2kRlst + chunk500Rlst
        if chunk5kres >= 1000:  # remaining between 1000 and 2000
            chunk1klst, chunk1kRlst = self.encodefixedchunks(seq[-chunk5kres:], 1000)
            chunk500lst, chunk500Rlst = self.encodefixedchunks(seq[-chunk1kres:], 500)
            return chunk5klst + chunk1klst + chunk500lst, chunk5kRlst + chunk1kRlst + chunk500Rlst
        if (chunk5kres >= 500): # remaining bewtween 500 and 1000 
            chunk500lst, chunk500Rlst = self.encodefixedchunks(seq[-chunk5kres:], 500)
            return chunk5klst + chunk500lst, chunk5kRlst + chunk500Rlst 
        else:
            return chunk5klst, chunk5kRlst         

    # batch encode a fasta file: main encode routine
    def encode_routine(self, fileName, contigLength, contigType, outDir):
        contigLengthk = contigLength / 1000
        if contigLengthk.is_integer():
            contigLengthk = int(contigLengthk)
        names=[]
        counts = []
        NCBIName = os.path.splitext((os.path.basename(fileName)))[0]
        with open(fileName, 'r') as faLines :
            code = []
            # codeR = []
            seqname = []
            head = ''
            lineNum = 0
            seqCat = ''
            flag = 0
            count = 0
            for line in faLines :
                if flag == 0 and line[0] == '>' :
                    lineNum += 1
                    head = line.strip()
                    names.append(head)
                    continue
                elif line[0] != '>' :
                    seqCat = seqCat + line.strip()
                    flag += 1
                    lineNum += 1
                elif flag > 0 and line[0] == '>' :
                    lineNum += 1
                    pos = 0
                    posEnd = pos + contigLength
                    while posEnd <= len(seqCat) :
                        contigName = contigType+"#"+NCBIName+"#"+str(contigLengthk)+"k#"+head.split('/')[-1]+"#"+str(pos)+"#"+str(posEnd)
                        seq = seqCat[pos:posEnd]
                        countN = seq.count("N")
                        if countN/len(seq) <= 0.3 : 
                            seqname.append(">"+contigName)
                            seqname.append(seq)
                            seq_code = self.encodeSeq(seq)
                            code.append(seq_code)
                            # seqR = Seq(seq).reverse_complement()
                            # seqR_code = self.encodeSeq(seqR)
                            # codeR.append(seqR_code)
                            count += 1
                        pos = posEnd
                        posEnd = pos + contigLength
                    counts.append(count)
                    count = 0
                    flag = 0
                    seqCat = ''
                    head = line.strip()
                    names.append(head)
                
            if flag > 0 :
                lineNum += 1
                count = 0
                pos = 0
                posEnd = pos + contigLength
                while posEnd <= len(seqCat) :
                    contigName = contigType+"#"+NCBIName+"#"+str(contigLengthk)+"k#"+head.split('/')[-1]+"#"+str(pos)+"#"+str(posEnd)
                    seq = seqCat[pos:posEnd]

                    countN = seq.count("N")
                    if countN/len(seq) <= 0.3 : 
                        seqname.append(">"+contigName)
                        seqname.append(seq)
                        seq_code = self.encodeSeq(seq)
                        code.append(seq_code)
                        # seqR = Seq(seq).reverse_complement()
                        # seqR_code = self.encodeSeq(seqR)
                        # codeR.append(seqR_code)
                        count += 1
                    pos = posEnd
                    posEnd = pos + contigLength
                counts.append(count)

        if len(code) > 0 :
            print("In block: last")
            codeFileNamefw = contigType+"#"+NCBIName+"#"+str(contigLengthk)+"k_seq"+str(len(code))+"_codefw.npy"
            # codeFileNamebw = contigType+"#"+NCBIName+"#"+str(contigLengthk)+"k_seq"+str(len(code))+"_codebw.npy"
            print("encoded sequences are saved in {}".format(codeFileNamefw))
            
            np.save( os.path.join(outDir, codeFileNamefw), np.array(code))
            # with pgzip.open(os.path.join(outDir, f'{codeFileNamefw}.gz'), mode='w', thread=16) as handle:
            #     np.save(handle, np.stack(code))

            # np.save( os.path.join(outDir, codeFileNamebw), np.array(codeR) )
        print("total number of whole sequences encoded: {}".format(len(counts)))
        print("total number of sequence chunks: {}".format(sum(counts)))
        return names,counts
