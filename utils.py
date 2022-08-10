import numpy as np

def getBackwardName(forwardName):
    if (forwardName[-6:] != "fw.npy"):
        print("invalid input name")
        exit(0)
    backwardName = forwardName[:-6] + "bw.npy"
    return backwardName

def nt2index(nucleotide: str):
    nucleotide = nucleotide.upper()
    if nucleotide == "A":
        return 0
    elif nucleotide == "C":
        return 1
    elif nucleotide == "G":
        return 2
    elif nucleotide == "T":
        return 3
    else:
        return 4


def seq2intseq(seq, contig_length=None):
    """Convert genome sequence into indexed sequence

    :param seq: genome sequence
    :type seq: str
    :param contig_length: If designated, sample a subsequence with length contig_length from the original sequence, defaults to None
    :type contig_length: int, optional
    :return: A numpy array
    :rtype: np.ndarray
    """
    if contig_length:
        if len(seq) <= contig_length:
            seq_int = np.zeros(contig_length)
            seq_int[: len(seq)] = list(map(nt2index, str(seq)))
        else:
            start_idx = np.random.randint(len(seq) - contig_length)
            seq_int = np.array(list(map(nt2index, str(seq[start_idx : start_idx + contig_length]))))
        return seq_int
    else:
        return np.array(list(map(nt2index, str(seq))))


def int2onehot(array):
    n = np.max(array) + 1
    return np.eye(int(n))[array.astype(int)]


def seq2onehot(seq, contig_length=None):
    assert isinstance(seq, str), "Input should be str"
    int_seq = seq2intseq(seq, contig_length)
    onehot = int2onehot(int_seq)
    return onehot.astype(np.uint8)