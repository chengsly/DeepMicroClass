import numpy as np
import torch

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


def int2onehot(array, num_classes=None):
    if num_classes is None:
        n = np.max(array) + 1
    else:
        n = num_classes
    output = np.eye(int(n))[array.astype(int)]
    output[output[:, 4] == 1, :] = 0.25
    return output


def seq2onehot(seq, contig_length=None):
    """Convert genome sequence into one-hot matrix

    :param seq: genome sequence
    :type seq: str or Bio.Seq.Seq
    :param contig_length: If designated, sample a subsequence with length contig_length from the original sequence, defaults to None
    :type contig_length: int, optional
    :return: A one-hot matrix
    :rtype: np.ndarray
    """
    seq = str(seq)
    int_seq = seq2intseq(seq, contig_length)
    onehot = int2onehot(int_seq)
    return onehot.astype(np.uint8)

def sample_onehot(genome_onehot, length):
    """
    Sample a given length fragment from a genome onehot matrix, if the genome length is shorter than the given length,
    return a 0-padded matrix
    :param genome_onehot:
    :param length:
    :return:
    """
    genome_length = genome_onehot.shape[0]
    if genome_length > length:
        start_index = np.random.randint(0, genome_length - length)
        end_index = start_index + length
        onehot_sampled = genome_onehot[start_index:end_index, :]
    elif genome_length == length:
        onehot_sampled = genome_onehot
    else:
        # onehot_sampled = np.vstack((genome_onehot, np.zeros((length - genome_length, genome_onehot.shape[1]))))
        onehot_sampled = torch.nn.functional.pad(genome_onehot, (0, 0, 0, length - genome_length), "constant", 0)
    return onehot_sampled

def mutate_onehot(org_onehot, mu, delta):
    from numpy.random import default_rng
    onehot = np.copy(org_onehot)
    if len(onehot.shape) == 3:
        onehot = np.squeeze(onehot)
        org_onehot = np.squeeze(org_onehot)
    # print(onehot.shape)
    rng = default_rng()
    possible_mutation = [
        [1, 2, 3],
        [0, 2, 3],
        [0, 1, 3],
        [0, 1, 2]
    ]

    for i in range(4):
        mutation_indicator = np.logical_and(org_onehot[:, i]==1, rng.random(size=org_onehot.shape[0]) < delta)
        mutation_target = rng.choice(possible_mutation[i], size=np.sum(mutation_indicator))
        onehot[mutation_indicator, :] = 0
        onehot[np.where(mutation_indicator)[0], mutation_target] = 1
    
    if mu != 0:
        insertion_indicator = rng.random(size=onehot.shape[0]) < 0.5 * mu
        if np.sum(insertion_indicator) > 0:
            insertion_sequence = rng.integers(0, 4, np.sum(insertion_indicator))
            insertion_sequence = np.eye(4)[insertion_sequence]
            onehot = np.insert(onehot, np.where(insertion_indicator)[0], insertion_sequence, axis=0)

        deletion_indicator = rng.random(size=onehot.shape[0]) > 0.5 * mu
        onehot = onehot[deletion_indicator, :]
    return onehot