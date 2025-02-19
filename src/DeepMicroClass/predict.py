import os, sys, optparse

# from encoding_model import EncodingScheme
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import DeepMicroClass.constants as constants

import pytorch_lightning as pl
from .model.DeepMicroClass import DeepMicroClass, LightningDMC, DMFTransformer
import torch
import torch.nn as nn
import torch.nn.functional as F
from . import utils
from pathlib import Path
import pkgutil
import pandas as pd
import logging


def predict(input_path, model_path, output_dir, encoding="one-hot", mode="hybrid", device="cuda", cpu_thread=0):
    if torch.cuda.is_available():
        device = torch.device(device)
        print(f"Run DeepMicroClass with {device}.")
    elif not torch.cuda.is_available():
        print("GPU is not available on your device. Run DeepMicroClass with CPU.")
        device = torch.device("cpu")

    if cpu_thread > 0:
        torch.set_num_threads(cpu_thread)

    print("Step 1/3: Loading models from {}".format(model_path))

    model = LightningDMC.load_from_checkpoint(model_path, model=DeepMicroClass(), map_location=device)
    model.to(device)
    model.eval()

    model_cpu = LightningDMC.load_from_checkpoint(model_path, model=DeepMicroClass(), map_location='cpu')
    model_cpu.to("cpu")
    model_cpu.eval()

    if not os.path.exists(output_dir):
        logging.info("Creating output directory {}".format(output_dir))
        os.makedirs(output_dir)

    print("Step 2/3: Loading input sequences from {}".format(input_path))

    seq_records = SeqIO.parse(input_path, "fasta")

    print("Step 3/3: Predicting the class of the input sequences")

    input_fn = Path(input_path).name
    result_file_name = Path(f"{input_fn}_pred_{encoding}_{mode}.tsv")
    result_file_path = Path(output_dir) / result_file_name

    prediction_scores = []

    with torch.no_grad():
        for record in seq_records:
            seq = str(record.seq)
            n_percent = seq.upper().count("N") / len(seq)
            if n_percent > constants.N_THRESHOLD:
                print(f"Skipping {record.id} due to too many Ns ({n_percent:.2f})")
                score = [record.description] + [0] * 5
                prediction_scores.append(score)
                continue
            if len(seq) < constants.MIN_SEQ_LEN:
                print(f"Skipping {record.id} due to too short length ({len(seq)})")
                score = [record.description] + [0] * 5
                prediction_scores.append(score)
                continue

            if mode == "hybrid":
                score = predict_hybrid(seq, model, model_cpu, encoding, prediction_scores, device=device)
                score_list = score.tolist()
                prediction_scores.append([record.id] + score.tolist() + [str(int(score_list.index(max(score_list))+1))])
            elif mode == "single":
                pass
    result_df = pd.DataFrame(
        prediction_scores,
        columns=[constants.NAME, constants.EUK, constants.EUKVIR, constants.PLASMID, constants.PROK, constants.PROKVIR, "best_choice"],
    )
    result_df.to_csv(result_file_path, sep="\t", index=False)
    print(f"Prediction result saved to {result_file_path}")


def predict_hybrid(seq, model, model_cpu, encoding, prediction_score, device="cuda"):
    onehot_fw = utils.seq2onehot(seq, num_classes=5)[:, :4]
    start_idx = 0
    onehot_chunks = []
    score = np.zeros(5)
    for l in constants.POSSIBLE_LEN:
        if len(seq) < l:
            onehot_chunks.append(np.array([]).reshape(0, l, 4))
            continue
        onehot_chunks.append(onehot_fw[start_idx : start_idx + (len(seq) - start_idx) // l * l].reshape(-1, l, 4))
        start_idx += (len(seq) - start_idx) // l * l
        # onehot_chunks.append(onehot_fw[:len(seq)//l*l].reshape(-1, l, 4))
    # onehot_chunks.append(onehot_fw[start_idx:].reshape(-1, len(seq)-start_idx, 4))

    for i, chunk in enumerate(onehot_chunks):
        if chunk.shape[0] == 0:
            continue
        chunk = torch.tensor(chunk).float()[:, None, :, :].to(device)
        chunk_score = model(chunk)
        chunk_score = F.softmax(chunk_score, dim=1).cpu().numpy()
        chunk_score = chunk_score.sum(axis=0) * constants.SCORE_FACTOR[i]
        score += chunk_score

    # return score

    # min_value = np.min(score)
    # max_value = np.max(score)
    # normalized_score = (score - min_value) / (max_value - min_value)
    # return normalized_score

    # Normalize the cumulative scores into probabilities using softmax
    final_probabilities = F.softmax(torch.tensor(score),dim=0).numpy()
    return score
