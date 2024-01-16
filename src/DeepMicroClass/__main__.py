import argparse
from .predict import predict
from importlib.resources import files
import requests


def download_file(url, destination):
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        with open(destination, 'wb') as file:
            for chunk in response.iter_content(chunk_size=1024):
                file.write(chunk)
    else:
        raise Exception(f"Failed to download file: status code {response.status_code}")
    

def predict_wrapper(input, model, output_dir, encoding="one-hot", mode="hybrid", single_len=1000, device="cuda"):
    if not model:
        url = "https://github.com/chengsly/DeepMicroClass/raw/master/src/DeepMicroClass/model.ckpt"
        local_checkpoint_file = "model.ckpt"

        # Download the file
        print(f"Downloading model from {url}")
        download_file(url, local_checkpoint_file)

        model_path = local_checkpoint_file
    else:
        model_path = model

    predict(input, model_path, output_dir, encoding, mode, device)


def main():
    parser = argparse.ArgumentParser(
        prog="DeepMicroClass",
        description="DeepMicroClass: A deep learning framework for classifying metagenomic sequences",
    )
    subparsers = parser.add_subparsers()

    parser_train = subparsers.add_parser("train", help="Train the model")
    parser_train.add_argument("--input", "-i", dest="input", help="Path to the input fasta file")
    parser_train.add_argument(
        "--log_prefix", "-l", dest="log_prefix", default="log", help="Prefix for the log directory"
    )

    parser_predict = subparsers.add_parser("predict", help="Predict the class of a sequence")
    parser_predict.add_argument("--input", "-i", dest="input", help="Path to the input fasta file", required=True)
    parser_predict.add_argument(
        "--output_dir", "-o", dest="output_dir", help="Path to the output directory", required=True
    )
    parser_predict.add_argument("--model", "-m", dest="model", help="Path to the trained model")
    parser_predict.add_argument(
        "--encoding", "-e", dest="encoding", help="Encoding method", choices=["onehot", "embedding"], default="one-hot"
    )
    parser_predict.add_argument(
        "--mode", "-md", dest="mode", help="Prediction mode", choices=["hybrid", "single"], default="hybrid"
    )
    parser_predict.add_argument(
        "--single-len", "-sl", dest="single_len", help="Length to use in the single mode", type=int, default=1000
    )
    parser_predict.add_argument(
        "--device", "-d", dest="device", help="Device to use", choices=["cpu", "cuda"], default="cuda"
    )
    parser_predict.set_defaults(func=predict_wrapper)

    args = parser.parse_args()
    # Convert args to a dictionary
    args_dict = vars(args)

    # Remove 'func' from the dictionary
    func = args_dict.pop("func")

    # Call the function with the remaining arguments
    func(**args_dict)


if __name__ == "__main__":
    main()
