from DeepMicroClass.predict import predict
import os
from importlib.resources import files

__all__ = (["DeepMicroClass"],)


def test_predict():
    input_path=files("DeepMicroClass").joinpath("demo/test.fa")
    output_dir=files("DeepMicroClass").joinpath("demo/")
    assert (
        os.system(f'DeepMicroClass predict -i {input_path} -o {output_dir}')
        == 0
    )
