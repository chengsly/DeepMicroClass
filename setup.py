from setuptools import setup
from setuptools import find_packages

setup(
    name='DeepMicrobeFinder',
    version='1.0.0',
    description='DeepMicrobeFinder',
    author='Siliangyu Cheng, Shenwei Hou',
    author_email='siliangc@usc.edu',
    url='https://github.com/chengsly/DeepMicrobeFinder',
    license='MIT',
    install_requires=[
        'tensorflow==1.15', 
        'keras==2.2.4',
        'numpy', 
        'scipy==1.4.1',
        'pandas', 
        'sklearn', 
        'biopython',
        'h5py==2.10.0'
    ],
    packages=find_packages())
'D
