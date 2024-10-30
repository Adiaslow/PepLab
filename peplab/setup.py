from setuptools import setup, find_packages
import os

setup(
    name="peplab",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "rdkit-pypi",
        "pandas",
        "numpy"
    ],
    author="Adam Murray",
    author_email="admmurra@ucsc.edu",
    description="A peptide library generation toolkit",
    long_description=open("README.md").read() if os.path.exists("README.md") else "",
    long_description_content_type="text/markdown",
    url="https://github.com/Adiaslow/PepLab",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.10",
)
