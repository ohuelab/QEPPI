import os
from glob import glob
from setuptools import setup

exec(open("QEPPI/version.py").read())

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="QEPPI",
    version=__version__,
    author="blacktanktop",
    author_email="blacktanktopme@gmail.com",
    description="Calculation module of QEPPI",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ohuelab/QEPPI/",
    license="MIT",
    packages=["QEPPI"],
    install_requires=[
        "rdkit-pypi>=2021.3.1.5",
        "numpy>=1.19.5",
        "pandas>=1.1.5",
    ],
    test_suite="tests",
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
