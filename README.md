# QEPPI
Quantitative Estimate of Protein-Protein Interaction Targeting Druglikeness

[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

## Calculation QEPPI with using Google Colab
We have made it so that you can use Google Colab to calculate QEPPI from SMILES without creating your own environment.   
If you have a lot of SMILES to calculate, please convert the SMILES to SDF files.  

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](http://colab.research.google.com/github/ohuelab/QEPPI/blob/main/notebook/QEPPI.ipynb)

## Mininal Environment Setup (Only QEPPI calculation)
We setup it on a Mac(macOS10.15.7), but I'm sure it will run fine on other platforms such as Linux.  

```
miniconda3-4.3.30
Python 3.7.9
rdkit=2020.09.1.0
numpy=1.20.1
pandas=1.2.3
```
We also confirmed that QEPPI works with the installation of RDKit by ```pip install rdkit-pypi``` instead of the conda environment. (see [notebook](https://github.com/ohuelab/QEPPI/blob/main/notebook/QEPPI.ipynb))

## Test
Test it when you are done with the setup.
If the test passes, the QEPPI calculation has been successfully performed.
(We used pytest version is 6.2.2)
```
pytest -v
```

## QEPPI Calculaton example
```
# for .sdf
python calc_QEPPI.py --sdf PATH_TO_YOUR_COMPOUND.sdf --out PATH_TO_OUTPUT.csv
```
```
# for .csv ("A column name of "SMILES" is required.")
python calc_QEPPI.py --csv PATH_TO_YOUR_COMPOUND.csv --out PATH_TO_OUTPUT.csv
```
If you want to implement your own, the following is a sample code.

```
import qeppi as ppi
from rdkit.Chem import SDMolSupplier

q = ppi.QEPPI_Calculator()
q.load("./model/QEPPI.model")

ppi_s = SDMolSupplier("PATH_TO_SDF/YOUR_COMPOUND.sdf")
ppi_mols = [mol for mol in ppi_s if mol is not None]
result = list(map(q.qeppi, ppi_mols))
```

## Reference
Kosugi T, Ohue M. **Quantitative estimate of protein-protein interaction targeting drug-likeness**. In _Proceedings of The 18th IEEE International Conference on Computational Intelligence in Bioinformatics and Computational Biology (CIBCB 2021)_. (accepted)
ChemRxiv, Preprint. 2021. [doi:10.33774/chemrxiv-2021-psqq4-v2](https://doi.org/10.33774/chemrxiv-2021-psqq4-v2)
