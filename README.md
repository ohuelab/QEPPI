<p align="center">
  <img src="https://user-images.githubusercontent.com/7370243/135420088-f616adc8-1e92-4d9b-8b53-0b863497244d.png"  width="400px">
</p>

# QEPPI
Quantitative Estimate Index for Compounds Targeting Protein-Protein Interactions

[![License](https://img.shields.io/badge/license-MIT-green?style=flat-square)](LICENSE)
[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Fohuelab%2FQEPPI&count_bg=%238EC9EE&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=Hits&edge_flat=true)](https://hits.seeyoufarm.com/)
[![GitHub Clones](https://img.shields.io/badge/dynamic/json?color=success&label=Clone&query=count&url=https://github.com/ohuelab/QEPPI?raw=True&logo=github)](https://github.com/ohuelab/QEPPI/)

## Calculation QEPPI with using Google Colab
We have made it so that you can use Google Colab to calculate QEPPI from SMILES without creating your own environment.   
If you have a lot of SMILES to calculate, please convert the SMILES to SDF files.  

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](http://colab.research.google.com/github/ohuelab/QEPPI/blob/main/notebook/QEPPI.ipynb)

## Mininal Environment Setup (Only QEPPI calculation)
We setup it on a Mac (macOS10.15.7), but I'm sure it will run fine on other platforms such as Linux.  

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
- Kosugi T, Ohue M. **Quantitative estimate of protein-protein interaction targeting drug-likeness**. In _Proceedings of The 18th IEEE International Conference on Computational Intelligence in Bioinformatics and Computational Biology (CIBCB 2021)_. (accepted)
ChemRxiv, Preprint. 2021. [doi:10.33774/chemrxiv-2021-psqq4-v2](https://doi.org/10.33774/chemrxiv-2021-psqq4-v2)
- Kosugi T, Ohue M. **Development of a quantitative estimate index for early-stage screening of compounds targeting protein-protein interactions**. (under revision)
