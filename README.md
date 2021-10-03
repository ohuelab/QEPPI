<p align="center">
  <img src="https://user-images.githubusercontent.com/7370243/135420088-f616adc8-1e92-4d9b-8b53-0b863497244d.png"  width="400px">
</p>

# QEPPI
Quantitative Estimate Index for Compounds Targeting Protein-Protein Interactions

[![License](https://img.shields.io/badge/license-MIT-green?style=flat-square)](LICENSE)
[![GitHub Clones](https://img.shields.io/badge/dynamic/json?style=flat-square?color=success&label=Clones_in_14days&query=count&url=https://github.com/ohuelab/QEPPI/blob/main/clone.json?raw=True&logo=github)](https://github.com/ohuelab/QEPPI/)
[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Fohuelab%2FQEPPI&count_bg=%238EC9EE&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=Hits&edge_flat=true)](https://hits.seeyoufarm.com/)
![PyPI](https://img.shields.io/pypi/v/QEPPI?style=flat-square)
[![Python Versions](https://img.shields.io/pypi/pyversions/QEPPI.svg)](https://pypi.org/project/QEPPI/)
[![tests](https://github.com/ohuelab/QEPPI/actions/workflows/tests.yml/badge.svg)](https://github.com/ohuelab/QEPPI)

## Calculation QEPPI with using Google Colab
We have made it so that you can use Google Colab to calculate QEPPI from SMILES without creating your own environment.   
If you have a lot of SMILES to calculate, please convert the SMILES to SDF files.  

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](http://colab.research.google.com/github/ohuelab/QEPPI/blob/main/notebook/QEPPI.ipynb)

## Mininal environment setup
We setup it on a Mac (macOS10.15.7), but I'm sure it will run fine on other platforms such as Linux.  

```bash
# Python 3.7 | 3.8
# dependencies
pip install rdkit-pypi # >= 2021.3.1.5
pip install numpy      # >= 1.19.5
pip install pandas     # >= 1.1.5
```

We also confirmed that QEPPI works with Colab. (see [notebook](https://github.com/ohuelab/QEPPI/blob/main/notebook/QEPPI.ipynb))

### Clone QEPPI
Clone QEPPI repository when you are done with the setup.

```bash
git clone https://github.com/ohuelab/QEPPI.git
```

### Test
Test it after git clone the QEPPI repository. If the test passes, the QEPPI calculation has been successfully performed. (We used pytest version is 6.2.2)  
```bash
cd QEPPI
pytest -v
```

## QEPPI calculation example
```bash
# for .sdf
python calc_QEPPI.py --sdf PATH_TO_YOUR_COMPOUND.sdf --out PATH_TO_OUTPUT.csv
```
```bash
# for .csv ("A column name of "SMILES" is required.")
python calc_QEPPI.py --csv PATH_TO_YOUR_COMPOUND.csv --out PATH_TO_OUTPUT.csv
```

## Instalation using pip install
You can install it with ```pip install QEPPI```. First, you need to install with the dependencies (see [Mininal environment setup](https://github.com/ohuelab/QEPPI#mininal-environment-setup)). The following sample code is available as an implementation example.  
```bash
# QEPPI
pip install QEPPI
```

```python
import QEPPI as ppi
from rdkit import Chem
from rdkit.Chem import SDMolSupplier

q = ppi.QEPPI_Calculator()
q.read()

# SMILES
smiles = "COC1=CC(=CC=C1NC(=O)[C@@H]1N[C@@H](CC(C)(C)C)[C@@](C#N)([C@H]1C1=CC=CC(Cl)=C1F)C1=CC=C(Cl)C=C1F)C(O)=O"
mol = Chem.MolFromSmiles(smiles)
print(q.qeppi(mol))
# 0.7862842663145835

# SDF
ppi_s = SDMolSupplier("PATH_TO_SDF/YOUR_COMPOUND.sdf")
ppi_mols = [mol for mol in ppi_s if mol is not None]
result = list(map(q.qeppi, ppi_mols))
```

## Reference
- Kosugi T, Ohue M. **Quantitative estimate of protein-protein interaction targeting drug-likeness**. In _Proceedings of The 18th IEEE International Conference on Computational Intelligence in Bioinformatics and Computational Biology (CIBCB 2021)_. (accepted)
ChemRxiv, Preprint. 2021. [doi:10.33774/chemrxiv-2021-psqq4-v2](https://doi.org/10.33774/chemrxiv-2021-psqq4-v2)
- Kosugi T, Ohue M. **Development of a quantitative estimate index for early-stage screening of compounds targeting protein-protein interactions**. (under revision)
