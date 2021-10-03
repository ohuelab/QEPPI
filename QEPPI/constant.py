from rdkit.Chem import rdMolDescriptors, Crippen
from rdkit import Chem

CHEM_LABELS = ["MW", "ALOGP", "HBD", "HBA", "TPSA", "ROTB", "AROM"]

DESC_FUNC_DICT = {
    "MW": rdMolDescriptors.CalcExactMolWt,
    "ALOGP": Crippen.MolLogP,
    "HBD": rdMolDescriptors.CalcNumHBD,
    "HBA": rdMolDescriptors.CalcNumHBA,
    "TPSA": rdMolDescriptors.CalcTPSA,
    "ROTB": rdMolDescriptors.CalcNumRotatableBonds,
    "AROM": rdMolDescriptors.CalcNumAromaticRings,
}