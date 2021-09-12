import pytest
import qeppi as ppi
from rdkit.Chem import SDMolSupplier


def test_qeppi_value_approx():
    correct_values = [
        0.14138488397092042,
        0.20844419255972127,
        0.20349457838077856,
        0.8374457753324647,
        0.40286064451685927,
        0.1037465060768894,
        0.08835509613475574,
        0.5713916001329409,
        0.7897952756770426,
        0.5757902101651236,
        0.759799221073293,
        0.5237059407202563,
        0.5713814680494992,
        0.6021230144891282,
    ]

    q = ppi.QEPPI_Calculator()
    q.load("./model/QEPPI.model")

    ppi_s = SDMolSupplier("./compunds/test.sdf")
    ppi_mols = [mol for mol in ppi_s if mol is not None]
    for i, j in zip(ppi_mols, correct_values):
        assert q.qeppi(i) == pytest.approx(j)