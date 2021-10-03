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

WEIGHT = {
    "params":
    {
        "MW":
            {
                "coefficients": 
                [0.04560031139224022, 1.7064244366404526, 494.04717332780626, 154.79350325269954, 77.10110879171509, 65.8193273039995]
            },
        "ALOGP":
            {
                "coefficients":
                [0.016181719609428855, 1.0652885460460373, 4.667743179259679, 5.065361777530944, 0.8433343264178254, 0.7385272147814145]
            },
        "HBD":
            {
                "coefficients":
                [0.01658054730902974, 1.1048232815489003, 1.6820912395701388, 2.7489254244785815, 0.4491737335600519, 0.5287174552710909]
            },
        "HBA":
            {
                "coefficients":
                [0.027008607260699185, 3.626899936363813, 4.311787382094175, 2.0241291464133954e-05, 0.7128127070334573, 1.2488426349682793]
            },
        "TPSA":
            {
                "coefficients":
                [-0.02291565292846544, 3.363039638274869, 61.95246699929325, 1.8933659551307073e-05, 11.596413833566858, 32.792265238065255]
            },
        "ROTB":
            {
                "coefficients":
                [0.10174790710720112, 0.9321790090551001, 5.795246533579942, 5.2517652130122645, 0.8612779744315607, 0.4733348902767623]
            },
        "AROM":
            {
                "coefficients":
                [-0.006806482079154438, 3.583035986979902, 2.392000709393099, 0.00020979756196675626, 0.47348016464327786, 1.0012588955406478]
        },
    },
    "weights": 
    {
        "max": [0.5, 0.0, 1.0, 1.0, 0.25, 0.5, 1.0],
        "max_ents": 380.89604800900497,
        "mo": [0.47025, 0.10025, 0.819, 0.80725, 0.36725, 0.5295, 0.89025],
        "mo_ents": 378.69627297573555
    }
}
