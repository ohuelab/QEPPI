import argparse

import QEPPI as ppi
from rdkit import Chem
from rdkit.Chem import SDMolSupplier
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sdf", help="Path to SDF file", type=str)
    parser.add_argument(
        # "The column name "SMILES" is required."
        "--csv",
        help="Path to CSV file",
        type=str,
    )
    parser.add_argument(
        "--smiles",
        help="Input valid SMILES",
        type=str,
    )
    parser.add_argument(
        "--out",
        help="Path to output CSV file",
        type=str,
        default="./result.csv",
    )
    args = parser.parse_args()

    q = ppi.QEPPI_Calculator()
    print("QEPPI.model LOADING...")
    q.load("./model/QEPPI.model")
    if not args.smiles:
        if args.sdf and args.csv:
            print("You need set either a SDF or a CSV file.")
            exit()
        elif args.sdf:
            print("SDF LOADING...")
            ppi_s = SDMolSupplier(args.sdf)
            ppi_mols = [mol for mol in ppi_s if mol is not None]
        elif args.csv:
            print("CSV LOADING...")
            smiles_df = pd.read_csv(args.csv)
            ppi_s = [Chem.MolFromSmiles(smiles) for smiles in smiles_df["SMILES"]]
            ppi_mols = [mol for mol in ppi_s if mol is not None]
        else:
            print("You need to set a SMILES or a SDF file or a CSV file")
            exit()
        print("Caculating QEPPI...")
        ppi_qeppi = list(map(q.qeppi, ppi_mols))
        smiles = list(map(Chem.MolToSmiles, ppi_mols))
        df = pd.DataFrame(list(zip(smiles, ppi_qeppi)), columns=["SMILES", "QEPPI"])
        print("Print head 5 rows")
        print(df.head())
        df.to_csv(args.out, index=False)
    if args.smiles:
        smiles = args.smiles
        mol = Chem.MolFromSmiles(smiles)
        if not mol == None:
            print(f"{smiles} : {q.qeppi(mol)}")
        else:
            print(
                f"{smiles} is not valid, because it can't be converted to MOL. Please confirm it"
            )


if __name__ == "__main__":
    main()