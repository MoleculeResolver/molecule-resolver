from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from moleculeresolver import MoleculeResolver

if __name__ == "__main__":
    smiles_list = [
        "CS(C)=O",
        "O=S1CCCC1",
        "O=S1CCOCC1",
        "CCS(=O)CC",
        "CCCS(C)=O",
        "CCCS(=O)CCC",
        "O=S(c1ccccc1)c1ccccc1",
        "O=S1CCCC1",
        "CCCS(=O)CCC",
        "CCS(=O)CC",
        "CCCS(C)=O",
        "O=S1CCCC1",
        "O=S1CCOCC1"
    ]

    for smiles in smiles_list:
        print('-'*70)
        mr = MoleculeResolver()

        # Normalize molecule with rdkit
        mol = Chem.MolFromSmiles(smiles)
        normalized_mol = rdMolStandardize.Normalize(mol)
        normalized_smiles = Chem.MolToSmiles(normalized_mol)
        print(f"Original: {smiles} -> RDKit Normalized: {normalized_smiles}")
        
        # Convert zwitterionic form back to sulfynil group
        corrected_mol = mr.convert_zwitterion_to_sulfynil(normalized_mol)
        corrected_smiles = Chem.MolToSmiles(corrected_mol)
        print(f"Original: {smiles} -> Corrected       : {corrected_smiles}")

        # with MolResolver
        mr_smiles = mr.standardize_SMILES(smiles)
        print(f"Original: {smiles} -> MR smiles       : {mr_smiles}")

        assert smiles == mr_smiles == corrected_smiles
        

        




