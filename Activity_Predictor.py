import numpy as np
import pickle
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


class Activity_Predictor:
    def __init__(self, model_path: str, radius: int, nBits: int):
        self.model = pickle.load(open(model_path, "rb"))
        self.radius = radius
        self.nBits = nBits

    def run_predictor(self, smiles):
        vector = []
        if isinstance(smiles, str):
            smiles = [smiles]
        for s in smiles:
            s = np.array(
                rdMolDescriptors.GetMorganFingerprintAsBitVect(
                    Chem.MolFromSmiles(s.strip()), radius=self.radius, nBits=self.nBits
                )
            )
            vector.append(s)
        return (
            self.model.predict(vector).tolist(),
            self.model.predict_proba(vector)[:, 1].tolist(),
        )
