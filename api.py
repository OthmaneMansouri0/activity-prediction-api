# -*- coding: utf-8 -*-
"""
Activity Predictor
"""

from flask import Flask, request
from rdkit import Chem

from Activity_Predictor import Activity_Predictor

ACTIVITY_PREDICTOR = Activity_Predictor(model_path="model.pkl", radius=4, nBits=512)

app = Flask(__name__)


@app.route("/score", methods=["POST", "GET"])
def score():
    if request.method == "POST":
        smiles = request.get_json().get(list(request.get_json().keys())[0])
        if not (
            isinstance(smiles, str)
            or (
                isinstance(smiles, list)
                and set([type(item) for item in smiles]) == {str}
            )
        ):
            return {"Error": "Smiles must be a string or a list of strings"}
        else:
            # check if molecules exist
            wrong_smiles = []
            if isinstance(smiles, str):
                smiles = [smiles]
            for s in smiles:
                if Chem.MolFromSmiles(s.strip()) is None:
                    wrong_smiles.append(s)
            if wrong_smiles == []:
                return {
                    "scores": [
                        dict(smiles=r[0], label=r[1], probability=r[2])
                        for r in zip(
                            smiles,
                            ACTIVITY_PREDICTOR.run_predictor(smiles)[0],
                            ACTIVITY_PREDICTOR.run_predictor(smiles)[1],
                        )
                    ]
                }
            else:
                return {"Error": f"Wrong smiles : {', '.join(wrong_smiles)}"}


if __name__ == "__main__":
    app.run()
