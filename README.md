# Activity_Predictor

This API takes as an input the Smiles of one or N molecules, and computes their activity using random forest classifier. 

### Install
pipenv install

### Run
pipenv run python api.py

### Input
A JSON with either :
- one smiles : {"smiles": "C"}
- N smiles : {"smiles" : ["C", "O"]}
