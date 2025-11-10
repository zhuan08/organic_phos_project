'''Script to run all molecules with ORCA'''
import csv
import pandas as pd
from fairchem.core import pretrained_mlip, FAIRChemCalculator
from organic_phos_script import optimize_geometry, st_gap_calculate

MODEL = 'uma-m-1p1'  # Models: uma-s-1p1, uma-m-1p1
DEVICE = 'cpu'  # or 'cpu'

# Load the pretrained model
predictor = pretrained_mlip.get_predict_unit(MODEL, device=DEVICE) # Models: uma-s-1p1, uma-m-1p1

def uma_singlet_setup(atoms):
    atoms.info['charge'] = 0
    atoms.info['spin'] = 1
    atoms.calc = FAIRChemCalculator(predictor, task_name='omol')

def uma_triplet_setup(atoms):
    atoms.info['charge'] = 0
    atoms.info['spin'] = 3
    atoms.calc = FAIRChemCalculator(predictor, task_name='omol')

# READ INPUT
intbl = pd.read_csv('smiles_energy.csv')
# id is in column mol_id, smiles is in column smiles
mol_ids = intbl['mol_id'].tolist()
smiles_strings = intbl['smiles'].tolist()

# RUN CALCULATION
outrows = []
for mol_id, smile in zip(mol_ids, smiles_strings):
    geom_path = f'new_organic_phos_geometries/uma_{mol_id}.xyz'
    atom = optimize_geometry(geom_path = geom_path,
            mol_id=mol_id, smile=smile, triplet_setup=uma_triplet_setup)
    st_gap = st_gap_calculate(atom,
            singlet_setup=uma_singlet_setup, triplet_setup=uma_triplet_setup)
    outrow = [mol_id, float(st_gap)]
    outrows.append(outrow)

# WRITE RESULTS
outfilename = 'uma_st_gap.csv'
with open(outfilename, 'w') as f:
    header = ['mol_id', 'st_gap']
    writer = csv.writer(f)
    writer.writerow(header)
    for outrow in outrows:
        writer.writerow(outrow)
