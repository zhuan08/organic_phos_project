'''Script to run all molecules in parallel'''
import csv
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd
from organic_phos_script import tblite_calc_singlet, tblite_calc_triplet, optimize_geometry, st_gap_calculate

max_workers = 20
intbl = pd.read_csv('smiles_energy.csv')
# id is in column mol_id, smiles is in column smiles
mol_ids = intbl['mol_id'].tolist()
smiles_strings = intbl['smiles'].tolist()
def run(calc_name, calc_singlet, calc_triplet, mol_id, smile):
    try:
        geom_path = f'new_organic_phos_geometries/{calc_name}_{mol_id}.xyz'
        optimize_geometry(geom_path = geom_path,
                mol_id=mol_id, smile=smile, calc_triplet=calc_triplet)
        st_gap = st_gap_calculate(path=geom_path,
                calc_singlet=calc_singlet, calc_triplet=calc_triplet)
        outrow = {'mol_id': mol_id, 'smiles': smile,
                  'calc_name': calc_name, 'st_gap': st_gap,
                  'error': None}
    except Exception as e:
        outrow = {'mol_id': mol_id, 'smiles': smile,
                  'calc_name': calc_name, 'st_gap': None,
                  'error': str(e)}
    return outrow
tasks = [('tblite', tblite_calc_singlet, tblite_calc_triplet, mol_id, smiles) \
        for mol_id, smiles in zip(mol_ids, smiles_strings)]
outrows = []
with ProcessPoolExecutor(max_workers=max_workers) as executor:
    futures = [executor.submit(run, *task) for task in tasks]
    for future in as_completed(futures):
        outrow = future.result()
        outrows.append(outrow)

pd.DataFrame(outrows).to_csv('st_gaps.csv', index=False)
