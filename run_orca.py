'''Script to run all molecules with ORCA'''
import sys
import os
import csv
import pandas as pd
from ase.calculators.orca import ORCA, OrcaProfile
from organic_phos_script import optimize_geometry, st_gap_calculate

# READ INPUT
target_row = int(sys.argv[1])
intbl = pd.read_csv('smiles_energy.csv')
# id is in column mol_id, smiles is in column smiles
mol_ids = intbl['mol_id'].tolist()
smiles_strings = intbl['smiles'].tolist()
mol_id = mol_ids[target_row]
smile = smiles_strings[target_row]
geom_path = f'new_organic_phos_geometries/orca_{mol_id}.xyz'

# CONSTRUCT CALCULATORS
# Simple input line
simpleinput = '{mode} wB97M-V def2-TZVPPD def2/J RIJCOSX DefGrid3 TightSCF EnGrad'
singlet_simpleinput = simpleinput.format(mode='RKS')
triplet_simpleinput = simpleinput.format(mode='UKS')

# Blocks
blocks = '''%pal
  nprocs 24
end

%maxcore 3000

%scf
  MaxIter 150
end'''

# Orca directory is in the EBROOTORCA environment variable
orca_dir_path = os.environ.get('EBROOTORCA')
orca_exec_path = os.path.join(orca_dir_path, 'orca')
profile = OrcaProfile(command=orca_exec_path)

orca_calc_singlet = ORCA(
        profile=profile,
        orcasimpleinput=singlet_simpleinput,
        orcablocks=blocks,
        directory=f'orca_singlet_{mol_id}',
        label=f'orca_singlet_{mol_id}',
        charge=0, mult=1)
orca_calc_triplet = ORCA(
        profile=profile,
        orcasimpleinput=triplet_simpleinput,
        orcablocks=blocks,
        directory=f'orca_triplet_{mol_id}',
        label=f'orca_triplet_{mol_id}',
        charge=0, mult=3)

# RUN CALCULATION
atom = optimize_geometry(geom_path = geom_path,
        mol_id=mol_id, smile=smile, calc_triplet=orca_calc_triplet)
st_gap = st_gap_calculate(atom, mol_id=mol_id,
        calc_singlet=orca_calc_singlet, calc_triplet=orca_calc_triplet)
outrow = [mol_id, float(st_gap)]

# WRITE RESULTS
outfilename = 'orca_st_gap.csv'
if os.path.exists(outfilename):
    with open(outfilename, 'a') as f:
        writer = csv.writer(f)
        writer.writerow(outrow)
else:
    with open(outfilename, 'w') as f:
        header = ['mol_id', 'st_gap']
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerow(outrow)
