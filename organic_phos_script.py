import os
import ase.io
from ase import Atoms
from ase.optimize import BFGS
from tblite.ase import TBLite
from ase.calculators.gamess_us import GAMESSUS
from ase.calculators.psi4 import Psi4
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDistGeom
import pandas as pd
IPythonConsole.ipython_3d = True

CHARGE = 0
METHOD = 'b3lyp'
BASIS = 'def2-svp'
psi4_calc_singlet = Psi4(method=METHOD, basis=BASIS, charge=CHARGE, multiplicity=1)
psi4_calc_triplet = Psi4(method=METHOD, basis=BASIS, charge=CHARGE,
                  multiplicity=3, reference='uks')
tblite_calc_singlet = TBLite(multiplicity=1)
tblite_calc_triplet = TBLite(multiplicity=3)
gamess_calc_singlet = GAMESSUS(contrl={'mult': 1}, label='molecule',
                    command='rungms PREFIX.inp 30Jun2020R1 > PREFIX.log 2> PREFIX.err')
gamess_calc_triplet = GAMESSUS(contrl={'mult': 3}, label='molecule',
                    command='rungms PREFIX.inp 30Jun2020R1 > PREFIX.log 2> PREFIX.err')

geom_dir_name = 'new_organic_phos_geometries'
try:
    os.mkdir(geom_dir_name)
except FileExistsError:
    pass

m_id = 'abp'
abp_smile = "Nc1ccccc1C(=O)c1ccccc1"

def optimize_geometry(mol_id, smile, calc_tiplet):
    geom_path = os.path.join(geom_dir_name, f'{mol_id}.xyz')
    # First, check if there's already a geometry saved, and if so, just load it
    if os.path.exists(geom_path):
        atom = ase.io.read(filename=geom_path)
    else:
        molecule = Chem.MolFromSmiles(smile)
        if molecule is None: 
            raise ValueError(f'MolFromSmiles returned None on {mol_id}')
        molecule = Chem.AddHs(molecule)
        rdDistGeom.EmbedMolecule(molecule)
        conf_mol = molecule.GetConformer()
        pos_mol = conf_mol.GetPositions()
        atom_mol = molecule.GetAtoms()
        atom_sym = []
        for char in atom_mol:
            atom_sym.append(char.GetSymbol())
        atom = Atoms(atom_sym, positions=pos_mol)
        atom.calc = calc_tiplet
        opt = BFGS(atom, logfile=None, trajectory=None)
        opt.run(fmax=0.05)
        ase.io.write(filename=geom_path, images=atom)

def st_gap_calculate(path, calc_singlet, calc_triplet):
    atom = ase.io.read(filename=path)
    diff_energy = 0
    for multiplicity in ['singlet', 'triplet']:
        if multiplicity == 'singlet':
            atom.calc = calc_singlet
            energy = atom.get_potential_energy()
            diff_energy -= energy
        if multiplicity == 'triplet':
            atom.calc = calc_triplet
            energy = atom.get_potential_energy()
            diff_energy += energy
    print(diff_energy)

# Test using abp.xyz
optimize_geometry(mol_id=m_id, smile=abp_smile, calc_tiplet=tblite_calc_triplet)
st_gap_calculate(path="new_organic_phos_geometries/abp.xyz", calc_singlet=tblite_calc_singlet, calc_triplet=tblite_calc_triplet)
