import os
import ase.io
from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.orca import ORCA, OrcaProfile
from tblite.ase import TBLite
from rdkit import Chem
from rdkit.Chem import rdDistGeom

# TBLITE CALCULATOR
tblite_calc_singlet = TBLite(multiplicity=1)
tblite_calc_triplet = TBLite(multiplicity=3)

# ORCA CALCULATOR
# Simple input line
simpleinput = '{mode} wB97M-V def2-TZVPPD def2/J RIJCOSX DefGrid3 TightSCF EnGrad'
singlet_simpleinput = simpleinput.format(mode='RKS')
triplet_simpleinput = simpleinput.format(mode='UKS')

# Blocks
blocks = '''%pal
  nprocs 192
end

%maxcore 3000

%scf
  MaxIter 150
end'''

# Orca directory is in the EBROOTORCA environment variable
orca_dir_path = os.environ.get('EBROOTORCA')
orca_exec_path = os.path.join(orca_dir_path, 'orca')
profile = OrcaProfile(command=orca_exec_path)

orca_singlet_calc = ORCA(
        profile=profile,
        orcasimpleinput=singlet_simpleinput,
        orcablocks=blocks,
        charge=0, mult=1)
orca_triplet_calc = ORCA(
        profile=profile,
        orcasimpleinput=triplet_simpleinput,
        orcablocks=blocks,
        charge=0, mult=3)

geom_dir_name = 'new_organic_phos_geometries'
try:
    os.mkdir(geom_dir_name)
except FileExistsError:
    pass

test_mol_id = 'abp'
test_smile = "Nc1ccccc1C(=O)c1ccccc1"

def optimize_geometry(geom_path, mol_id, smile, calc_triplet):
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
        atom.calc = calc_triplet
        atom.calc.set_label(f'{mol_id}_opt')
        opt = BFGS(atom, logfile=f'{mol_id}.opt', trajectory=f'{mol_id}.traj')
        opt.run(fmax=0.05)
        ase.io.write(filename=geom_path, images=atom)
    return atom

def st_gap_calculate(atom, mol_id, calc_singlet, calc_triplet):
    diff_energy = 0
    for multiplicity in ['singlet', 'triplet']:
        if multiplicity == 'singlet':
            atom.calc = calc_singlet
            atom.calc.set_label(f'{mol_id}_singlet')
            energy = atom.get_potential_energy()
            diff_energy -= energy
        if multiplicity == 'triplet':
            atom.calc = calc_triplet
            atom.calc.set_label(f'{mol_id}_triplet')
            energy = atom.get_potential_energy()
            diff_energy += energy
    return(diff_energy)

# Test using abp.xyz
if __name__ == "__main__":
    geom_path = os.path.join(geom_dir_name, f'{test_mol_id}.xyz')
    atom = optimize_geometry(geom_path=geom_path, mol_id=test_mol_id, smile=test_smile, calc_triplet=tblite_calc_triplet)
    st_gap = st_gap_calculate(atom, mol_id=test_mol_id, calc_singlet=tblite_calc_singlet, calc_triplet=tblite_calc_triplet)
    print(st_gap)
