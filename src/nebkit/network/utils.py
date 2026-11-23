from rdkit import Chem 
from rdkit.Chem import AllChem
from rdkit import rdBase
from ase import Atoms  
from ase.neighborlist import NeighborList, natural_cutoffs
import numpy as np 
from copy import deepcopy

def mol_to_ase_atoms(mol:Chem.Mol):
    new_mol = deepcopy(mol)
    with rdBase.BlockLogs():
        new_mol.UpdatePropertyCache(strict=False)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42 
        AllChem.EmbedMolecule(new_mol, params)
    atomic_numbers = [atom.GetAtomicNum() for atom in new_mol.GetAtoms()] 
    conf = new_mol.GetConformer()
    positions = conf.GetPositions()
    atoms = Atoms(numbers=atomic_numbers, positions=positions)
    return atoms 

def get_adsorbate_binding_atom_index(atoms:Atoms) -> tuple | None:
    max_coord = {
        1:(1, 0), #H, 最大配位数1，0对孤对电子
        6:(4, 0), #C

        7:(3, 1), #N，最大配位数3，1对孤对电子
        8:(2, 2), #O
        9:(1, 0), #F

        15:(3, 1), #P
        16:(2, 2), #S
        17: (1, 0), #Cl
    }

    cutoffs = natural_cutoffs(atoms)
    nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(atoms)
    numbers = atoms.get_atomic_numbers()

    scores = []
    for i in range(len(atoms)):
        indices, offsets = nl.get_neighbors(i)
        remain_coord = max_coord[numbers[i]][0] - len(indices)
        score = remain_coord if remain_coord > 0 else 0
        if remain_coord < 0:
            score = score + max_coord[numbers[i]][1] + remain_coord * 0.5 
        else:
            score = score + max_coord[numbers[i]][1]
        scores.append(score)
    
    scores = np.array(scores)
    index = np.argsort(scores)

    if scores[index[-1]] == 0: #最大的不饱和度是0，没有键合原子
        return None 
    if len(index) == 1: # 只有一个原子，直接返回 0
        return (0, )
    if scores[index[-2]] == 0: #有一个以上的原子，但第二大的不饱和度是0 
        return (index[-1], )
    return (index[-1], index[-2])
    
        


def pattern_to_mol(symbols:str, edges:str):
    symbols_list = symbols.split(",")
    edges_list = edges.split(",")
    mol = Chem.RWMol()
    for s in symbols_list:
        mol.AddAtom(Chem.Atom(s))
    if len(symbols_list) == 1:
        return mol 
    
    for edge in edges_list:
        a,b = edge.split('-')
        mol.AddBond(int(a), int(b), Chem.BondType.SINGLE)
    return mol.GetMol() 

def ase_atoms_to_mol(atoms:Atoms, skin:float=0.1, ignoreH:bool = False) -> Chem.Mol:
    '''
    Convert `ase.Atoms` to `Chem.Mol`. All the bonds are set to be single bond.
    Args:
        atoms: `ase.Atoms` object.
        skin: For atom `A` and `B`, if their distance satisfy $d(A-B) < skin*2 + r(A) + r(B)$, consider them to be bonded. Here r(A) and r(B) are the covalent radii.
        ignoreH: ignore all H atoms.
    '''
    if ignoreH:
        atoms = atoms.copy()
        numbers = atoms.get_atomic_numbers()
        flag = numbers != 1
        atoms = Atoms(numbers=numbers[flag], positions = atoms.positions[flag])

    nl = NeighborList(natural_cutoffs(atoms), self_interaction=False, bothways=True, skin=skin)
    nl.update(atoms)
    natoms = len(atoms)
    mol = Chem.RWMol()
    symbols = atoms.get_chemical_symbols()
    for s in symbols:
        atom = Chem.Atom(s)
        atom.SetNoImplicit(True)
        mol.AddAtom(atom)
    for i in range(natoms):
        index, _ = nl.get_neighbors(i)
        for j in index:
            if j > i:
                mol.AddBond(i, int(j), Chem.BondType.SINGLE)
    with rdBase.BlockLogs():
        mol.UpdatePropertyCache(strict=False)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42 
        AllChem.EmbedMolecule(mol, params)
    return mol.GetMol()

