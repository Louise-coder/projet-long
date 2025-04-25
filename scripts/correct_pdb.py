from Bio.PDB import PDBParser, PDBIO, StructureBuilder, Select
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue

class NoHydrogens(Select):
    def accept_atom(self, atom):
        return atom.element != "H"

def correct_pdb(pdb_file: str, inf: int, sup: int, new_inf: int):
    parser = PDBParser(QUIET=True)
    original = parser.get_structure("prot", pdb_file)

    new_structure = Structure("new")
    new_model = Model(0)
    new_chain = Chain("A")

    p = new_inf
    for model in original:
        for chain in model:
            for residue in chain:
                res_id = residue.id[1]
                if inf <= res_id <= sup:
                    new_residue = Residue((' ', p, ' '), residue.resname, residue.segid)
                    for atom in residue:
                        if atom.element != "H":
                            new_residue.add(atom.copy())
                    new_chain.add(new_residue)
                    p += 1

    new_model.add(new_chain)
    new_structure.add(new_model)

    io = PDBIO()
    io.set_structure(new_structure)
    io.save("structure_modifiée.pdb")

if __name__ == '__main__':
    PDB_FILE = "aligned-pred.pdb"
    correct_pdb(PDB_FILE, 261, 316, 283)