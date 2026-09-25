"""A tiny synthetic mmCIF exercising every rule family, written as literal text.

Chain A, entity 1 (8 residues): G1 H2 MSE3 H4 G5 G6 G7 H8; residues 1 (N-terminal) and 5
(internal) unobserved. Zn coordinated tetrahedrally by His2/His4/His8 NE2 and a water.
A sulfate (no rule) and a water with two altlocs complete it.
"""

from simprep.structure.geometry import Point
from tests.unit.builders import TETRAHEDRAL, scaled

ZN_LIGAND_ANGSTROM = 2.05


def _positions() -> dict[str, Point]:
    t = [scaled(d, ZN_LIGAND_ANGSTROM) for d in TETRAHEDRAL]
    return {
        "H2_NE2": t[0],
        "H4_NE2": t[1],
        "H8_NE2": t[2],
        "W_O": t[3],
        "H2_CA": scaled(TETRAHEDRAL[0], 6.0),
        "M3_CA": (8.0, 0.0, 0.0),
        "M3_SE": (9.5, 0.0, 0.0),
        "H4_CA": scaled(TETRAHEDRAL[1], 6.0),
        "G6_CA": (0.0, -7.0, -2.0),
        "G7_CA": (0.0, -9.0, 1.0),
        "H8_CA": scaled(TETRAHEDRAL[2], 6.0),
        "SO4_S": (20.0, 20.0, 20.0),
        "W2_A": (15.0, 0.0, 0.0),
        "W2_B": (15.6, 0.0, 0.0),
    }


ATOM_ROWS = [
    # group, id, element, atom, alt, comp, asym, entity, label_seq, key, occ, auth_seq
    ("ATOM", "C", "CA", ".", "HIS", "A", 1, 2, "H2_CA", 1.0, 2),
    ("ATOM", "N", "NE2", ".", "HIS", "A", 1, 2, "H2_NE2", 1.0, 2),
    ("HETATM", "C", "CA", ".", "MSE", "A", 1, 3, "M3_CA", 1.0, 3),
    ("HETATM", "SE", "SE", ".", "MSE", "A", 1, 3, "M3_SE", 1.0, 3),
    ("ATOM", "C", "CA", ".", "HIS", "A", 1, 4, "H4_CA", 1.0, 4),
    ("ATOM", "N", "NE2", ".", "HIS", "A", 1, 4, "H4_NE2", 1.0, 4),
    ("ATOM", "C", "CA", ".", "GLY", "A", 1, 6, "G6_CA", 1.0, 6),
    ("ATOM", "C", "CA", ".", "GLY", "A", 1, 7, "G7_CA", 1.0, 7),
    ("ATOM", "C", "CA", ".", "HIS", "A", 1, 8, "H8_CA", 1.0, 8),
    ("ATOM", "N", "NE2", ".", "HIS", "A", 1, 8, "H8_NE2", 1.0, 8),
    ("HETATM", "ZN", "ZN", ".", "ZN", "B", 2, None, None, 1.0, 901),
    ("HETATM", "S", "S", ".", "SO4", "C", 3, None, "SO4_S", 1.0, 902),
    ("HETATM", "O", "O", ".", "HOH", "D", 4, None, "W_O", 1.0, 1001),
    ("HETATM", "O", "O", "A", "HOH", "D", 4, None, "W2_A", 0.6, 1002),
    ("HETATM", "O", "O", "B", "HOH", "D", 4, None, "W2_B", 0.4, 1002),
]

HEADER = """data_MINI
_entry.id MINI
#
loop_
_entity.id
_entity.type
1 polymer
2 non-polymer
3 non-polymer
4 water
#
loop_
_entity_poly_seq.entity_id
_entity_poly_seq.num
_entity_poly_seq.mon_id
_entity_poly_seq.hetero
1 1 GLY n
1 2 HIS n
1 3 MSE n
1 4 HIS n
1 5 GLY n
1 6 GLY n
1 7 GLY n
1 8 HIS n
#
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_auth_asym_id
_struct_conn.ptnr1_auth_seq_id
_struct_conn.pdbx_ptnr1_PDB_ins_code
_struct_conn.ptnr1_auth_comp_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.pdbx_ptnr1_label_alt_id
_struct_conn.ptnr2_auth_asym_id
_struct_conn.ptnr2_auth_seq_id
_struct_conn.pdbx_ptnr2_PDB_ins_code
_struct_conn.ptnr2_auth_comp_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.pdbx_ptnr2_label_alt_id
_struct_conn.pdbx_dist_value
metalc1 metalc A 901 ? ZN ZN ? A 2 ? HIS NE2 ? 2.05
#
loop_
_pdbx_struct_mod_residue.id
_pdbx_struct_mod_residue.auth_asym_id
_pdbx_struct_mod_residue.auth_seq_id
_pdbx_struct_mod_residue.PDB_ins_code
_pdbx_struct_mod_residue.auth_comp_id
_pdbx_struct_mod_residue.label_comp_id
_pdbx_struct_mod_residue.parent_comp_id
_pdbx_struct_mod_residue.details
1 A 3 ? MSE MSE MET SELENOMETHIONINE
#
loop_
_pdbx_unobs_or_zero_occ_residues.id
_pdbx_unobs_or_zero_occ_residues.PDB_model_num
_pdbx_unobs_or_zero_occ_residues.polymer_flag
_pdbx_unobs_or_zero_occ_residues.occupancy_flag
_pdbx_unobs_or_zero_occ_residues.auth_asym_id
_pdbx_unobs_or_zero_occ_residues.auth_comp_id
_pdbx_unobs_or_zero_occ_residues.auth_seq_id
_pdbx_unobs_or_zero_occ_residues.PDB_ins_code
_pdbx_unobs_or_zero_occ_residues.label_asym_id
_pdbx_unobs_or_zero_occ_residues.label_comp_id
_pdbx_unobs_or_zero_occ_residues.label_seq_id
1 1 Y 1 A GLY 1 ? A GLY 1
2 1 Y 1 A GLY 5 ? A GLY 5
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
"""


def mini_cif_text(models: int = 1) -> str:
    positions = _positions() | {None: (0.0, 0.0, 0.0)}
    rows = []
    serial = 0
    for model in range(1, models + 1):
        for group, element, atom, alt, comp, asym, entity, seq, key, occ, auth in ATOM_ROWS:
            serial += 1
            x, y, z = positions[key]
            rows.append(
                f"{group} {serial} {element} {atom} {alt} {comp} {asym} {entity} "
                f"{seq if seq is not None else '.'} ? {x:.3f} {y:.3f} {z:.3f} {occ} 20.00 "
                f"{auth} A {model}"
            )
    return HEADER + "\n".join(rows) + "\n#\n"
