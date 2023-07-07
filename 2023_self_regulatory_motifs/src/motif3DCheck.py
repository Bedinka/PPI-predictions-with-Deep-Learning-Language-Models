#!/usr/bin/env python3

import numpy as np
import os, csv
from Bio.PDB import PDBParser, Select
from Bio.PDB.ccealign import run_cealign
from Bio.PDB.PDBExceptions import PDBException
from Bio.PDB.qcprot import QCPSuperimposer
import Bio.PDB.PDBIO as io


_RESID_SORTER = lambda r: r.id[1]  # noqa: E731

# CustomSelect class for selective output
class CustomSelect(Select):
    def __init__(self, start_index, end_index):
        self.start_index = start_index
        self.end_index = end_index

    def accept_residue(self, residue):
        residue_index = residue.get_id()[1]
        if self.start_index <= residue_index <= self.end_index:
            return 1  # Accept the residue
        else:
            return 0  # Reject the residue

class CEAligner2:
    """Protein Structure Alignment by Combinatorial Extension."""

    def __init__(self, window_size=8, max_gap=30):
        """Superimpose one set of atoms onto another using structural data.

        Structures are superimposed using guide atoms, CA and C4', for protein
        and nucleic acid molecules respectively.

        Parameters
        ----------
        window_size : float, optional
            CE algorithm parameter. Used to define paths when building the
            CE similarity matrix. Default is 8.
        max_gap : float, optional
            CE algorithm parameter. Maximum gap size. Default is 30.
        """
        assert window_size > 0, "window_size must be greater than 0"
        assert max_gap >= 0, "max_gap must be positive (or zero)"

        self.window_size = window_size
        self.max_gap = max_gap

        self.rms = None

    def get_guide_coord_from_structure(self, structure):
        """Return the coordinates of guide atoms in the structure.

        We use guide atoms (C-alpha and C4' atoms) since it is much faster than
        using all atoms in the calculation without a significant loss in
        accuracy.
        """
        coords = []
        # CE algorithm is sensitive to atom ordering. To reproduce Pymol
        # results, sort atoms by chain and then residue number.
        for chain in sorted(structure.get_chains()):
            for resid in sorted(chain, key=_RESID_SORTER):
                if "CA" in resid:
                    coords.append(resid["CA"].coord.tolist())
                elif "C4'" in resid:
                    coords.append(resid["C4'"].coord.tolist())
        if not coords:
            msg = f"Structure {structure.id} does not have any guide atoms."
            raise PDBException(msg)
        return coords

    def set_reference(self, structure):
        """Define a reference structure onto which all others will be aligned."""
        self.refcoord = self.get_guide_coord_from_structure(structure)

        if len(self.refcoord) < self.window_size * 2:
            n_atoms = len(self.refcoord)
            msg = (
                f"Too few atoms in the reference structure ({n_atoms}). "
                "Try reducing the window_size parameter."
            )
            raise PDBException(msg)

    def align(self, structure, transform=True):
        """Align the input structure onto the reference structure.

        Parameters
        ----------
        transform: bool, optional
            If True (default), apply the rotation/translation that minimizes
            the RMSD between the two structures to the input structure. If
            False, the structure is not modified but the optimal RMSD will
            still be calculated.
        """
        self.rms = None  # clear before aligning

        coord = self.get_guide_coord_from_structure(structure)

        if len(coord) < self.window_size * 2:
            n_atoms = len(coord)
            msg = (
                f"Too few atoms in the mobile structure ({n_atoms}). "
                "Try reducing the window_size parameter."
            )
            raise PDBException(msg)

        # Run CEAlign
        # CEAlign returns the best N paths, where each path is a pair of lists
        # with aligned atom indices. Paths are not guaranteed to be unique.
        paths = run_cealign(self.refcoord, coord, self.window_size, self.max_gap)
        unique_paths = {(tuple(pA), tuple(pB)) for pA, pB in paths}

        # Iterate over unique paths and find the one that gives the lowest
        # corresponding RMSD. Use QCP to align the molecules.
        best_rmsd, best_u = 1e6, None
        for u_path in unique_paths:
            idxA, idxB = u_path

            coordsA = np.array([self.refcoord[i] for i in idxA])
            coordsB = np.array([coord[i] for i in idxB])

            aln = QCPSuperimposer()
            aln.set(coordsA, coordsB)
            aln.run()
            if aln.rms < best_rmsd:
                best_rmsd = aln.rms
                best_u = (aln.rot, aln.tran)

        if best_u is None:
            raise RuntimeError("Failed to find a suitable alignment.")

        if transform:
            # Transform all atoms
            rotmtx, trvec = best_u
            for chain in structure.get_chains():
                for resid in chain.get_unpacked_list():
                    for atom in resid.get_unpacked_list():
                        atom.transform(rotmtx, trvec)
        self.best_u = best_u
        self.rms = best_rmsd
        
    def align_from_best_u(self, structure, transform=True):
        if transform:
            # Transform all atoms
            rotmtx, trvec = self.best_u
            for chain in structure.get_chains():
                for resid in chain.get_unpacked_list():
                    for atom in resid.get_unpacked_list():
                        atom.transform(rotmtx, trvec)
        
trgDomain = "/home/roger/2023_self_regulatory_motifs/SwissProt/domainStructures/PF18883.3/best_model.pdb"

# pdb2 = "/home/roger/2023_self_regulatory_motifs/SwissProt/domainStructures/PF18883.3/AF-Q9JMS3-F1-model_v4_903-1009_PF18883.3.pdb"
# pdb3 = "/home/roger/2023_self_regulatory_motifs/SwissProt/data/AF-Q9JMS3-F1-model_v4.pdb"
pdb2 = "/home/roger/2023_self_regulatory_motifs/SwissProt/domainStructures/PF18883.3/AF-P39180-F1-model_v4_573-690_PF18883.3.pdb"
pdb3 = "/home/roger/2023_self_regulatory_motifs/SwissProt/data/AF-P39180-F1-model_v4.pdb"


def AlignPDBWithDomain( srcDomain, srcPDB, trgDomain ):
    reference_structure = PDBParser().get_structure( "reference", trgDomain )
    structure = PDBParser().get_structure("structure", srcDomain )
    structure2 = PDBParser().get_structure("structure", srcPDB )
    aligner1 = CEAligner2()

    aligner1.set_reference(reference_structure)
    aligner1.align(structure) 
    aligner1.align_from_best_u(structure2)
    
    # io.set_structure(reference_structure)
    # io.save("/home/roger/2023_self_regulatory_motifs/trgDomain2.pdb")
    # io.set_structure(structure)
    # io.save("/home/roger/2023_self_regulatory_motifs/srcDomain2.pdb")
    # io.set_structure(structure2)
    # io.save("/home/roger/2023_self_regulatory_motifs/srcPDB2.pdb")
    return structure2
# AlignPDBWithDomain( pdb2, pdb3, trgDomain )



motifData = dict()
table1 = open('SwissProt/SwissProt_results.csv', 'r')
reader1 = csv.reader(table1)
next(reader1)
for line in reader1:
    Species, PDB, Domain, domainName, motifSeq, domainStrart, domainEnd, motifStart, motifEnd, distanceDM = line
    domainPos = domainEnd+'-'+domainStrart
    if Domain not in motifData.keys():
        motifData[Domain] = dict()
    if PDB not in motifData[Domain].keys():
        motifData[Domain][PDB] = dict()
    if domainPos not in motifData[Domain][PDB].keys():
        motifData[Domain][PDB][domainPos] = [(motifStart, motifEnd)]
    else:
        motifData[Domain][PDB][domainPos].append((motifStart, motifEnd))


pdbdir = '/home/roger/2023_self_regulatory_motifs/SwissProt'
c = 0
# Process directories one by one
for entry in os.scandir(os.path.join(pdbdir, 'domainStructures')):
    if entry.is_dir():
        domain = entry.path.split('/')[-1]
        for section in os.scandir(entry.path):
            s = section.path.split('/')[-1]
            if s=='best_model.pdb': continue
            protein = s.split('v4_')[0]+'v4'
            domainPos = s.split(protein)[-1].split('_')[1]
            if (domain not in motifData) or (protein not in motifData[domain]) or (domainPos not in motifData[domain][protein]):
                continue
            motifStart, motifEnd = motifData[domain][protein][domainPos]
            print(motifStart, motifEnd)
            print(motifData[domain][protein][domainPos])
            
            srcDomain = section.path
            srcPDB = pdbdir+'/data/'+protein+'.pdb'
            trgDomain = os.path.join(entry.path, "best_model.pdb")
            wholeProt = AlignPDBWithDomain(srcDomain, srcPDB, trgDomain)
            
            # select = CustomSelect(start_index=motifStart, end_index=motifEnd)
            
        # continue
        # # Load the reference structure
        # reference_structure = PDBParser().get_structure("reference", os.path.join(entry.path, "best_model.pdb"))
        # # for loop that saves the positions of the motifs after alignment
        # for idx, motif in enumerate(motifDict[entry.name]):
        #     print(motif)
        #     # AF-Q9ZW31-F1-model_v4_49-152_PF00011.24.pdb
        #     structure_entry= motif[1] +'_'+ motif[5] +'-'+ motif[6] +'_'+ motif[2] + '.pdb'
        #     print(structure_entry)
        #     # Load the structure to align
        #     structure = PDBParser().get_structure("structure", os.path.join(entry.path, structure_entry))
            
        #     # Align structure to the reference
        #     aligner1 = CEAligner()
        #     aligner1.set_reference(reference_structure)
        #     aligner1.align(structure)
        #     # Save the translation and rotation (code line is missing, add appropriate code here)
        #     moves = aligner1.get_guide_coord_from_structure(structure)
        #     print(moves)
        #     # Align full structure to moved structure
        #     # full_structure = PDBParser().get_structure("full_structure", os.path.join('/home/roger/2023_self_regulatory_motifs/SwissProt/data/', structure_entry.split('_')[0]+'_'+structure_entry.split('_')[1]+'.pdb'))
        #     # aligner2 = CEAligner()
        #     # aligner2.set_reference(structure)
        #     # aligner2.align(full_structure)
        #     # should_be_moves = aligner2.refcoord
        #     # assert(should_be_moves==moves)
        #     # motifStart = int(motif[7])
        #     # motifEnd = int(motif[8])
        #     # motifPosition = aligner2.get_guide_coord_from_structure(full_structure)[motifStart:motifEnd]
        #     # motifDict[entry.name][idx].append(motifPosition)
        #     print('\n----------------------------')
        # break
