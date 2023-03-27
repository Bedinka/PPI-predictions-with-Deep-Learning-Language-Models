#!/usr/bin/env python3

import os
import numpy
import math
import networkx as nx
from pyvis.network import Network
import shelve
import gzip
from IPython.core.display import display, HTML
import matplotlib.pyplot as plt


# PDB READER
# https://cupnet.net/pdb-format/

LIGANDS = [ "NDP", "SPM" ]

class Atom:
    def __init__( self, line, atomname, pos, x, y, z ):
        self.line = line
        self.atomname = atomname
        self.pos = pos
        self.x = x
        self.y = y
        self.z = z

    def __str__( self ):
        return "%d,%s,%f,%f,%f" % ( self.pos, self.atomname, self.x, self.y, self.z )

class Residue:
    def __init__( self, pos, resname ):
        self.pos = pos
        self.atoms = []
        self.center = None
        self.CA = None
        self.resname = resname
    
    def __str__(self):
        return self.resname+str(self.pos)

    def get_center( self ):
        if self.center == None:
            x = numpy.average( [ atom.x for atom in self.atoms ] )
            y = numpy.average( [ atom.y for atom in self.atoms ] )
            z = numpy.average( [ atom.z for atom in self.atoms ] )
            self.center = Atom( "CENTER", "X", -1, x, y, z )
        return self.center

    def get_CA( self ):
        if self.CA == None:
            CAs = [ atom for atom in self.atoms if atom.atomname == "CA"]
            self.CA = CAs[0]
        return self.CA
     

class Chain:
    def __init__( self, chain_id ):
        self.chain_id = chain_id
        self.residues = {}
        self.residues_list = []

    def is_contacting( self, chain, dist_cutoff = 6.0 ):
        for atom1 in self.atoms:
            for atom2 in chain.atoms:
                dist = distance( atom1, atom2 )
                #print atom1, atom2, dist
                if dist <= dist_cutoff:
                    return True
        return False

    def all_to_all_distance( self, chain ):
        output = []
        for atom1 in self.atoms:
            for atom2 in chain.atoms:
                dist = distance( atom1, atom2 )
                output.append( [ dist, atom1, atom2 ] )
        return output

    def all_to_all_CA_distance( self ):
        output = []
        for residue1 in self.residues_list:
            for residue2 in self.residues_list:
                dist = distance( residue1.get_CA(), residue2.get_CA() )
                output.append( [ dist, residue1, residue2 ] )
        return output

    def get_self_contacting_residue_network( self, domains, dist_cutoff = 5.0 ):
        # CA distance <= 5.0
        aGraph = nx.Graph()
        #print(len(self.residues_list))
        for residue1 in self.residues_list:
            for residue2 in self.residues_list:
                if residue1 == residue2: continue
                dist = distance( residue1.get_CA(), residue2.get_CA() )
                if dist <= dist_cutoff:
                    aGraph.add_edge(residue1.pos, residue2.pos)
        for residue_pos in aGraph.nodes:
            aGraph.nodes[residue_pos]['label'] = "%d - %s" % (residue_pos, self.residues[residue_pos].resname)

        colors = ['red', 'yellow', 'green', 'cyan', 'magenta']
        c = 0
        for dom in domains:
            #print('dom\n', dom)
            for i in range(dom[0], dom[1]+1):
                try:
                    aGraph.nodes[i]['color'] = colors[c]
                except:
                    continue
            if c==5: c=0
            else: c += 1
        
        return aGraph

    def add_atom( self, resname, atomname, pos, x, y, z ):
        aAtom = Atom( resname, atomname, pos, x, y, z )
        if pos not in self.residues:
            aResidue = Residue( pos, resname )
            self.residues[ pos ] = aResidue
            self.residues_list.append( aResidue )
        self.residues[ pos ].atoms.append( aAtom )


class Ligand( Chain ):
    def get_center( self ):
        x = numpy.average( [ atom.x for atom in self.atoms ] )
        y = numpy.average( [ atom.y for atom in self.atoms ] )
        z = numpy.average( [ atom.z for atom in self.atoms ] )
        aAtom = Atom( "", "", x, y, z )
        return aAtom 

    def get_contacting_residues( self, chain, dist_cutoff = 5.0 ):
        pos_set = set()
        for atom in self.residues.values()[0].atoms:
            for residue in chain.residues.values():
                center = residue.get_center()
                dist = distance( center, atom )
                if dist <= dist_cutoff:
                    pos_set.add( atom.pos )
        return pos_set

class PDB:
    def __init__( self ):
        self.chains = {}
        self.ligands = {}
    pass
 
    def get_all_contact_chains( self ):
        chains = self.chains.keys()
        for i in range( len(chains ) ):
            chain1 = self.chains[ chains[i] ]
            for j in range( i+1, len(chains) ):
                chain2 = self.chains[ chains[i] ]
                if chain1.is_contacting( chain2 ):
                    print( chain1, chain2 )

    def get_all_self_dist( self, chain_ID ):
        chain = self.chains[ chain_ID ]
        return chain.all_to_all_distance()

    def get_all_ligand_contacts( self ):
        for name, pos in self.ligands:
            ligand = self.ligands[ (name, pos) ]
            for chain in self.chains.values():
                print( name, pos, "--", chain.chain_id, ligand.get_contacting_residues( chain ) )


def distance( atom1, atom2 ):
    return math.sqrt( (atom1.x-atom2.x)**2 + (atom1.y-atom2.y)**2 + (atom1.z-atom2.z)**2 )

def parsePDB( filepath ):
    aPDB = PDB()

    initialPos = dict()
    
    f = gzip.open(filepath, mode='rb')
    for line in f:
        line = line.decode('ascii')
        record = line[:6]
        if record == "ATOM  ":
            chain_id = line[21]
            atomname = line[12:16].strip()
            resname = line[17:20] 
            x = float( line[30:38] )
            y = float( line[38:46] )
            z = float( line[46:54] )
            pos = int( line[22:26] )
            try:
                pos = int( line[22:27] )
            except:
                return 0

            if chain_id not in initialPos.keys():
                initialPos[chain_id] = pos - 1 
            pos -= initialPos[chain_id]

            if chain_id not in aPDB.chains:
                aPDB.chains[ chain_id ] = Chain( chain_id )
            if atomname == 'CA':
               aPDB.chains[ chain_id ].add_atom( resname, atomname, pos, x, y, z )

    f.close()

    return aPDB

def domains(qrdict, name):
    '''
    returns a list with all the domains found in the protein
    '''
    l = []
    # for hit in queryresult
    for h in qrdict[name]:
        # for hsp in hit
        for hsp in h:
            # add 1 as it is a zero-based and half-open interval
            l.append(tuple([hsp.query_range[0]+1, hsp.query_range[1]]))
    return l

def find_unknown_res(aPDB, chainID, domdict, protname):
    # finding unknown residues
    unknown_residues = set( range( 1, len(aPDB.chains[chainID].residues_list) + 1 ) )
    
    for h in domdict[(protname.upper()+':'+chainID)]:
        for hsp in h:
            unknown_residues = unknown_residues - set(range(hsp.query_range[0]+1, hsp.query_range[1]+1))
    return unknown_residues


def save_dom_info(domdict, protname):
    filepath = "compressed_D/pdb"+protname+".ent.gz"
    aPDB = parsePDB( filepath )
    # if the protein redidues have insertion codes do not run the program
    if aPDB == 0:
        return aPDB
    domain_unknown_residues_interactions = shelve.open("domain_interactions.db")
    
    for chainID in aPDB.chains:
        if (protname.upper()+':'+chainID) not in domdict.keys():
            continue
        D = domains(domdict, (protname.upper()+':'+chainID))
        
        nx_graph = aPDB.chains[chainID].get_self_contacting_residue_network(D)
        #interactive_graph(nx_graph)
        
        unknown_residues = find_unknown_res(aPDB, chainID, domdict, protname)
        
        for h in domdict[(protname.upper()+':'+chainID)]:
            for hsp in h:
                unique_domain_id = "%s__%s__%d__%d" % ( protname+':'+chainID, h.accession, hsp.query_range[0]+1, hsp.query_range[1] )
                interacting_residues = set()
                for i in range(hsp.query_range[0]+1, hsp.query_range[1]+1):
                    if i in nx_graph:
                        for j in nx_graph[i]:
                            # if j is not inside the domain --> add to interacting residues    
                            if j not in range(hsp.query_range[0]+1, hsp.query_range[1]+1):
                                interacting_residues.add(j)
                if len(interacting_residues) != 0:
                    domain_unknown_residues_interactions[unique_domain_id] = list(interacting_residues)
                    domain_unknown_residues_interactions[unique_domain_id].sort()
        print(unique_domain_id)
        print(interacting_residues)

    domain_unknown_residues_interactions.close()

def static_graph(nx_graph):
    pos = nx.spring_layout(nx_graph, seed=7) #, iterations =300)  # positions for all nodes - seed for reproducibility
    
    # nodes
    nx.draw_networkx_nodes(nx_graph, pos, node_size=300)

    # edges
    nx.draw_networkx_edges(nx_graph, pos, width=5)
    
    # node labels
    nx.draw_networkx_labels(nx_graph, pos, font_size=8, font_family="sans-serif")
    # edge weight labels
    # edge_labels = nx.get_edge_attributes(G, "weight")
    #nx.draw_networkx_edge_labels(G, pos, edge_labels)

    ax = plt.gca()
    ax.margins(0.08)
    plt.axis("off")
    plt.tight_layout()
    plt.show()

def interactive_graph(nx_graph, x_size = '1000px', y_size = '1000px', html_file = 'nx.html'):
    #nt = Network(x_size, y_size)
    nt = Network(height="1000px", width="100%", bgcolor="#FFFFFF", font_color="#000000", notebook = True, heading = '') #, select_menu=True, filter_menu = True)

    nt.from_nx(nx_graph)
    nt.show_buttons(filter_=['nodes', 'edges', 'physics'])
    nt.show(html_file)

    display(HTML(html_file))

if __name__ == "__main__":
    from datetime import datetime
    init = datetime.now()
    domdict = shelve.open('domdict.shelve')

    c = 1
    protlist = [item[:-2].lower() for item in domdict.keys()]
    total = len(set(protlist))
    for f in os.listdir('data'):
        protname = f[3:-4]
        print(protname, '\t', c, '/', total, sep='')
        if protname in protlist:
            save_dom_info(domdict, protname)
        c += 1
    domdict.close()
    print(init, datetime.now(), sep='\t')
