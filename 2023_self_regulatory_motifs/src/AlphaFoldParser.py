import sys
import numpy
import math
import networkx as nx
from pyvis.network import Network

from IPython.core.display import display, HTML
import matplotlib.pyplot as plt

from Bio import SearchIO


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

    def get_center( self ):
        if self.center == None:
            x = numpy.average( [ atom.x for atom in self.atoms ] )
            y = numpy.average( [ atom.y for atom in self.atoms ] )
            z = numpy.average( [ atom.z for atom in self.atoms ] )
            self.center = Atom( "CENTER", "X", -1, x, y, z )
        return self.center

    def get_CA( self ):
        if self.CA == None:
            #print ("###",self.pos, self.atoms)
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

    def get_self_contacting_residue_network( self, dist_cutoff = 5.0 ):
        # CA distance <= 5.0
        aGraph = nx.Graph()
        #print(self.residues_list)
        for residue1 in self.residues_list:
            for residue2 in self.residues_list:
                if residue1 == residue2: continue
                dist = distance( residue1.get_CA(), residue2.get_CA() )
                if dist <= dist_cutoff:
                    aGraph.add_edge(residue1.pos, residue2.pos)
        for residue_pos in aGraph.nodes:
            aGraph.nodes[residue_pos]['label'] = "%d - %s" % (residue_pos, self.residues[residue_pos].resname)

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


def parsePDB( filename ):
    # ATOM   2692  N   LYS G  12      55.763  32.551  37.146  1.00 55.62           N
    # ATOM   2693  CA  LYS G  12      56.776  32.217  36.152  1.00 55.72           C
    # ATOM   2694  C   LYS G  12      56.408  31.077  35.204  1.00 55.02           C
    # 31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms.                       
    # 39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms.                            
    # 47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms.                            
    # X1 = ( 1.000000)*Xorig + ( 0.000000)*Yorig + ( 0.000000)*Zorig + (    0.000000)
    # Y1 = (-0.000000)*Xorig + ( 1.000000)*Yorig + ( 0.000000)*Zorig + (   -0.000000)
    # Z1 = (-0.000000)*Xorig + ( 0.000000)*Yorig + ( 1.000000)*Zorig + (   -0.000000)

    aPDB = PDB()

    chain = ""
    f = gzip.open(filepath, mode='rb')
    for line in f:
        line = line.decode('ascii')
        record = line[:6]
        print( line )
        #print (record + "#")
        if record == "ATOM  ":
            chain_id = line[21]
            atomname = line[12:16].strip()
            resname = line[17:20] 
            x = float( line[30:38] )
            y = float( line[38:46] )
            z = float( line[46:54] )
            pos = int( line[22:26] )
            if resname in LIGANDS:
                if ( resname, pos ) not in aPDB.ligands:
                    aPDB.ligands[ ( resname, pos ) ] = Ligand( chain_id )
                aPDB.ligands[ (resname, pos) ].add_atom( resname, atomname, pos, x, y, z )
                continue
            if chain_id not in aPDB.chains:
                aPDB.chains[ chain_id ] = Chain( chain_id )
            aPDB.chains[ chain_id ].add_atom( resname, atomname, pos, x, y, z )
        if record == "HETATM":
            chain_id = line[21]
            resname = line[17:20] 
            atomname = line[12:16].strip()
            x = float( line[30:38] )
            y = float( line[38:46] )
            z = float( line[46:54] )
            pos = int( line[22:26] )
            if resname in LIGANDS:
                if ( resname, pos ) not in aPDB.ligands:
                    aPDB.ligands[ ( resname, pos ) ] = Ligand( chain_id )
                aPDB.ligands[ (resname, pos) ].add_atom( resname, atomname, pos, x, y, z )
                continue
            if chain_id not in aPDB.chains:
                aPDB.chains[ chain_id ] = Chain( chain_id )
            aPDB.chains[ chain_id ].add_atom( resname, atomname, pos, x, y, z )
    f.close()

    return aPDB

def overlap(a,b):
    if a is None or b is None:
        return False
    return a[0] <= b[0] <= a[1] or b[0] <= a[0] <= b[1]

def domtbl_parser(file, E = 10e-6):
    qrdict = dict()
    evalue_filter = lambda hsp: hsp.evalue < E

    with open(file, 'r') as f:
        # transforms to zero-based and half-open intervals
        domtbl = SearchIO.parse(f, 'hmmscan3-domtab')
        for qresult in domtbl:
            for hit in qresult.hits:
                # filter the hsp on the hit by evalue
                filtered_hit = hit.filter(evalue_filter)
                # if no hsps of the hit pass the threshold pop the hit from the QueryResult
                if filtered_hit is None:
                    qresult.pop(hit.id)
                else:
                    prevhsp = None
                    index = 0
                    for hsp in filtered_hit:
                        # if two hsps overlap pop the worst one (they are sorted by evalue)
                        if overlap(prevhsp, hsp.query_range):
                            filtered_hit.pop(index)
                        prevhsp = hsp.query_range
                        index += 1
                    #save the "cleaned" hit to the QueryResult
                    qresult[hit.id] = filtered_hit
            # Only save QueryResults that contain information
            if len(qresult) > 0:
                qrdict[qresult.id] = qresult
        return qrdict
    
def domains(qrdict, name):
    '''
    returns a list with all the domains found in the protein
    '''
    l = []
    for h in qrdict[name]:
        for hsp in h:
            l.append(hsp.query_range)
    return l


def static_graph(nx_graph, domains):
    pos = nx.spring_layout(nx_graph, seed=7) #, iterations =300)  # positions for all nodes - seed for reproducibility
    
    color_map = ['blue']*(nx_graph.number_of_nodes())
    colors = ['red', 'green', 'purple']
    c = 0
    for dom in domains:
        for i in range(dom[0], dom[1]):
            color_map[i] = colors[c]
        c += 1
    
    # nodes
    nx.draw_networkx_nodes(nx_graph, pos, node_color=color_map, node_size=300)

    # edges
    nx.draw_networkx_edges(nx_graph, pos, width=5)
    
    # node labels
    nx.draw_networkx_labels(nx_graph, pos, font_size=8, font_family="sans-serif")
    # edge weight labels
    #edge_labels = nx.get_edge_attributes(G, "weight")
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
    import gzip
    wd = "/home/roger/YangLabIntern/2023_self_regulatory_motifs/"
    protname = "AF-A0A024RBG1-F1-model_v2"
    filepath = wd+"compressed_D/"+protname+".pdb.gz"

    domdict = domtbl_parser(wd+'domtbl.out')
    

    aPDB = parsePDB( filepath )

    #aPDB.get_all_ligand_contacts()
    print ( aPDB.chains )
    nx_graph = aPDB.chains['A'].get_self_contacting_residue_network()
    print ( nx_graph )
    
    D = domains(domdict, protname)

    print(domains(domdict, protname))
    static_graph(nx_graph, D)
    interactive_graph(nx_graph)
    print(domains(domdict, protname))

    '''
    import plotly.graph_objects as go

    # Create Edges
    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = G.nodes[edge[0]]['pos']
        x1, y1 = G.nodes[edge[1]]['pos']
        edge_x.append(x0)
        edge_x.append(x1)
        edge_x.append(None)
        edge_y.append(y0)
        edge_y.append(y1)
        edge_y.append(None)

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines')

    node_x = []
    node_y = []
    for node in G.nodes():
        x, y = G.nodes[node]['pos']
        node_x.append(x)
        node_y.append(y)

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers',
        hoverinfo='text',
        marker=dict(
            showscale=True,
            # colorscale options
            #'Greys' | 'YlGnBu' | 'Greens' | 'YlOrRd' | 'Bluered' | 'RdBu' |
            #'Reds' | 'Blues' | 'Picnic' | 'Rainbow' | 'Portland' | 'Jet' |
            #'Hot' | 'Blackbody' | 'Earth' | 'Electric' | 'Viridis' |
            colorscale='YlGnBu',
            reversescale=True,
            color=[],
            size=10,
            colorbar=dict(
                thickness=15,
                title='Node Connections',
                xanchor='left',
                titleside='right'
            ),
            line_width=2))



    # Color Node Points
    node_adjacencies = []
    node_text = []
    for node, adjacencies in enumerate(G.adjacency()):
        node_adjacencies.append(len(adjacencies[1]))
        node_text.append('# of connections: '+str(len(adjacencies[1])))

    node_trace.marker.color = node_adjacencies
    node_trace.text = node_text

    # Create Network Graph
    fig = go.Figure(data=[edge_trace, node_trace],
             layout=go.Layout(
                title='<br>Network graph made with Python',
                titlefont_size=16,
                showlegend=False,
                hovermode='closest',
                margin=dict(b=20,l=5,r=5,t=40),
                annotations=[ dict(
                    text="Python code: <a href='https://plotly.com/ipython-notebooks/network-graphs/'> https://plotly.com/ipython-notebooks/network-graphs/</a>",
                    showarrow=False,
                    xref="paper", yref="paper",
                    x=0.005, y=-0.002 ) ],
                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                )
    fig.show()
    '''

