import networkx

"""
Coding convention: http://google.github.io/styleguide/pyguide.html
"""

class HuRI:
    def __init__(self):
        self.proteins = {}
        self.g = networkx.Graph() # protein-protein interaction information

    def init(self):
        self.read_huri("./Y2H/data/HuRI/HuRI.pmi")
        self.read_pfam("./Y2H/data/HuRI/HuRI.fa.pfam")

    def summary(self):
        print("number of proteins:%s\n" % len(self.proteins))

        cnt = 0
        for name in self.proteins:
            aProtein = self.proteins[name]
            aProtein.summary()
            if cnt >= 10: break
            cnt += 1

    def read_huri(self, filepath):
        # Read HuRI network file and fill the information of protein and network.
        # uniprotkb:Q4VC05-2 uniprotkb:Q92624 ensembl:ENST00000538010.5|ensembl:ENSP00000445868.1|ensembl:ENSG00000110987.8 ensembl:ENST00000083182.7|ensembl:ENSP00000083182.3|ensembl:ENSG00000062725.9 human orfeome collection:53866(author assigned name) human orfeome collection:820(author assigned name) psi-mi:MI:0397(two hybrid array) Luck et al.(2019) unassigned1304 taxid:9606(Homo Sapiens) taxid:9606(Homo Sapiens) psi-mi:MI:0915(physical association) - - author score: 0.876896302095 - - - psi-mi:MI:0496(bait) psi-mi:MI:0496(prey) psi-mi:MI:0326(protein) psi-mi:MI:0326(protein) - - - comment:vector name: pDEST-DB|comment:centromeric vector|comment:yeast strain: Y8930 comment:vector name: pDEST-AD|comment:centromeric vector|comment:yeast strain: Y8800 comment:Found in screens 1. taxid:4932(Saccharomyces cerevisiae) - 2019/07/25 - - - - - gal4 dna binding domain:n-n (DB domain (n-terminal)) gal4 activation domain:n-n (AD domain (n-terminal)) - - psi-mi:MI1180(partial DNA sequence identification) psi-mi:MI1180(partial DNA sequence identification)
        # uniprotkb:A0A0A0MT20 uniprotkb:Q8IY31-3 ensembl:ENST00000433140.1|ensembl:ENSP00000411201.1|ensembl:ENSG00000138080.13 ensembl:ENST00000357896.7|ensembl:ENSP00000350570.3|ensembl:ENSG00000109083.13 human orfeome collection:8967(author assigned name) human orfeome collection:4088(author assigned name) psi-mi:MI:0397(two hybrid array) Luck et al.(2019) unassigned1304 taxid:9606(Homo Sapiens) taxid:9606(Homo Sapiens) psi-mi:MI:0915(physical association) - - author score: 0.904058209056 - - - psi-mi:MI:0496(bait) psi-mi:MI:0496(prey) psi-mi:MI:0326(protein) psi-mi:MI:0326(protein) - - - comment:vector name: pDEST-DB|comment:centromeric vector|comment:yeast strain: Y8930 comment:vector name: pAR68|comment:2u vector|comment:yeast strain: Y8800 comment:Found in screens 8. taxid:4932(Saccharomyces cerevisiae) - 2019/07/25 - - - - - gal4 dna binding domain:n-n (DB domain (n-terminal)) gal4 activation domain:c-c (AD domain (c-terminal)) - - psi-mi:MI1180(partial DNA sequence identification) psi-mi:MI1180(partial DNA sequence identification)
        f = open(filepath)
        for line in f.readlines():
            fields = line.split()
            enst1 = fields[2].split("|")[0][8:] # ENST00000538010.5
            enst2 = fields[3].split("|")[0][8:] # ENST00000083182.7
            self.g.add_edge(enst1, enst2)
        f.close()
        
        for enst in self.g.nodes():
            aProtein = Protein(enst)
            self.proteins[enst] = aProtein

    def read_pfam(self, filepath):
        # fills domain information for each protein
        # <seq id> <alignment start> <alignment end> <envelope start> <envelope end> <hmm acc> <hmm name> <type> <hmm start> <hmm end> <hmm length> <bit score> <E-value> <significance> <clan>
        # ENST00000428680.6      14     60     14     61 PF14570.7   zf-RING_4         Domain     1    47    48     71.5   3.5e-20   1 CL0229   
        f = open(filepath)
        for line in f.readlines():
            if line[0] == "#": continue
            fields = line.split()
            enst = fields[0]
            domain = fields[6]
            if enst in self.proteins:
                aProtein = self.proteins[enst]
                aProtein.add_domain(domain)
        pass

class Protein:
    def __init__(self, name):
        self.name = name
        self.domains = [] # domain list associated with this protein
        self.linearmotifs = [] # linear motif list associated with this protein

    def add_domain(self, domain):
        self.domains.append(domain)

    def summary(self):
        print("Protein Name:", self.name)
        print("Domains:", self.domains)
        print("LinearMotifs:", self.linearmotifs)

class DDIInfo:
    # Information about DDI pairs
    def __init__(self):
        self.DDIs = {}

    def read_3did(self):
        # fill self.DDIs 
        # For example "DOMAIN_A" have a partner "DOMAIN_B"
        # self.DDIs[ "DOMAIN_A" ]  = set()
        # self.DDIs[ "DOMAIN_A" ].add( "DOMAIN_B" )
        # Also fill reverse information, so that DDI information can be easily explored
        # self.DDIs[ "DOMAIN_B" ]  = set()
        # self.DDIs[ "DOMAIN_B" ].add( "DOMAIN_A" )
        pass
    
class DLIInfo:
    # Information about DLI pairs
    def __init__(self):
        self.DLIs = {}
        self.LDIs = {}

    def read_elm(self):
        # fill self.DDIs 
        # For example "DOMAIN_A" have a partner "LinearMotif_B"
        # self.DLIs[ "DOMAIN_A" ]  = set()
        # self.DLIs[ "DOMAIN_A" ].add( "LinearMotif_B" )
        # Also fill reverse information in different dictionary, so that DDI information can be easily explored
        # self.LDIs[ "LinearMotif_B" ]  = set()
        # self.LDIs[ "LinearMotif_B" ].add( "DOMAIN_A" )
        pass
    

if __name__ == "__main__":
    aHuRI = HuRI()
    aHuRI.init()
    
    aHuRI.summary()

    aDDIInfo = DDIInfo()
    aDLIInfo = DLIInfo()

    