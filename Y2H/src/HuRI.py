import networkx

"""
Coding convention: http://google.github.io/styleguide/pyguide.html
"""

class HuRI:
    def __init__(self):
        self.proteins = {}
        self.g = networkx.Graph() # protein-protein interaction information

    def init(self):
        self.read_huri("./data/HuRI/HuRI.pmi")
        self.read_pfam("./data/HuRI/HuRI.fa.pfam")

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
        import os
        print(os.path.abspath(filepath))
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

    def get_PPIs_known_DDI(self, aDDIInfo):
        _g = networkx.Graph()

        for enstA in self.g.nodes():
            proteinA = self.proteins[enstA]

            #print(proteinA.domains)
            #print(aDDIInfo.get_known_DDI_Domains(proteinA))
            #print("\n")

            for enstB in self.g[enstA]:
                proteinB = self.proteins[enstB]
                if aDDIInfo.have_known_DDIs(proteinA, proteinB) == True:
                    _g.add_edge( enstA, enstB )
        return _g

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
        self.DDIs = networkx.Graph()
        self.read_3did("./data/3DID/3did_flat_Apr_10_2020.1.dat")

    def read_3did(self, filepath):
        f = open(filepath)
        for line in f.readlines():
            if line[:4] != "#=ID": continue
            fields = line.split("\t")
            domainA = fields[1]
            domainB = fields[2]
            self.DDIs.add_edge( domainA, domainB )
        f.close()
    
    def have_known_DDIs(self, proteinA, proteinB):
        for domainA in proteinA.domains:
            for domainB in proteinB.domains:
                if self.DDIs.has_edge(domainA, domainB):
                    return True
        return False

    def count_known_DDIs(self, proteinA, proteinB):
        known_ddi_cnt = 0
        for domainA in proteinA.domains:
            for domainB in proteinB.domains:
                if self.DDIs.has_edge(domainA, domainB):
                    known_ddi_cnt += 1
        return known_ddi_cnt

    def get_known_DDI_Domains(self, proteinA):
        return [ domain for domain in proteinA.domains if domain in self.DDIs ]

    def summary(self):
        print("number of domains:%s\n" % len(self.DDIs))

        cnt = 0
        for domainA in self.DDIs:
            print(domainA)
            print(self.DDIs[domainA])
            print("\n")
            if cnt >= 10: break
            cnt += 1

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
    

def ANALYSIS_1_DDI_MEDIATED_INTERACTIONS( aHuRI, aDDIInfo ):
    print("\n#====================================")
    print("ANALYSIS_1_DDI_MEDIATED_INTERACTION")
    print("#====================================")
    _g = aHuRI.get_PPIs_known_DDI(aDDIInfo)
    print("Total HuRI proteins in the network", len(aHuRI.g.nodes()))
    print("Total HuRI interactions", len(aHuRI.g.edges()))
    print("Total HuRI proteins with known DDIs",len(_g.nodes()))
    print("Total HuRI interactions with known DDIs",len(_g.edges()))


if __name__ == "__main__":
    aHuRI = HuRI()
    aHuRI.init()
    
    aHuRI.summary()

    aDDIInfo = DDIInfo()
    aDDIInfo.summary()

    aDLIInfo = DLIInfo()


    ANALYSIS_1_DDI_MEDIATED_INTERACTIONS( aHuRI, aDDIInfo )



    