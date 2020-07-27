import networkx

class HuRI:
    def __init__(self):
        self.proteins = []
        self.g = networkx.Graph() # protein-protein interaction information

    def ReadHuRI(self):
        # Read HuRI network file and fill the information of protein and network.
        pass
    
    def ReadPFam(self):
        # fills domain information for each protein
        pass

class Protein:
    def __init__(self, name):
        self.name = name
        self.domains = [] # domain list associated with this protein
        self.linearmotifs = [] # linear motif list associated with this protein

class DDIInfo:
    # Information about DDI pairs
    def __init__(self):
        self.DDIs = {}

    def Read3DID(self):
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

    def ReadELM(self):
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
    aHuRI.ReadPFam()
    aHuRI.ReadHuRI()

    aDDIInfo = DDIInfo()
    aDLIInfo = DLIInfo()

    