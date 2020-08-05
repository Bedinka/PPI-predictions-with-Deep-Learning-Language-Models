#!/usr/bin/env python
# coding: utf-8

# In[28]:


import networkx

"""
Coding convention: http://google.github.io/styleguide/pyguide.html
"""

class HuRI:
    def __init__(self):
        self.proteins = {}
        self.g = networkx.Graph() # protein-protein interaction information

    def init(self):
        self.read_huri(r"C:\Users\Usuario\Desktop\CRAG\data\HuRI.pmi")
        self.read_pfam(r"C:\Users\Usuario\Desktop\CRAG\data\HuRI.fa.pfam")

    def summary(self):
        print("number of proteins:%s\n" % len(self.proteins))
        cnt = 0
        for name in self.proteins:
            aProtein = self.proteins[name]
            aProtein.summary()
            if cnt >= 0: break
            cnt += 1

    def read_huri(self, filepath):

        f = open(filepath)
        for line in f.readlines():
            fields = line.split()
            enst1 = fields[2].split("|")[0][8:] # ENST00000538010.5
            enst2 = fields[3].split("|")[0][8:] # ENST00000083182.
            self.g.add_edge(enst1, enst2)
        f.close()
        
        for enst in self.g.nodes():
            aProtein = Protein(enst)
            self.proteins[enst] = aProtein
            aProtein = Protein(enst1)
            self.proteins[enst1] = aProtein

        f = open(filepath)
        for line in f.readlines():
            fields = line.split()
            enst = fields[2].split("|")[0][8:] # ENST00000538010.5
            enst2 = fields[3].split("|")[0][8:] # ENST00000083182.
            if enst in self.proteins:
                aProtein = self.proteins[enst]
                aProtein.add_interaction(enst2)
        
    def read_pfam(self, filepath):
        f = open(filepath)
        for line in f.readlines():
            if line[0] == "#": continue
            fields = line.split()
            enst = fields[0]
            hmm = fields[5][:7]
            domain = fields[6]
            if enst in self.proteins:
                aProtein = self.proteins[enst]
                aProtein.add_domain(domain)
                aProtein.add_hmm(hmm)
    
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

    def get_PPIs_known_DLI(self, aDLIInfo):
        _g = networkx.Graph()

        for enstA in self.g.nodes():
            proteinA = self.proteins[enstA]

            #print(proteinA.domains)
            #print(aDLIInfo.get_known_DLI_Domains(proteinA))
            #print("\n")

            for enstB in self.g[enstA]:
                proteinB = self.proteins[enstB]
                if aDLIInfo.have_known_DLIs(proteinA, proteinB) == True:
                    _g.add_edge( enstA, enstB )
        return _g  
    
class DDIInfo:
    # Information about DDI pairs
    def __init__(self):
        self.DDIs = networkx.Graph()
        self.completeDDIs = {}
        self.read_3did_1(r"C:\Users\Usuario\Desktop\CRAG\data\3did_flat_Apr_10_2020.1.dat")
        self.read_3did_2(r"C:\Users\Usuario\Desktop\CRAG\data\3did_flat_Apr_10_2020.2.dat")
   
    def read_3did_1(self, filepath):
        f = open(filepath)
        for line in f.readlines():
            if "#=ID" in line:
                fields = line.split()
                Domain_B = fields[2]
                Domain_A = fields[1]
                self.DDIs.add_edge(Domain_A, Domain_B)
                if Domain_A in self.completeDDIs and Domain_B in self.completeDDIs:
                    self.completeDDIs[Domain_A].add(Domain_B)
                    self.completeDDIs[Domain_B].add(Domain_A)
                elif Domain_A in self.completeDDIs and Domain_B not in self.completeDDIs:
                    self.completeDDIs[Domain_A].add(Domain_B)
                    self.completeDDIs[Domain_B] = set()
                    self.completeDDIs[Domain_B].add(Domain_A)
                elif Domain_A not in self.completeDDIs and Domain_B in self.completeDDIs:
                    self.completeDDIs[Domain_A] = set()
                    self.completeDDIs[Domain_A].add(Domain_B)
                    self.completeDDIs[Domain_B].add(Domain_A)
                elif Domain_A not in self.completeDDIs and Domain_B not in self.completeDDIs:
                    self.completeDDIs[Domain_A] = set()
                    self.completeDDIs[Domain_A].add(Domain_B)
                    self.completeDDIs[Domain_B] = set()
                    self.completeDDIs[Domain_B].add(Domain_A)

    
    def read_3did_2(self, filepath):
        f = open(filepath)
        for line in f.readlines():
            if "#=ID" in line:
                fields = line.split()
                Domain_B = fields[2]
                Domain_A = fields[1]
                self.DDIs.add_edge(Domain_A, Domain_B)
                if Domain_A in self.completeDDIs and Domain_B in self.completeDDIs:
                    self.completeDDIs[Domain_A].add(Domain_B)
                    self.completeDDIs[Domain_B].add(Domain_A)
                elif Domain_A in self.completeDDIs and Domain_B not in self.completeDDIs:
                    self.completeDDIs[Domain_A].add(Domain_B)
                    self.completeDDIs[Domain_B] = set()
                    self.completeDDIs[Domain_B].add(Domain_A)
                elif Domain_A not in self.completeDDIs and Domain_B in self.completeDDIs:
                    self.completeDDIs[Domain_A] = set()
                    self.completeDDIs[Domain_A].add(Domain_B)
                    self.completeDDIs[Domain_B].add(Domain_A)
                elif Domain_A not in self.completeDDIs and Domain_B not in self.completeDDIs:
                    self.completeDDIs[Domain_A] = set()
                    self.completeDDIs[Domain_A].add(Domain_B)
                    self.completeDDIs[Domain_B] = set()
                    self.completeDDIs[Domain_B].add(Domain_A)
                    
    def have_known_DDIs(self, proteinA, proteinB):
        for Domain_A in proteinA.domains:
            for Domain_B in proteinB.domains:
                if self.DDIs.has_edge(Domain_A, Domain_B):
                    return True
        return False
    
    
    def count_known_DDIs(self, proteinA, proteinB):
        known_ddi_cnt = 0
        for Domain_A in proteinA.domains:
            for Domain_B in proteinB.domains:
                if self.DDIs.has_edge(Domain_A, Domain_B):
                    known_ddi_cnt += 1
        return known_ddi_cnt
    
    def get_known_DDI_Domains(self, proteinA):
        return [ domain for domain in proteinA.domains if domain in self.DDIs ]  
    
    def summary(self):
        #print("Complete list of DDIs:", self.completeDDIs) 
        print("------------------------------------------------------------------------------------------------ ")
        print("DOMAIN DOMAIN INTERACTIONS:")
        print("number of domains:%s\n" % len(self.DDIs))
        cnt = 0
        for Domain_A in self.DDIs:
            print(Domain_A)
            print(self.DDIs[Domain_A])
            print("\n")
            if cnt >= 10: break
            cnt += 1
            
class DLIInfo:

    def __init__(self):
        self.completeDLIs = {}
        self.completeLDIs = {}
        self.completeDMIs = {}
        self.completeMDIs = {} 
        self.DLIs = networkx.Graph()
        self.DMIs = networkx.Graph()
        self.read_elm(r"C:\Users\Usuario\Desktop\CRAG\data\elm_interaction_domains.tsv")
        self.read_DMI_3DID(r"C:\Users\Usuario\Desktop\CRAG\data\3did_DMI_flat_Apr_10_2020.dat")
    
    def read_elm(self, filepath):
        f = open(filepath)
        for line in f.readlines():
            fields = line.split()
            LinearMotif_B = fields[0]
            Domain_A = fields[2]
            Hmm = fields[1]
            self.DLIs.add_edge(Domain_A, LinearMotif_B)
            if Domain_A in self.completeDLIs and LinearMotif_B in self.completeLDIs:
                self.completeDLIs[Domain_A].add(LinearMotif_B)
                self.completeLDIs[LinearMotif_B].add(Domain_A)
            elif Domain_A in self.completeDLIs and LinearMotif_B not in self.completeLDIs:
                self.completeDLIs[Domain_A].add(LinearMotif_B)
                self.completeLDIs[LinearMotif_B] = set()
                self.completeLDIs[LinearMotif_B].add(Domain_A)
            elif Domain_A not in self.completeDLIs and LinearMotif_B in self.completeLDIs:
                self.completeDLIs[Domain_A] = set()
                self.completeDLIs[Domain_A].add(LinearMotif_B)
                self.completeLDIs[LinearMotif_B].add(Domain_A)
            elif Domain_A not in self.completeDLIs and LinearMotif_B not in self.completeLDIs:
                self.completeDLIs[Domain_A] = set()
                self.completeDLIs[Domain_A].add(LinearMotif_B)
                self.completeLDIs[LinearMotif_B] = set()
                self.completeLDIs[LinearMotif_B].add(Domain_A)
    
    def read_DMI_3DID(self, filepath):
        f = open(filepath)
        for line in f.readlines():
            if "#=ID" in line:
                fields = line.split()
                Motif_B = fields[3]
                Domain_A = fields[1]
                self.DMIs.add_edge(Domain_A, Motif_B)
                if Domain_A in self.completeDMIs and Motif_B in self.completeMDIs:
                    self.completeDMIs[Domain_A].add(Motif_B)
                    self.completeMDIs[Motif_B].add(Domain_A)
                elif Domain_A in self.completeDMIs and Motif_B not in self.completeMDIs:
                    self.completeDMIs[Domain_A].add(Motif_B)
                    self.completeMDIs[Motif_B] = set()
                    self.completeMDIs[Motif_B].add(Domain_A)
                elif Domain_A not in self.completeDMIs and Motif_B in self.completeMDIs:
                    self.completeDMIs[Domain_A] = set()
                    self.completeDMIs[Domain_A].add(Motif_B)
                    self.completeMDIs[Motif_B].add(Domain_A)
                elif Domain_A not in self.completeDMIs and Motif_B not in self.completeMDIs:
                    self.completeDMIs[Domain_A] = set()
                    self.completeDMIs[Domain_A].add(Motif_B)
                    self.completeMDIs[Motif_B] = set()
                    self.completeMDIs[Motif_B].add(Domain_A)
                    
    def have_known_DLIs(self, proteinA, proteinB):
        for Domain_A in proteinA.domains:
            for (u,v) in self.DLIs.edges():
                if Domain_A in self.DLIs.nodes:
                    return True
        return False
        
    def count_known_DLIs(self, proteinA, proteinB):
        known_dli_cnt = 0
        for Domain_A in proteinA.domains:
            for (u,v) in self.DLIs.edges():
                if Domain_A in self.DLIs.nodes:
                    known_dli_cnt += 1
        return known_dli_cnt
    
    def get_known_DLI_Domains(self, proteinA):
        return [ domain for domain in proteinA.domains if domain in self.DLIs ]    

    def have_known_DMIs(self, proteinA, proteinB):
        for Domain_A in proteinA.domains:
            for (u,v) in self.DMIs.edges():
                if Domain_A in self.DMIs.nodes:
                    return True
        return False
        
    def count_known_DMIs(self, proteinA, proteinB):
        known_dmi_cnt = 0
        for Domain_A in proteinA.domains:
            for (u,v) in self.DMIs.edges():
                if Domain_A in self.DMIs.nodes:
                    known_dmi_cnt += 1
        return known_dmi_cnt
    
    def get_known_DMI_Domains(self, proteinA):
        return [ domain for domain in proteinA.domains if domain in self.DMIs ] 
    
    def summary(self):
        #print("Complete list of DLIs:", self.DLIs)
        #print("Complete list of LDIs", self.LDIs)
        #print("Complete list of DMIs:", self.DMIs)
        #print("Complete list of MDIs:", self.MDIs)

        print("------------------------------------------------------------------------------------------------ ")
        print("DOMAIN LINEAR-MOTIF INTERACTIONS:")
        print("number of domains:%s\n" % len(self.DLIs))
        cnt2 = 0
        for Domain_A in self.DLIs:
            print(Domain_A)
            print(self.DLIs[Domain_A])
            print("\n")
            if cnt2 >= 10: break
            cnt2 += 1

        print("------------------------------------------------------------------------------------------------ ")
        print("DOMAIN MOTIF INTERACTIONS:")
        print("number of domains:%s\n" % len(self.DMIs))
        cnt3 = 0
        for Domain_A in self.DMIs:
            print(Domain_A)
            print(self.DMIs[Domain_A])
            print("\n")
            if cnt3 >= 10: break
            cnt3 += 1
        
class Protein:
    def __init__(self, name):
        self.name = name
        self.domains = [] # domain list associated with this protein
        self.linearmotifs = [] # linear motif list associated with this protein
        self.hmms = [] #hmm acc associated with this protein
        self.interactions = []
        
    def add_domain(self, domain):
        self.domains.append(domain) 
        
    def add_linearmofif(self, linearmotif):
        self.linearmotifs.append(linearmotif)

    def add_hmm(self, hmm):         
        self.hmms.append(hmm)   
        
    def add_interaction(self, enst2):
        if enst2 not in self.interactions:
            self.interactions.append(enst2)
            
    def summary(self):
        print("Protein Name:", self.name)
        print("Domains:", self.domains)
        print("LinearMotifs:", self.linearmotifs)
        print("Interacts with:", self.interactions) 
        
#def ANALYSIS_1_DDI_MEDIATED_INTERACTIONS( aHuRI, aDDIInfo ):
    #print("\n#====================================")
    #print("ANALYSIS_1_DDI_MEDIATED_INTERACTION")
    #print("#====================================")
    #_g = aHuRI.get_PPIs_known_DDI(aDDIInfo)
   # print("Total HuRI proteins in the network", len(aHuRI.g.nodes()))
  #  print("Total HuRI interactions", len(aHuRI.g.edges()))
 #   print("Total HuRI proteins with known DDIs",len(_g.nodes()))
#    print("Total HuRI interactions with known DDIs",len(_g.edges()))
    
#def ANALYSIS_2_DLI_MEDIATED_INTERACTIONS( aHuRI, aDLIInfo ):
 #   print("\n#====================================")
  #  print("ANALYSIS_2_DLI_MEDIATED_INTERACTION")
   # print("#====================================")
    #_g = aHuRI.get_PPIs_known_DLI(aDLIInfo)
    #print("Total HuRI proteins in the network", len(aHuRI.g.nodes()))
    #print("Total HuRI interactions", len(aHuRI.g.edges()))
    #print("Total HuRI proteins with known DLIs",len(_g.nodes()))
    #print("Total HuRI interactions with known DLIs",len(_g.edges()))
    
if __name__ == "__main__":
    aHuRI = HuRI()
    aHuRI.init()
    
    aHuRI.summary()

    aDDIInfo = DDIInfo()
    aDDIInfo.summary()
    
    aDLIInfo = DLIInfo()
    aDLIInfo.summary()
    
    #ANALYSIS_1_DDI_MEDIATED_INTERACTIONS( aHuRI, aDDIInfo )
    #ANALYSIS_2_DLI_MEDIATED_INTERACTIONS( aHuRI, aDLIInfo )
    


# In[ ]:





# In[ ]:




