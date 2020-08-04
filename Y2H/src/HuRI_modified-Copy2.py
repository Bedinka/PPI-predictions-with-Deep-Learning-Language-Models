#!/usr/bin/env python
# coding: utf-8

# In[11]:


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
            if cnt >= 1: break
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

        
        
class DLIInfo:

    def __init__(self):
        self.DLIs = {}
        self.LDIs = {}
        self.DMIs = {}
        self.MDIs = {} 
        
    def init(self):
        self.read_elm(r"C:\Users\Usuario\Desktop\CRAG\data\elm_interaction_domains.tsv")
        self.read_DMI_3DID(r"C:\Users\Usuario\Desktop\CRAG\data\3did_DMI_flat_Apr_10_2020.dat")
        
    def read_elm(self, filepath):
        f = open(filepath)
        for line in f.readlines():
            fields = line.split()
            LinearMotif_B = fields[0]
            Domain_A = fields[2]
            Hmm = fields[1]
            if Domain_A in self.DLIs and LinearMotif_B in self.LDIs:
                self.DLIs[Domain_A].add(LinearMotif_B)
                self.LDIs[LinearMotif_B].add(Domain_A)
            elif Domain_A in self.DLIs and LinearMotif_B not in self.LDIs:
                self.DLIs[Domain_A].add(LinearMotif_B)
                self.LDIs[LinearMotif_B] = set()
                self.LDIs[LinearMotif_B].add(Domain_A)
            elif Domain_A not in self.DLIs and LinearMotif_B in self.LDIs:
                self.DLIs[Domain_A] = set()
                self.DLIs[Domain_A].add(LinearMotif_B)
                self.LDIs[LinearMotif_B].add(Domain_A)
            elif Domain_A not in self.DLIs and LinearMotif_B not in self.LDIs:
                self.DLIs[Domain_A] = set()
                self.DLIs[Domain_A].add(LinearMotif_B)
                self.LDIs[LinearMotif_B] = set()
                self.LDIs[LinearMotif_B].add(Domain_A)
    
    def read_DMI_3DID(self, filepath):
        f = open(filepath)
        for line in f.readlines():
            if "#=ID" in line:
                fields = line.split()
                Motif_B = fields[3]
                Domain_A = fields[1]
                if Domain_A in self.DMIs and Motif_B in self.MDIs:
                    self.DMIs[Domain_A].add(Motif_B)
                    self.MDIs[Motif_B].add(Domain_A)
                elif Domain_A in self.DMIs and Motif_B not in self.MDIs:
                    self.DMIs[Domain_A].add(Motif_B)
                    self.MDIs[Motif_B] = set()
                    self.MDIs[Motif_B].add(Domain_A)
                elif Domain_A not in self.DMIs and Motif_B in self.MDIs:
                    self.DMIs[Domain_A] = set()
                    self.DMIs[Domain_A].add(Motif_B)
                    self.MDIs[Motif_B].add(Domain_A)
                elif Domain_A not in self.DMIs and Motif_B not in self.MDIs:
                    self.DMIs[Domain_A] = set()
                    self.DMIs[Domain_A].add(Motif_B)
                    self.MDIs[Motif_B] = set()
                    self.MDIs[Motif_B].add(Domain_A)
                    
    def summary(self):
        #print("Complete list of DLIs:", self.DLIs)
        #print("Complete list of LDIs", self.LDIs)
        #print("Complete list of DMIs:", self.DMIs)
        #print("Complete list of MDIs:", self.MDIs)
        pass
    
class DDIInfo:
    # Information about DDI pairs
    def __init__(self):
        self.DDIs = {}
        self.particularDDIs = {}
        
    def init(self):
        self.read_3did_1(r"C:\Users\Usuario\Desktop\CRAG\data\3did_flat_Apr_10_2020.1.dat")
        self.read_3did_2(r"C:\Users\Usuario\Desktop\CRAG\data\3did_flat_Apr_10_2020.2.dat")
        
    def read_3did_1(self, filepath):
        f = open(filepath)
        for line in f.readlines():
            if "#=ID" in line:
                fields = line.split()
                Domain_B = fields[2]
                Domain_A = fields[1]
                if Domain_A in self.DDIs and Domain_B in self.DDIs:
                    self.DDIs[Domain_A].add(Domain_B)
                    self.DDIs[Domain_B].add(Domain_A)
                elif Domain_A in self.DDIs and Domain_B not in self.DDIs:
                    self.DDIs[Domain_A].add(Domain_B)
                    self.DDIs[Domain_B] = set()
                    self.DDIs[Domain_B].add(Domain_A)
                elif Domain_A not in self.DDIs and Domain_B in self.DDIs:
                    self.DDIs[Domain_A] = set()
                    self.DDIs[Domain_A].add(Domain_B)
                    self.DDIs[Domain_B].add(Domain_A)
                elif Domain_A not in self.DDIs and Domain_B not in self.DDIs:
                    self.DDIs[Domain_A] = set()
                    self.DDIs[Domain_A].add(Domain_B)
                    self.DDIs[Domain_B] = set()
                    self.DDIs[Domain_B].add(Domain_A)

    
    def read_3did_2(self, filepath):
        f = open(filepath)
        for line in f.readlines():
            if "#=ID" in line:
                fields = line.split()
                Domain_B = fields[2]
                Domain_A = fields[1]
                if Domain_A in self.DDIs and Domain_B in self.DDIs:
                    self.DDIs[Domain_A].add(Domain_B)
                    self.DDIs[Domain_B].add(Domain_A)
                elif Domain_A in self.DDIs and Domain_B not in self.DDIs:
                    self.DDIs[Domain_A].add(Domain_B)
                    self.DDIs[Domain_B] = set()
                    self.DDIs[Domain_B].add(Domain_A)
                elif Domain_A not in self.DDIs and Domain_B in self.DDIs:
                    self.DDIs[Domain_A] = set()
                    self.DDIs[Domain_A].add(Domain_B)
                    self.DDIs[Domain_B].add(Domain_A)
                elif Domain_A not in self.DDIs and Domain_B not in self.DDIs:
                    self.DDIs[Domain_A] = set()
                    self.DDIs[Domain_A].add(Domain_B)
                    self.DDIs[Domain_B] = set()
                    self.DDIs[Domain_B].add(Domain_A)

    def add_ddi(self,domain):
        for k,v in self.domain.items():
            if domain in self.domain.keys():
                self.particularDDIs[domain] = set()
                self.particularDDIs[domain] = self.domain.get(domain)
            
    def summary(self):
        #print("Complete list of DDIs:", self.DDIs) 
        print(self.particularDDIs)
        
        
class Protein:
    def __init__(self, name):
        self.name = name
        self.domains = [] # domain list associated with this protein
        self.linearmotifs = [] # linear motif list associated with this protein
        self.hmms = [] #hmm acc associated with this protein
        self.interactions = []
        
    def add_domain(self, domain):
        self.domains.append(domain)      

    def add_hmm(self, hmm):         
        self.hmms.append(hmm)   
        
    def add_interaction(self, enst2):
        if enst2 not in self.interactions:
            self.interactions.append(enst2)
            
    def domain_domain(self):
        for domain in self.domains:
            aDDIInfo = self.domains[domain]
            aDDIInfo.add_ddi(domain)
            
            
    #def add_dmi:
        
    #def add_dli:
        
    def summary(self):
        print("Protein Name:", self.name)
        print("Domains:", self.domains)
        #print("LinearMotifs:", self.linearmotifs)
        print("Interacts with:", self.interactions)   

if __name__ == "__main__":
    aHuRI = HuRI()
    aHuRI.init()
    
    aHuRI.summary()

    aDDIInfo = DDIInfo()
    aDLIInfo = DLIInfo()
    aDLIInfo.init()
    aDDIInfo.init()
    
    aDLIInfo.summary()
    aDDIInfo.summary()
    


# In[ ]:





# In[ ]:




