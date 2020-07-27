#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np


# In[2]:


data = pd.read_table("Mapping_output.txt", sep=";")


# In[5]:


df = pd.DataFrame(data)


# In[12]:


df= df.drop_duplicates()


# In[14]:


df["INTERACTING_PROTEINS"] = df["PROTEIN 1 (uniprotkb)"] + "   " + df["PROTEIN 2 (uniprotkb)"]


# In[15]:


df


# In[18]:


count_1 = df.groupby(["PROTEIN 1 (uniprotkb)","PROTEIN 2 (uniprotkb)"]).size()


# In[19]:


count_1


# In[20]:


count_1.to_csv("number_interactions_for_pair.csv")


# In[21]:


count_2 = pd.crosstab(df.INTERACTING_PROTEINS, df.INTERACTION)


# In[22]:


print(count_2)


# In[23]:


count_2.to_csv("detailed_interactions_for_pair.csv")


# In[ ]:




