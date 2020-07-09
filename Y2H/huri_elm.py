#!/usr/bin/env python
# coding: utf-8

# In[49]:


import pandas as pd


# In[50]:


df1 = pd.read_csv ("huri.csv", sep=';')
print(df)


# In[51]:


df2 =  pd.read_csv ("elm_interaction_domains.csv", sep=';')
print(df2)         

#note that column with INTERACTION DOMAIN DESCRIPTION IS DUPLICATED WITH NAME <hmm name> in order to compare with column <hmm name> of the HuRI file 


# In[52]:


import numpy as np


# In[53]:


df_merge = pd.merge(df, df2, how='inner')   #to obtain common values between both files


# In[54]:


df_merge


# In[55]:


df_merge.to_csv(r'C:\Users\Marina\Desktop\common_huri_elm.csv')


# In[13]:


#OUTPUT FILE CONTAINS THE HURI INTERACTIONS THAT ARE EXPLAINED BY ELM WITH COMPLETE INFORMATION

