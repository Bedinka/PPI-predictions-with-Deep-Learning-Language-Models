#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd


# In[2]:


df1 = pd.read_csv ("huri.csv", sep=';')
print(df1)


# In[13]:


df2 = pd.read_csv ("DDI_for_map.txt", sep='\t')
print(df2)


# In[14]:


df_merge = pd.merge(df1, df2, how='inner')


# In[15]:


df_merge


# In[16]:


df_merge.to_csv(r'C:\Users\Marina\Desktop\common_huri_ddi.csv')


# In[ ]:




