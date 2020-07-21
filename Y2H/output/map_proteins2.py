#!/usr/bin/env python
# coding: utf-8

# In[7]:


import pandas as pd


# In[8]:


df1 = pd.read_csv ("proteins.csv", sep=';')
print(df1)


# In[9]:


df2 = pd.read_csv ("all_interactions2.csv", sep=';')
print(df2)


# In[10]:


df_merge = pd.merge(df1, df2, how='inner')


# In[11]:


df_merge


# In[12]:


df_merge.to_csv(r'C:\Users\Usuario\Desktop\CRAG\mapping2.csv')


# In[ ]:




