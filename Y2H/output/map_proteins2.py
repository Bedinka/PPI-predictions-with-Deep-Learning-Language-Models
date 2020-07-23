#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd


# In[3]:


df1 = pd.read_csv ("proteins2.csv", sep=';')
print(df1)


# In[4]:


df2 = pd.read_csv ("all_interactions.csv", sep=';')
print(df2)


# In[5]:


df_merge = pd.merge(df1, df2, how='inner')


# In[6]:


df_merge


# In[7]:


df_merge.to_csv(r'C:\Users\Usuario\Desktop\CRAG\mapping2.csv')


# In[ ]:




