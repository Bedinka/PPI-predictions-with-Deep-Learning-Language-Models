#!/usr/bin/env python
# coding: utf-8

# In[5]:


import pandas as pd


# In[6]:


df1 = pd.read_csv ("huri.csv", sep=';')
print(df1)


# In[31]:


df2 = pd.read_csv ("DMI.csv", sep=';')
print(df2)


# In[32]:


df_merge = pd.merge(df1, df2, how='inner')   #to obtain common values between both files


# In[33]:


df_merge


# In[34]:


df_merge.to_csv(r'C:\Users\Marina\Desktop\common_huri_dmi.csv')


# In[ ]:




