#!/usr/bin/env python
# coding: utf-8

# In[57]:


from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
import pandas as pd

from rdkit.Chem import Draw


# In[58]:


data = pd.read_csv('test_ssp.csv', index_col=None)
print(data)


# In[59]:


#Proff and make a list of Smiles and id
c_smiles = []
count = 0
for index, row in data.iterrows():
    try:
      cs = Chem.CanonSmiles(row['SMILES'])
      c_smiles.append([row['ID_Name'], cs])
    except:
      count = count + 1
      print('Count Invalid SMILES:', count, row['ID_Name'], row['SMILES'])

print(c_smiles)


# In[60]:


# make a list of id, smiles, and mols
ms = []
df =pd.DataFrame(c_smiles,columns=['ID_Name','SMILES'])
for index, row in df.iterrows():
    mol = Chem.MolFromSmiles(row['SMILES'])
    ms.append([row['ID_Name'], row['SMILES'], mol])

 


# In[61]:


# make a list of id, smiles, mols, and fingerprints (fp)
fps = []
df_fps = pd.DataFrame(ms,columns=['ID_Name','SMILES', 'mol'])
df_fps.head

for index, row in df_fps.iterrows():
   fps_cal = FingerprintMols.FingerprintMol(row['mol'], minPath=1, maxPath=7, fpSize=2048,
                              bitsPerHash=2, useHs=True, tgtDensity=0.0,
                              minSize=128)
   fps.append([row['ID_Name'], fps_cal])


fps_2 = pd.DataFrame(fps,columns=['ID_Name','fps'])
fps_2 = fps_2[fps_2.columns[1]]
fps_2 = fps_2.values.tolist()


# In[62]:


# the list for the dataframe
qu, ta, sim = [], [], []


# In[63]:


# compare all fp pairwise without duplicates
for n in range(len(fps_2)): 
    s = DataStructs.BulkTanimotoSimilarity(fps_2[n], fps_2[n+1:])
    for m in range(len(s)):
        qu.append(c_smiles[n])
        ta.append(c_smiles[n+1:][m])
        sim.append(s[m])
print()


# In[64]:


# build the dataframe and sort it
d = {'query':qu, 'target':ta, 'Similarity':sim}
df_final = pd.DataFrame(data=d)
df_final = df_final.sort_values('Similarity', ascending=False)
print(df_final)


# In[65]:


# save as csv
df_final.to_csv('test.csv', index=False, sep=',')


# In[ ]:




