#!/usr/bin/env python
# coding: utf-8

# In[2]:


#install the package
pip install pypdb


# In[2]:


#alternatively install the development version from github: 
pip install git+git://github.com/williamgilpin/pypdb


# In[3]:


get_ipython().run_line_magic('pylab', 'inline')
from IPython.display import HTML

from pypdb import *

get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')


# In[4]:


#query database for x-ray structures
XRAY = Query('X-RAY DIFFRACTION', query_type='ExpTypeQuery').search()


# In[5]:


#not sure if these are needed but I downloaded just in case they were required for the new query below
from pypdb.clients.search.search_client import perform_search
from pypdb.clients.search.search_client import ReturnType
from pypdb.clients.search.operators import text_operators


# In[17]:


#filter the found X-ray structures for those under 4 angstroms
for IDs in XRAY:
    search_operator = text_operators.ComparisonOperator(
           value=4,
           attribute="rcsb_entry_info.resolution_combined",
           comparison_type=text_operators.ComparisonType.LESS)
    return_type = ReturnType.ENTRY

    results = perform_search(search_operator, return_type)
    
    #write to a txt file
    with open(r'C:\Users\cozyb\OneDrive\Documents\Grad School\PlantsPython\found_IDs.txt', 'w') as fp:
        for IDs in found_IDs:
            fp.write("%s\n" % IDs)
            print('Done')
    break
    


# In[ ]:





# In[ ]:





# In[ ]:




