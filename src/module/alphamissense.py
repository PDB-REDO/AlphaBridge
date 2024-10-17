# -*- coding: utf-8 -*-
"""
Created on Wed Feb 7 16:53:05 2024

@author: Dan_salv
"""

import os
import sys
import json
import pandas as pd
import requests
import re


module_dir = os.path.dirname(os.path.realpath(__file__))
WEBSITE_API = 'https://alphafold.ebi.ac.uk/api/prediction'
KEY_API = 'key=AIzaSyCeurAJz7ZGjPQUtEaerUkBZ3TaBkXrY94'

class ALPHAMISSENSE:
    
    def __init__(self, UNIPROT_ID):
        
        self.UNIPROT_ID = UNIPROT_ID
    
    
    def extract_am_annotation(self):
        UNIPROT_ID = self.UNIPROT_ID
        
        r = get_url(f"{WEBSITE_API}/{UNIPROT_ID}?{KEY_API}")

        data = r.json()[0]

        am_annotation = pd.read_csv(data['amAnnotationsUrl'])
        
        return am_annotation
    
    def get_pathogenicity_list(self):
        
        am_annotation = self.extract_am_annotation()
        
        am_annotation = am_annotation.apply(lambda row:split_column(row), axis=1)
        
        pathogenicity_tuple_list = sorted([(int(group_name[1]) ,df_group['am_pathogenicity'].mean()) for group_name, df_group in am_annotation.groupby(['REF','POS'],sort = False)], key=lambda a: a[0])

        pathogenicity_list = list(zip(*pathogenicity_tuple_list))[1]
        
        return pathogenicity_list
    
    
    

def get_url(url, **kwargs):
  response = requests.get(url, **kwargs);

  if not response.ok:
    print(response.text)
    response.raise_for_status()
    sys.exit()

  return response

def split_column(row):
    
    protein_variant = row.protein_variant

    REF, POS, ALT = tuple(re.findall(r'[A-Za-z]+|\d+', protein_variant))
    
    row['REF'] = REF
    row['POS'] = POS
    row['ALT'] = ALT
    return row