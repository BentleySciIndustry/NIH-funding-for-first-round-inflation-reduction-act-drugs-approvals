#!/usr/bin/env python
# coding: utf-8

# In[1]:


####### 10/31/23 last update by ED ###################
#NOTE: 
    #Folders must be on the same folder layer as code
    
    #function_data = RePORTER local link for grant data on PMID 1999-2020
    
    #search_term_inputs >> BrandName_Linker.csv = your SearchID and Brand_Name cols (KEEP AS 2 COL, So Brand_Name col not always unique)
    
    #search_term_inputs >> SearchTerm_Drug and SearchTerm_TargetOnly.csv = Input / SearchID / TIME_CUT
        #Input = PubMed Search term of intrest (DO NOT ADD "NOT REVIEW[pt] or TIME RANGEs", code adds them later)
        #SearchID = This should be the Drug## or target number of the inputs (or whatever you want them to be, DO NOT LEAVE BLANK)
            # for step 2 to work, Drug ID must start with 'drug' , and target ID must be all numbers (234 for example)
            # if step 2 is not needed, any unique ID will work
        #TIME_CUT = if you want to cut PMID download for each term at a set year (2021 = full capture ) 1980 > TIME_CUT
            #Year of APY and PMID pub also part of final output, so this can be done later
            #set at for year in range(1980, 2021) > auto add to search term at 5 download year intervals (bypass 10k per cycle cap) for full analysis


#NOTE: csv / os / time / pandas should be a part of most python interface tools pre-load list
#pip install biopython
#pip install beautifulsoup4
#pip install Bio


from Bio import Entrez
from Bio.Entrez import efetch
from bs4 import BeautifulSoup

import csv   
import os
import time
start_time = time.time()
import pandas as pd


# In[2]:


####### Step 1
## CORE PMID DOWNLOAD STEPS (Have Search term inputs folder populated as seen in example format / ensure all other dependent folders are in place)


# In[4]:


######### PubYear CUT code set (use TIME_CUT col in search term CSV to set year for year stop like approval years)

def pub_key_entz_target(Search_List):
    Entrez.email = "ezhou@bentley.edu"  # USE YOUR EMAIL IF YOU'RE USING THIS CODE!!!!

    # Initialize the output CSV with headers, if necessary
    if not os.path.exists('output.csv'):
        with open('output.csv', 'w') as f:
            f.write("PMID,search_term,Search_ID,Modified_Search_Term,TIME_CUT\n")

    for index, row in Search_List.iterrows():
        search_term = row['search_term']
        time_cut = int(row['TIME_CUT']) if not pd.isna(row['TIME_CUT']) else 2021
        search_id = row['SearchID']

        # Loop from 1980 to the time_cut value with a 1-year interval (bypass pubmed limited of 10k unique PMID per download request)
        for year in range(1980, time_cut+1):
            modified_search_term = f"{search_term} AND ({year}:{year} [pdat])"

            handle = Entrez.esearch(db="pubmed", term=modified_search_term, retmax="1500000")
            record = Entrez.read(handle)
            entz = record["IdList"]

            search_output = pd.DataFrame(entz, columns=['PMID'])
            search_output['search_term'] = search_term
            search_output['Search_ID'] = search_id
            search_output['Modified_Search_Term'] = modified_search_term    
            search_output['TIME_CUT'] = '' if pd.isna(row['TIME_CUT']) else row['TIME_CUT']

            search_output.to_csv('output.csv', mode='a', index=None, header=False)

            time.sleep(0.9)  # Remove this if your input list is small. Too many requests/min can result in auto timeout


# In[5]:


SearchTerm_TargetOnly=pd.read_csv('search_term_inputs/SearchTerm_TargetOnly.csv').astype(str)
SearchTerm_TargetOnly = SearchTerm_TargetOnly.rename(columns={'ID': 'SearchID'})
SearchTerm_TargetOnly=SearchTerm_TargetOnly[['SearchID','TIME_CUT']]
APP_YEAR_Target=SearchTerm_TargetOnly.drop_duplicates()



SearchTerm_Drug=pd.read_csv('search_term_inputs/SearchTerm_Drug.csv').astype(str)
SearchTerm_Drug = SearchTerm_Drug.rename(columns={'ID': 'SearchID'})
SearchTerm_Drug=SearchTerm_Drug[['SearchID','TIME_CUT']]
APP_YEAR_DRUG=SearchTerm_Drug.drop_duplicates()


APP_YEAR_Search=APP_YEAR_Target.merge(APP_YEAR_DRUG,how='outer')
APP_YEAR_Search=APP_YEAR_Search.dropna()
APP_YEAR_Search.to_csv('new_data/Search_ID_Years.csv', index=False)


# In[7]:


################ Target PMID Download +ID LINK +TIMECUT  #############################
#DO NOT ADD "AND (1998:2020 [pdat])" to search term, 
inputkey='Target'

mod=pd.read_csv('search_term_inputs/SearchTerm_TargetOnly.csv')

mod['A']="("
mod['C']=")"

mod['search_term']=mod['A']+mod['Input']+mod['C']

mod_fix=mod[['search_term','Input','TIME_CUT','SearchID']]
mod_fix.to_csv('B_'+inputkey+'.csv', index=None)
Term_list='B_'+inputkey+'.csv'

target_list=pd.read_csv(Term_list)


# In[9]:


pub_key_entz_target(target_list)


# In[10]:


#headerList = ['PMID','search_term','SearchID','Modified_Search_Term','TIME_CUT']

w0=pd.read_csv('output.csv')
w0=w0.drop_duplicates()
w0.to_csv('new_data/P_'+Term_list, index=False)


# In[11]:


w1=pd.read_csv('new_data/P_'+Term_list)

ID_CORE=w1[['PMID','Search_ID']]
ID_CORE=ID_CORE.drop_duplicates()
ID_CORE = ID_CORE.rename(columns={'Search_ID': 'Target_ID'})
ID_CORE.to_csv('new_data/ID_CORE_Target.csv', index=False)

os.remove('output.csv')
w1.to_csv('new_data/P_'+Term_list,index=False)
w1['PMID'].nunique()


# In[12]:


################ Drug PMID Download  #############################
#DO NOT ADD "AND (1998:2020 [pdat])" to search term, 
inputkey='Drug'

mod=pd.read_csv(('search_term_inputs/SearchTerm_Drug.csv'))
#mod=pd.read_csv('Full_input_drug_only.csv')

mod['A']="("
mod['C']=")"

mod['search_term']=mod['A']+mod['Input']+mod['C']

mod_fix=mod[['search_term','Input','TIME_CUT','SearchID']]
#mod_fix=mod[['search_term']]
mod_fix.to_csv('B_'+inputkey+'.csv',index=None)
Term_list='B_'+inputkey+'.csv'


target_list=pd.read_csv(Term_list)


# In[14]:


pub_key_entz_target(target_list)


# In[15]:


w0=pd.read_csv('output.csv')
w0=w0.drop_duplicates()
w0.to_csv('new_data/P_'+Term_list, index=False)


# In[16]:


w1=pd.read_csv('new_data/P_'+Term_list)

ID_CORE=w1[['PMID','Search_ID']]
ID_CORE=ID_CORE.drop_duplicates()
ID_CORE = ID_CORE.rename(columns={'Search_ID': 'Drug_ID'})
ID_CORE.to_csv('new_data/ID_CORE_Drug.csv', index=False)


w1.to_csv('new_data/P_'+Term_list,index=False)

w1['PMID'].nunique()


# In[17]:


ID_CORE_Total_Drug=pd.read_csv('new_data/ID_CORE_Target.csv')
ID_CORE_Total_Target=pd.read_csv('new_data/ID_CORE_Drug.csv')

ID_CORE_Total=ID_CORE_Total_Drug.merge(ID_CORE_Total_Target, how='outer')
ID_CORE_Total=ID_CORE_Total.drop_duplicates()
ID_CORE_Total=ID_CORE_Total.fillna(0)

ID_CORE_Total.to_csv('new_data/PMID_Search_ID.csv', index=False)


# In[18]:


print("--- %s min ---" % ((time.time() - start_time)/60))


# In[19]:


os.remove('B_Target.csv')
os.remove('B_Drug.csv')

os.remove('output.csv')

os.remove('new_data/ID_CORE_Drug.csv')
os.remove('new_data/ID_CORE_Target.csv')


# In[2]:


####### Step 2
##################### Grant Funding Link Code (Start run here if working with local P_B_files from step 1)
import csv   
import os
import time
start_time = time.time()
import pandas as pd


# In[21]:


import datetime
timestamp = pd.Timestamp(datetime.datetime(1994, 10, 10))
TIME      = timestamp.today()


# In[22]:


def NIH_Search_Drug(search_term_drug,search_term_target):

###################### Ref Table Load ##################    
    time_short = pd.read_csv('function_data/Reporter_pub_time_1980_2020_compact.csv',dtype=str)
    core_set=pd.read_csv('function_data/Reporter_project_1985_2020_compact_sum.csv')
    
    core_set['TOTAL_COST'] = pd.to_numeric(core_set['TOTAL_COST'], errors='coerce').fillna(0)

    
    pub_key = pd.read_csv('function_data/Reporter_Link_Table_1980_2021_compact.csv',dtype=str)
    
    #Inflation_key_V2018 = pd.read_csv('function_data/inf2018_key.csv')  ### DO not use if Step 2 is also used (inflation change in Step 2 code too)
    
    Inflation_key_V2018 = pd.read_csv('function_data/inf_None_key.csv') #### code to not inflation adjust any years
    
    IC_Key = pd.read_csv('function_data/institute_key.csv',dtype=str)
    
###################### Inflation ##################
    #core_set = core_set.loc[core_set['FY'] >=1985] 1985-2020 dataset loaded: 09/09/23
    core_set_2018=pd.merge(core_set,Inflation_key_V2018, how='outer')
    core_set_2018['inf_2018_costs']=core_set_2018['TOTAL_COST']*core_set_2018['inf_2018']
    core_set_2018.drop('TOTAL_COST', inplace=True, axis=1)
    core_set_2018.drop('inf_2018', inplace=True, axis=1)
    core_set_2018.rename(columns = {"inf_2018_costs": "TOTAL_COST"}, 
              inplace = True)
    core_set=core_set_2018
    core_set_2018=1
################################################# Drug vs Target PMID ############################### 
############### Drug #########
    read_pmid = pd.read_csv('new_data/P_B_Drug.csv')
    read_pmid_pmid_UQ=read_pmid[['PMID']]
    read_pmid_pmid_UQ=read_pmid_pmid_UQ.drop_duplicates()
    read_pmid_pmid_UQ['Data_Type']='Drug'
    read_pmid_pmid_UQ.to_csv('new_data/pmid_drug.csv',index=None)
    J_Drug_pmid_UQ = read_pmid_pmid_UQ
    #os.remove("new_data/P_B_Drug.csv")
    
###################### Target ##################    
    read_pmid = pd.read_csv('new_data/P_B_Target.csv')
    read_pmid_pmid_UQ=read_pmid[['PMID']]
    read_pmid_pmid_UQ=read_pmid_pmid_UQ.drop_duplicates()
    read_pmid_pmid_UQ['Data_Type']='Target'
    read_pmid_pmid_UQ.to_csv('new_data/pmid_target.csv',index=None)
    J_Target_pmid_UQ=read_pmid_pmid_UQ
    #os.remove("new_data/P_B_Target.csv")
    
##################### Drug \ Target  ##########################################    
    Drug_hold=pd.read_csv('new_data/pmid_drug.csv')
    Target_hold=pd.read_csv('new_data/pmid_target.csv')
    
    Drug_hold_full=Drug_hold
    Target_hold_full=Target_hold
    Drug_hold=Drug_hold[['PMID']]
    Target_hold=Target_hold[['PMID']]
    Target_hold['search']='Target'
    Drug_hold['search2']='Drug'
    Drug_UQ_FULL=pd.merge(Drug_hold,Target_hold, how='outer')
    Drug_UQ_FULL=Drug_UQ_FULL.fillna('only')
    Drug_UQ_FULL['search_term']=Drug_UQ_FULL['search2']+Drug_UQ_FULL['search']
    Drug_UQ_FULL_join=Drug_UQ_FULL[['PMID','search_term']]
    Drug_UQ_FULL_join=Drug_UQ_FULL_join.drop_duplicates()
    result = Drug_hold_full.append(Target_hold_full)
    Type_Fix_PMID=pd.merge(result,Drug_UQ_FULL_join, how='outer')
    Type_Fix_PMID=Type_Fix_PMID.drop_duplicates()
    Type_Fix_PMID['search_term'] = Type_Fix_PMID['search_term'].replace({'DrugTarget':'Drug', 'Drugonly':'Drug'})
    Type_Fix_PMID['search_term'] = Type_Fix_PMID['search_term'].replace({'onlyTarget':'TargetOnly', 'Drugonly':'Drug'})
    Type_Fix_Drug=Type_Fix_PMID.loc[Type_Fix_PMID['search_term'] == 'Drug']
    Type_Fix_TargetOnly=Type_Fix_PMID.loc[Type_Fix_PMID['search_term'] == 'TargetOnly']
    
###################### Drug NIH Funding ###########################################################################
     
    NIH_PMID=pd.merge(Type_Fix_Drug.astype(str),pub_key.astype(str), how='inner')
    NIH_PMID=NIH_PMID.drop_duplicates()
    NIH_PMID_time=pd.merge(NIH_PMID.astype(str),time_short.astype(str), how='inner')
    NIH_PMID_time_core=pd.merge(core_set,NIH_PMID_time, how='inner')
    NIH_PMID_time_core=NIH_PMID_time_core.fillna(0)
    NIH_PMID_time_core=NIH_PMID_time_core.drop_duplicates()
    
###################### Pubyear_fixer    
    Type_Fix_Drug_pubyear=NIH_PMID_time_core[['PMID','PUB_YEAR']]
    Type_Fix_Drug_pubyear=Type_Fix_Drug_pubyear.drop_duplicates()
    PB_d = Type_Fix_Drug_pubyear.groupby('PMID')['PUB_YEAR'].max().reset_index()
    PB_d.columns = ['PMID', 'PUB_YEAR']
    NIH_PMID_time_core.drop('PUB_YEAR', inplace=True, axis=1)
    NIH_PMID_time_core=NIH_PMID_time_core.drop_duplicates()
    #PB_d.to_csv('Pub_Year_SET1.csv')
    NIH_PMID_time_core=NIH_PMID_time_core.merge(PB_d, how='inner')
    
    nih_test_start=NIH_PMID_time_core
    nih_test_Last=NIH_PMID_time_core
    nih_test1=NIH_PMID_time_core
    
###################### NIH_Funding_APY ##################  
    nih_test_start2 = nih_test_start[nih_test_start.groupby('PROJECT_NUMBER')['FY'].transform('min') == nih_test_start['FY']]
    nih_test_start2 = nih_test_start2[['PROJECT_NUMBER','FY']]
    nih_test_start2=nih_test_start2.drop_duplicates()
    nih_test_start2.rename(columns = {"FY": "FY_Start"}, 
              inplace = True)
    nih_test_Last2 = nih_test_Last[nih_test_Last.groupby('PROJECT_NUMBER')['FY'].transform('max') == nih_test_start['FY']]
    nih_test_Last2 = nih_test_Last2[['PROJECT_NUMBER','FY']]
    nih_test_Last2=nih_test_Last2.drop_duplicates()
    nih_test_Last2.rename(columns = {"FY": "FY_Last"}, 
              inplace = True)

    nih_target_home = pd.merge(nih_test1, nih_test_start2,  how='left', left_on=['PROJECT_NUMBER'], right_on = ['PROJECT_NUMBER'])
    nih_target_home=nih_target_home.drop_duplicates()
    nih_target_home = pd.merge(nih_target_home, nih_test_Last2,  how='left', left_on=['PROJECT_NUMBER'], right_on = ['PROJECT_NUMBER'])
    nih_target_home=nih_target_home.drop_duplicates()
    
    nih_target_home_sl=nih_target_home[['PROJECT_NUMBER','FY_Start','FY_Last','PMID','PUB_YEAR','TOTAL_COST','FY']]
    nih_target_home_sl['FY_Last'] = nih_target_home_sl['FY_Last'].astype(int)
    nih_target_home_sl['FY_Start'] = nih_target_home_sl['FY_Start'].astype(int)
    nih_target_home_sl['PUB_YEAR'] = nih_target_home_sl['PUB_YEAR'].astype(int)
    nih_target_home_sl['Pub_V_Start']=nih_target_home_sl['PUB_YEAR'] - nih_target_home_sl['FY_Start']
    nih_target_home_sl['Pub_V_Last']=nih_target_home_sl['PUB_YEAR'] - nih_target_home_sl['FY_Last']
    nih_target_home_sl = nih_target_home_sl.loc[nih_target_home_sl['Pub_V_Last'] <=4]
    nih_target_home_sl = nih_target_home_sl.loc[nih_target_home_sl['Pub_V_Start'] >=0]
    nih_target_home_sl.drop('Pub_V_Last', inplace=True, axis=1)
    nih_target_home_sl.drop('Pub_V_Start', inplace=True, axis=1)
    
    APY_Cost_index=nih_target_home_sl[['PROJECT_NUMBER','TOTAL_COST','FY']]
    APY_Cost_index=APY_Cost_index.drop_duplicates()
    APY_Cost_index['TOTAL_COST'] = APY_Cost_index['TOTAL_COST'].astype(int)
    g = APY_Cost_index.groupby(['PROJECT_NUMBER','FY'])['TOTAL_COST'].sum()
    j = APY_Cost_index.groupby(['PROJECT_NUMBER','FY']).size().to_frame('count')
    NIH_base=pd.merge(g, j, left_index=True, right_index=True).reset_index()
    NIH_base.rename(columns = {"FY": "APY"}, 
              inplace = True)
    
    nih_test2=nih_target_home_sl
    nih_test2['PUB_YEAR'] = nih_test2['PUB_YEAR'].astype(int)
    nih_test2['FY_Last'] = nih_test2['FY_Last'].astype(int)
    nih_test2['APY'] = nih_test2[['PUB_YEAR','FY_Last']].min(axis=1)
    nih_test2['APY'] = nih_test2['APY'].astype(str)
    nih_test2['PROJECT_NUMBER'] = nih_test2['PROJECT_NUMBER'].astype(str)
    nih_test2["ACTUAL_PROJECT_YEAR"] = nih_test2["APY"] + nih_test2["PROJECT_NUMBER"]
    nih_test2.drop('TOTAL_COST', inplace=True, axis=1)
    nih_test2.drop('FY', inplace=True, axis=1)

    NIH_base['APY'] = NIH_base['APY'].astype(str)
    nih_test2=pd.merge(nih_test2,NIH_base, how='inner')
    nih_test2=nih_test2.drop_duplicates()  
    
    nih_test2.rename(columns = {"count": "Project_Count"}, 
          inplace = True)
    nih_test2.rename(columns = {"TOTAL_COST": "APY_COST_inf2018"}, 
              inplace = True)
    nih_test2=nih_test2.assign(Activity_Code=nih_test2['PROJECT_NUMBER'].str[:3])
    nih_test2=nih_test2.assign(Institute_Code=nih_test2['PROJECT_NUMBER'].str[3:5])
    nih_test2=pd.merge(nih_test2,IC_Key, how='inner')
    nih_test2=nih_test2.drop_duplicates()
    nih_test2['Search__ID']=search_term_drug
    nih_test2['Search_Type']='Drug'

    nih_test2=nih_test2[['Search__ID','Search_Type','PMID','PUB_YEAR','PROJECT_NUMBER','FY_Start','FY_Last','APY','ACTUAL_PROJECT_YEAR','APY_COST_inf2018','Activity_Code','Institute_Code','Acronym_institute_name','full_institute_name','Compressed Names','Project_Count']]
    nih_test2.to_csv('new_data/Drug_hold.csv',index=None)
    J_Drug_hold=nih_test2
    
    
###################### Target NIH Funding ########################################################################### 
     
    NIH_PMID=pd.merge(Type_Fix_TargetOnly.astype(str),pub_key.astype(str), how='inner')
    NIH_PMID=NIH_PMID.drop_duplicates()
    NIH_PMID_time=pd.merge(NIH_PMID.astype(str),time_short.astype(str), how='inner')
    NIH_PMID_time_core=pd.merge(core_set,NIH_PMID_time, how='inner')
    NIH_PMID_time_core=NIH_PMID_time_core.fillna(0)
    NIH_PMID_time_core=NIH_PMID_time_core.drop_duplicates()
    
###################### Pubyear_fixer    
    Type_Fix_Drug_pubyear=NIH_PMID_time_core[['PMID','PUB_YEAR']]
    Type_Fix_Drug_pubyear=Type_Fix_Drug_pubyear.drop_duplicates()
    PB_d = Type_Fix_Drug_pubyear.groupby('PMID')['PUB_YEAR'].max().reset_index()
    PB_d.columns = ['PMID', 'PUB_YEAR']
    NIH_PMID_time_core.drop('PUB_YEAR', inplace=True, axis=1)
    NIH_PMID_time_core=NIH_PMID_time_core.drop_duplicates()
    #PB_d.to_csv('Pub_Year_SET2.csv')
    NIH_PMID_time_core=NIH_PMID_time_core.merge(PB_d, how='inner')
    
    nih_test_start=NIH_PMID_time_core
    nih_test_Last=NIH_PMID_time_core
    nih_test1=NIH_PMID_time_core
    
###################### NIH_Funding_APY ##################  
    nih_test_start2 = nih_test_start[nih_test_start.groupby('PROJECT_NUMBER')['FY'].transform('min') == nih_test_start['FY']]
    nih_test_start2 = nih_test_start2[['PROJECT_NUMBER','FY']]
    nih_test_start2=nih_test_start2.drop_duplicates()
    nih_test_start2.rename(columns = {"FY": "FY_Start"}, 
              inplace = True)
    nih_test_Last2 = nih_test_Last[nih_test_Last.groupby('PROJECT_NUMBER')['FY'].transform('max') == nih_test_start['FY']]
    nih_test_Last2 = nih_test_Last2[['PROJECT_NUMBER','FY']]
    nih_test_Last2=nih_test_Last2.drop_duplicates()
    nih_test_Last2.rename(columns = {"FY": "FY_Last"}, 
              inplace = True)

    nih_target_home = pd.merge(nih_test1, nih_test_start2,  how='left', left_on=['PROJECT_NUMBER'], right_on = ['PROJECT_NUMBER'])
    nih_target_home=nih_target_home.drop_duplicates()
    nih_target_home = pd.merge(nih_target_home, nih_test_Last2,  how='left', left_on=['PROJECT_NUMBER'], right_on = ['PROJECT_NUMBER'])
    nih_target_home=nih_target_home.drop_duplicates()
    
    nih_target_home_sl=nih_target_home[['PROJECT_NUMBER','FY_Start','FY_Last','PMID','PUB_YEAR','TOTAL_COST','FY']]
    nih_target_home_sl['FY_Last'] = nih_target_home_sl['FY_Last'].astype(int)
    nih_target_home_sl['FY_Start'] = nih_target_home_sl['FY_Start'].astype(int)
    nih_target_home_sl['PUB_YEAR'] = nih_target_home_sl['PUB_YEAR'].astype(int)
    nih_target_home_sl['Pub_V_Start']=nih_target_home_sl['PUB_YEAR'] - nih_target_home_sl['FY_Start']
    nih_target_home_sl['Pub_V_Last']=nih_target_home_sl['PUB_YEAR'] - nih_target_home_sl['FY_Last']
    nih_target_home_sl = nih_target_home_sl.loc[nih_target_home_sl['Pub_V_Last'] <=4]
    nih_target_home_sl = nih_target_home_sl.loc[nih_target_home_sl['Pub_V_Start'] >=0]
    nih_target_home_sl.drop('Pub_V_Last', inplace=True, axis=1)
    nih_target_home_sl.drop('Pub_V_Start', inplace=True, axis=1)
    
    APY_Cost_index=nih_target_home_sl[['PROJECT_NUMBER','TOTAL_COST','FY']]
    APY_Cost_index=APY_Cost_index.drop_duplicates()
    APY_Cost_index['TOTAL_COST'] = APY_Cost_index['TOTAL_COST'].astype(int)
    g = APY_Cost_index.groupby(['PROJECT_NUMBER','FY'])['TOTAL_COST'].sum()
    j = APY_Cost_index.groupby(['PROJECT_NUMBER','FY']).size().to_frame('count')
    NIH_base=pd.merge(g, j, left_index=True, right_index=True).reset_index()
    NIH_base.rename(columns = {"FY": "APY"}, 
              inplace = True)
    
    nih_test2=nih_target_home_sl
    nih_test2['PUB_YEAR'] = nih_test2['PUB_YEAR'].astype(int)
    nih_test2['FY_Last'] = nih_test2['FY_Last'].astype(int)
    nih_test2['APY'] = nih_test2[['PUB_YEAR','FY_Last']].min(axis=1)
    nih_test2['APY'] = nih_test2['APY'].astype(str)
    nih_test2['PROJECT_NUMBER'] = nih_test2['PROJECT_NUMBER'].astype(str)
    nih_test2["ACTUAL_PROJECT_YEAR"] = nih_test2["APY"] + nih_test2["PROJECT_NUMBER"]
    nih_test2.drop('TOTAL_COST', inplace=True, axis=1)
    nih_test2.drop('FY', inplace=True, axis=1)

    NIH_base['APY'] = NIH_base['APY'].astype(str)
    nih_test2=pd.merge(nih_test2,NIH_base, how='inner')
    nih_test2=nih_test2.drop_duplicates()  
    
    nih_test2.rename(columns = {"count": "Project_Count"}, 
          inplace = True)
    nih_test2.rename(columns = {"TOTAL_COST": "APY_COST_inf2018"}, 
              inplace = True)
    nih_test2=nih_test2.assign(Activity_Code=nih_test2['PROJECT_NUMBER'].str[:3])
    nih_test2=nih_test2.assign(Institute_Code=nih_test2['PROJECT_NUMBER'].str[3:5])
    nih_test2=pd.merge(nih_test2,IC_Key, how='inner')
    nih_test2=nih_test2.drop_duplicates()
    nih_test2['Search__ID']=search_term_target
    nih_test2['Search_Type']='target'

    nih_test2=nih_test2[['Search__ID','Search_Type','PMID','PUB_YEAR','PROJECT_NUMBER','FY_Start','FY_Last','APY','ACTUAL_PROJECT_YEAR','APY_COST_inf2018','Activity_Code','Institute_Code','Acronym_institute_name','full_institute_name','Compressed Names','Project_Count']]
    nih_test2.to_csv('new_data/Target_hold.csv',index=None)
    J_Target_hold=nih_test2
    
###########################Drug vs Target Analysis########################################################

#################### Search Terms ############################
    Drug_hold=pd.read_csv('new_data/Drug_hold.csv')
    Target_hold=pd.read_csv('new_data/Target_hold.csv')
    
    Drug_hold_full=Drug_hold
    Target_hold_full=Target_hold
    Drug_hold=Drug_hold[['ACTUAL_PROJECT_YEAR']]
    Target_hold=Target_hold[['ACTUAL_PROJECT_YEAR']]
    Target_hold['search']='Target'
    Drug_hold['search2']='Drug'
    Drug_UQ_FULL=pd.merge(Drug_hold,Target_hold, how='outer')
    Drug_UQ_FULL=Drug_UQ_FULL.fillna('only')
    Drug_UQ_FULL['search_term']=Drug_UQ_FULL['search2']+Drug_UQ_FULL['search']
    Drug_UQ_FULL_join=Drug_UQ_FULL[['ACTUAL_PROJECT_YEAR','search_term']]
    Drug_UQ_FULL_join=Drug_UQ_FULL_join.drop_duplicates()
    result = Drug_hold_full.append(Target_hold_full)
    resultUQ_FULL=pd.merge(result,Drug_UQ_FULL_join, how='outer')
    resultUQ_FULL=resultUQ_FULL.drop_duplicates()
    resultUQ_FULL['search_term'] = resultUQ_FULL['search_term'].replace({'DrugTarget':'Drug', 'Drugonly':'Drug'})
    resultUQ_FULL['search_term'] = resultUQ_FULL['search_term'].replace({'onlyTarget':'TargetOnly', 'Drugonly':'Drug'})

############# Project/Grant ID ##############################
    Grant_code=pd.read_csv('function_data/Grant_Types.csv')
    Grant_code=Grant_code.assign(Activity_Code=Grant_code['Activity_Code'].str[:3])
    resultUQ_GT=pd.merge(resultUQ_FULL,Grant_code,how='outer')
    resultUQ_GT = resultUQ_GT.dropna(axis=0, subset=['Search_Type'])
    resultUQ_GT["Grant_Type_Name"][resultUQ_GT['Activity_Code'].str.contains("Z")] = "Intramural Programs"
    resultUQ_GT=resultUQ_GT.fillna('Others')
    resultUQ_GT=resultUQ_GT.drop_duplicates()
    resultUQ_FULL=resultUQ_GT
    
    
############# Cost per APY ##################################    
    resultUQ_FULL=resultUQ_FULL[['Search__ID','Search_Type','search_term','PMID','PUB_YEAR','PROJECT_NUMBER','FY_Start','FY_Last','APY','ACTUAL_PROJECT_YEAR','APY_COST_inf2018','Activity_Code','Institute_Code','Acronym_institute_name','full_institute_name','Compressed Names','Project_Count','Grant_Type_Name']]
    resultUQ_FULL.rename(columns = {"Search_Type": "Source_Search_Type"}, 
          inplace = True)
    resultUQ_FULL.rename(columns = {"search_term": "Search_Type"}, 
          inplace = True)

    resultUQ_FULL.to_csv('new_data/resultUQ_FULL.csv',index=None)
    resultUQ_APY_COST=resultUQ_FULL[['APY','ACTUAL_PROJECT_YEAR','APY_COST_inf2018','Search_Type']]
    resultUQ_APY_COST=resultUQ_APY_COST.drop_duplicates()
    resultUQ_APY_COST.to_csv('new_data/resultUQ_APY_COST.csv',index=None)
    J_resultUQ_APY_COST=resultUQ_APY_COST
    J_resultUQ_FULL=resultUQ_FULL


    
########### Debug Cleanup  ####################    
    core_set=1
    time = 1
    key = 1
    NIH_PMID=1
    time_short=1
    J_resultUQ_FULL=1
    J_resultUQ_APY_COST=1
    resultUQ_APY_COST=1
    resultUQ_FULL=1
    Drug_UQ_FULL=1
    Target_hold=1
    resultUQ_GT=1
    Grant_code=1
    result=1


# In[ ]:


NIH_Search_Drug('APY_drug','APY_target_only')


# In[24]:


print("--- %s min ---" % ((time.time() - start_time)/60))


# In[ ]:


print('#0################################################################Input_PMID#')
Grant_code=pd.read_csv('new_data/resultUQ_APY_COST.csv')
Drug=Grant_code.loc[Grant_code['Search_Type'] == 'Drug']
Target=Grant_code.loc[Grant_code['Search_Type'] == 'TargetOnly']
Grant_code = Grant_code[['ACTUAL_PROJECT_YEAR','APY_COST_inf2018']]

Grant_code_FULL=pd.read_csv('new_data/resultUQ_FULL.csv')
Drug_PMID=Grant_code_FULL.loc[Grant_code_FULL['Search_Type'] == 'Drug']
Drug_PMID=Drug_PMID[['PMID']]
Target_PMID=Grant_code_FULL.loc[Grant_code_FULL['Search_Type'] == 'TargetOnly']
Target_PMID=Target_PMID[['PMID']]
Grant_code_FULL=Grant_code_FULL[['PMID']]

Drug = Drug.drop_duplicates()
Drug_PMID = Drug_PMID.drop_duplicates()
Target_PMID = Target_PMID.drop_duplicates()
Target = Target.drop_duplicates()
Grant_code_FULL=Grant_code_FULL.drop_duplicates()

PMID_DRUG=pd.read_csv('new_data/pmid_drug.csv')
PMID_TARGET=pd.read_csv('new_data/pmid_target.csv')
PMID_DRUG = PMID_DRUG.drop_duplicates()
PMID_TARGET = PMID_TARGET.drop_duplicates()
PMID_Full=pd.merge(PMID_DRUG,PMID_TARGET, how='outer')
PMID_Full=PMID_Full.drop_duplicates()

################  PMID  ###################################

Drug_hold=pd.read_csv('new_data/pmid_drug.csv')
Drug_hold=Drug_hold[['PMID']]
Target_hold=pd.read_csv('new_data/pmid_target.csv')
Target_hold=Target_hold[['PMID']]
SET = pd.merge(Drug_hold, Target_hold, how='outer', indicator=True)
SET=SET[['PMID','_merge']]
SET = SET.drop_duplicates()

Target1 = SET.loc[SET._merge == 'right_only', ['PMID']]
Target1=Target1[['PMID']]
Target1 = Target1.drop_duplicates()

Drug1 = SET.loc[SET._merge == 'left_only', ['PMID']]
Drug1=Drug1[['PMID']]
Drug2 = SET.loc[SET._merge == 'both', ['PMID']]
Drug2=Drug2[['PMID']]
Drug_total = pd.merge(Drug1, Drug2, how='outer')
Drug_total = Drug_total.drop_duplicates()

Drug_Start_PMID_number=Drug_total.nunique()
print(('Drug_PMID_number::::') + str(Drug_Start_PMID_number))

Target_Start_PMID_number=Target1.nunique()
print(('Target_PMID_number::::') + str(Target_Start_PMID_number))
print('#1##############################################################Funded_PMID#')

################ Funded PMID  ###################################

Drug_hold=pd.read_csv('new_data/Drug_hold.csv')
Drug_hold=Drug_hold[['PMID']]
Target_hold=pd.read_csv('new_data/Target_hold.csv')
Target_hold=Target_hold[['PMID']]
SET = pd.merge(Drug_hold, Target_hold, how='outer', indicator=True)
SET=SET[['PMID','_merge']]
SET = SET.drop_duplicates()

Target1 = SET.loc[SET._merge == 'right_only', ['PMID']]
Target1=Target1[['PMID']]
Target1 = Target1.drop_duplicates()

Drug1 = SET.loc[SET._merge == 'left_only', ['PMID']]
Drug1=Drug1[['PMID']]
Drug2 = SET.loc[SET._merge == 'both', ['PMID']]
Drug2=Drug2[['PMID']]
Drug_total = pd.merge(Drug1, Drug2, how='outer')
Drug_total = Drug_total.drop_duplicates()

Drug_PMID_number=Drug_total.nunique()
print(('Drug_PMID_number::::') + str(Drug_PMID_number))

Target_PMID_number=Target1.nunique()
print(('Target_PMID_number::::') + str(Target_PMID_number))
print('#2#####################################################################APYs#')

################ APY  ##################################################

Target_APY=Target['ACTUAL_PROJECT_YEAR'].nunique()
print(('Target_APY::::') + str(Target_APY))

Drug_APY=Drug['ACTUAL_PROJECT_YEAR'].nunique()
print(('Drug_APY::::') + str(Drug_APY))
print('#3#########################################################NIH_Funding_Cost#')

############# Cost of funding ##########################################################

Target_cost=Target['APY_COST_inf2018'].sum(axis = 0, skipna = True)
print(('Target_cost::::') + str(Target_cost))

Drug_cost=Drug['APY_COST_inf2018'].sum(axis = 0, skipna = True)
print(('Drug_cost::::') + str(Drug_cost))
print('#4####################################################Total_Unique_analysis#')

################## Total  ####################

PMID_Full=PMID_Full[['PMID']]
PMID_Full.astype(int)

PMID_Full_number=PMID_Full.nunique()
print(('Grant_code_ALL_PMID::::') + str(PMID_Full_number))

Grant_code_FULL=Grant_code_FULL.nunique()
print(('Grant_code_Funded_PMID::::') + str(Grant_code_FULL))

Grant_code_APY=Grant_code['ACTUAL_PROJECT_YEAR'].nunique()
print(('Grant_code_APY::::') + str(Grant_code_APY))

Grant_code_cost=Grant_code['APY_COST_inf2018'].sum(axis = 0, skipna = True)
print(('Grant_code_cost::::') + str(Grant_code_cost))

Overview=pd.DataFrame(([['Drug_Input_PMID', int(Drug_Start_PMID_number)], 
                  ['Target_Input_PMID', int(Target_Start_PMID_number)], 
                  ['Drug_Funded_PMID', int(Drug_PMID_number)],
                  ['Target_Funded_PMID', int(Target_PMID_number)], 
                  ['Target_APY', Target_APY],
                  ['Drug_APY', Drug_APY], 
                  ['Target_cost', Target_cost],
                  ['Drug_cost', Drug_cost], 
                  ['Unique_PMID', int(PMID_Full_number)],
                  ['Unique_Funded_PMID', int(Grant_code_FULL)], 
                  ['Unique_APY', Grant_code_APY],
                  ['Unique_Cost', Grant_code_cost],
                  ['Analysis_Date', TIME]]), 
                columns=['Index', 'Amount'])
Overview.to_csv('new_data/Overview.csv',index=None)


# In[ ]:




