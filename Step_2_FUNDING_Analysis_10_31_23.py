#!/usr/bin/env python
# coding: utf-8

# In[28]:


####### 10/30/23 last update by ED ###################
#### NOTE:  This assumes you ran Step 1 first
#### NOTE:  This is only used if you want to breakdown output by drug (this will re-add all cross search duplicates = TargetOnly/Drug true by unique brandname only)
                # Remember with this applied and basic unique totals no longer exist if cross-search TargetOnly/Drug mismatch for any shared PMID

import pandas as pd
import numpy as np
import os

Grant_code=pd.read_csv('new_data/resultUQ_FULL.csv')
Search_Year=pd.read_csv('new_data/Search_ID_Years.csv')
BRAND_ID=pd.read_csv('search_term_inputs/BrandName_Linker.csv')

ID_CORE_Total_Drug=pd.read_csv('new_data/P_B_Drug.csv')
ID_CORE_Total_Target=pd.read_csv('new_data/P_B_Target.csv')
ID_CORE_Total_Target = ID_CORE_Total_Target[(ID_CORE_Total_Target != 0).all(axis=1)]
ID_CORE_Total_Drug = ID_CORE_Total_Drug[(ID_CORE_Total_Drug != 0).all(axis=1)]
ID_CORE_Total_Drug = ID_CORE_Total_Drug.rename(columns={'SearchID': 'Search_ID'})
ID_CORE_Total_Target = ID_CORE_Total_Target.rename(columns={'SearchID': 'Search_ID'})

Search_Year = Search_Year.rename(columns={'SearchID': 'Search_ID'})

BRAND_ID = BRAND_ID.rename(columns={'SearchID': 'Search_ID'})


# In[30]:


PMID_DOWNLOADED_COUNT_Drug=ID_CORE_Total_Drug[['PMID','Search_ID','TIME_CUT']]
PMID_DOWNLOADED_COUNT_Drug=PMID_DOWNLOADED_COUNT_Drug.drop_duplicates()
PMID_DOWNLOADED_COUNT_Target=ID_CORE_Total_Target[['PMID','Search_ID','TIME_CUT']]
ID_CORE_Total_Target=ID_CORE_Total_Target.drop_duplicates()

PMID_DOWNLOADED_COUNT=pd.concat([PMID_DOWNLOADED_COUNT_Drug,PMID_DOWNLOADED_COUNT_Target])


# In[31]:


PMID_DOWNLOADED_Analysis = PMID_DOWNLOADED_COUNT.groupby('Search_ID')['PMID'].nunique().reset_index()
PMID_DOWNLOADED_Analysis.columns = ['Search_ID', 'PMID_Downloaded']


# In[32]:


BRAND_ID['Search_ID']=BRAND_ID['Search_ID'].astype(str)
PMID_DOWNLOADED_Analysis['Search_ID']=PMID_DOWNLOADED_Analysis['Search_ID'].astype(str)


# In[33]:


def check_for_drug(x):
    if 'drug' in x.lower():
        return 'Drug'
    else:
        return 'TargetOnly'
PMID_DOWNLOADED_Analysis=PMID_DOWNLOADED_Analysis.merge(BRAND_ID, how='outer')
PMID_DOWNLOADED_Analysis['Data_Type_by_Drug'] = PMID_DOWNLOADED_Analysis['Search_ID'].apply(check_for_drug)
PMID_DOWNLOADED_Analysis=PMID_DOWNLOADED_Analysis[['Brand_Name','Data_Type_by_Drug','PMID_Downloaded']]


# In[35]:


ID_CORE_Total_Drug=ID_CORE_Total_Drug[['PMID','Search_ID']]
ID_CORE_Total_Drug=ID_CORE_Total_Drug.drop_duplicates()


ID_CORE_Total_Target=ID_CORE_Total_Target[['PMID','Search_ID']]
ID_CORE_Total_Target=ID_CORE_Total_Target.drop_duplicates()


# In[36]:


def remove_decimal_from_id(id_value):
    try:
        # Convert to float, then to int, then to string
        return str(int(float(id_value)))
    except ValueError:
        # If it's not a number, just return the original string
        return id_value
ID_CORE_Total_Drug['Search_ID'] = ID_CORE_Total_Drug['Search_ID'].apply(remove_decimal_from_id)
ID_CORE_Total_Target['Search_ID'] = ID_CORE_Total_Target['Search_ID'].apply(remove_decimal_from_id)
Search_Year['Search_ID'] = Search_Year['Search_ID'].apply(remove_decimal_from_id)


# In[37]:


Grant_code_drug=ID_CORE_Total_Drug.merge(Grant_code,how='inner')
Grant_code_Target=ID_CORE_Total_Target.merge(Grant_code,how='inner')
Grant_code=Grant_code_Target.merge(Grant_code_drug,how='outer')
Grant_code=Grant_code.merge(Search_Year,how='outer')
Grant_code=Grant_code.drop_duplicates()

Grant_code['PMID'].nunique()


# In[38]:


Grant_code = Grant_code[Grant_code['Search_ID'] != "237"]


# In[39]:


Grant_code.drop('Search__ID', axis=1, inplace=True)
Grant_code.drop('Source_Search_Type', axis=1, inplace=True)
Grant_code.drop('Search_Type', axis=1, inplace=True)


# In[40]:


################ UQ_COMBO (PMID|Brand) Maker
BRAND_ID['Search_ID'] = BRAND_ID['Search_ID'].apply(remove_decimal_from_id)
Grant_code=Grant_code.merge(BRAND_ID, how='outer')
Grant_code['PMID'] = Grant_code['PMID'].apply(remove_decimal_from_id)
#Grant_code['UQ_COMBO']=Grant_code['PMID'].astype(str)+"|"+Grant_code['Brand_Name']

Grant_code['UQ_COMBO_APY']=Grant_code['ACTUAL_PROJECT_YEAR'].astype(str)+"|"+Grant_code['Brand_Name']


# In[41]:


################ UQ_COMBO (PMID|Brand) Maker
BRAND_ID['Search_ID'] = BRAND_ID['Search_ID'].apply(remove_decimal_from_id)
Grant_code=Grant_code.merge(BRAND_ID, how='outer')
Grant_code['PMID'] = Grant_code['PMID'].apply(remove_decimal_from_id)
#Grant_code['UQ_COMBO']=Grant_code['PMID'].astype(str)+"|"+Grant_code['Brand_Name']

Grant_code['UQ_COMBO_PMID']=Grant_code['PMID'].astype(str)+"|"+Grant_code['Brand_Name']


# In[42]:


def COMBO_Linker(search_id):
    if 'drug' in search_id:
        return 'Drug'
    elif search_id.isdigit():
        return 'TargetOnly'
    else:
        return 'Error'

Grant_code['Data_Type_by_Drug'] = Grant_code['Search_ID'].apply(COMBO_Linker)

# Identify UQ_COMBO values which have at least one Drug
uq_combo_with_drug = Grant_code[Grant_code['Data_Type_by_Drug'] == 'Drug']['UQ_COMBO_APY'].unique()

# Update rows with those UQ_COMBO values to "Drug"
Grant_code.loc[Grant_code['UQ_COMBO_APY'].isin(uq_combo_with_drug), 'Data_Type_by_Drug'] = 'Drug'


# In[43]:


# Identify UQ_COMBO APY values which have at least one Drug
Grant_code['Data_Type_by_APY'] = Grant_code['Search_ID'].apply(COMBO_Linker)
uq_combo_with_drug = Grant_code[Grant_code['Data_Type_by_APY'] == 'Drug']['ACTUAL_PROJECT_YEAR'].unique()

# Update rows with those UQ_COMBO values to "Drug"
Grant_code.loc[Grant_code['ACTUAL_PROJECT_YEAR'].isin(uq_combo_with_drug), 'Data_Type_by_APY'] = 'Drug'


# In[44]:


# Identify UQ_COMBO PMID values which have at least one Drug
Grant_code['Data_Type_by_PMID'] = Grant_code['Search_ID'].apply(COMBO_Linker)
uq_combo_with_drug = Grant_code[Grant_code['Data_Type_by_PMID'] == 'Drug']['UQ_COMBO_PMID'].unique()

# Update rows with those UQ_COMBO values to "Drug"
Grant_code.loc[Grant_code['UQ_COMBO_PMID'].isin(uq_combo_with_drug), 'Data_Type_by_PMID'] = 'Drug'


# In[45]:


Grant_code = Grant_code.dropna(subset=['PMID'])
Grant_code = Grant_code[Grant_code['PMID'].ne('') & Grant_code['PMID'].ne(0)]


# In[46]:


####### TIME CUT +1 Year Correction (adjust for download year offset back to True approval cut year)
Grant_code.rename(columns={'TIME_CUT': 'TIME_CUT(Download_Used)'}, inplace=True)
Grant_code['TIME_CUT']=Grant_code['TIME_CUT(Download_Used)'].astype(int) - 1

Grant_code['TIME_CUT']=Grant_code['TIME_CUT'].apply(remove_decimal_from_id)

Grant_code['TIME_CUT']=Grant_code['TIME_CUT'].astype(str)


# In[47]:


#### 2018USD inflation adjustment

Grant_code.rename(columns = {"APY": "FY"}, 
          inplace = True)

Grant_code.rename(columns = {"APY_COST_inf2018": "Original_COST"}, 
          inplace = True)

Inflation_key_V2018 = pd.read_csv('function_data/inf2018_key.csv')
core_set_2018=pd.merge(Grant_code,Inflation_key_V2018, how='outer')
core_set_2018['APY_COST_inf2018']=core_set_2018['Original_COST']*core_set_2018['inf_2018']
core_set_2018.drop('inf_2018', inplace=True, axis=1)

Grant_code=core_set_2018

Grant_code.rename(columns = {"FY": "APY"}, 
          inplace = True)


# In[48]:


Grant_code = Grant_code.dropna(subset=['PMID'])
Grant_code = Grant_code[Grant_code['PMID'].ne('') & Grant_code['PMID'].ne(0)]


# In[51]:


######## Time cut by PMID PubYear to approval year
Grant_code['TIME_CUT_FIC'] = Grant_code['TIME_CUT']

#### First in class for targets adjustment ( Use for Basic/Target search ID that have an earlier FIC drug approval as cut year )
# Grant_code.loc[Grant_code['Search_ID'] == "135", 'TIME_CUT_FIC'] = "2013"
# Grant_code.loc[Grant_code['Search_ID'] == "63", 'TIME_CUT_FIC'] = "2011"

Grant_code['Years_from_app'] = Grant_code['APY'].astype(int) - Grant_code['TIME_CUT'].astype(int)
Grant_code['Years_from_app_FIC_TARGET'] = Grant_code['APY'].astype(int) - Grant_code['TIME_CUT_FIC'].astype(int)

Grant_code['FY_START_Years_from_app_FIC_TARGET'] = Grant_code['FY_Start'].astype(int) - Grant_code['TIME_CUT_FIC'].astype(int)


Grant_code['Pre_Approval_APY'] = Grant_code['Years_from_app'] <= 0
Grant_code['Pre_Approval_APY_FIC'] = Grant_code['Years_from_app_FIC_TARGET'] <= 0


# In[52]:


#### Use for data from 4+ Years after approval (not part of main method)
# Grant_code['4YEAR_RULE_1'] = Grant_code['Years_from_app_FIC_TARGET'].apply(lambda x: 1 if 4 >= x >= -12 else 0)
# Grant_code['4YEAR_RULE_2'] = Grant_code['FY_START_Years_from_app_FIC_TARGET'].apply(lambda x: 1 if 0 >= x else 0)

# Grant_code['APY_Start_PreApproval_4Year'] = np.where(Grant_code['4YEAR_RULE_1'] + Grant_code['4YEAR_RULE_2'] == 2, 1, 0)

# Grant_code.drop('4YEAR_RULE_1', axis=1, inplace=True)
# Grant_code.drop('4YEAR_RULE_2', axis=1, inplace=True)

####  Last year of total costs change >2% per year for 10 Drug dataset, 9_12_23
Grant_code['0_12_Year_Rule_FIC'] = Grant_code['Years_from_app_FIC_TARGET'].apply(lambda x: 1 if 0 >= x >= -12 else 0)
Grant_code=Grant_code.fillna(0)
Grant_code_output=Grant_code


# In[53]:


cols_order = ["Brand_Name", 
"Search_ID", 
"PMID", 
"PUB_YEAR", 
"PROJECT_NUMBER", 
"FY_Start", 
"FY_Last", 
"APY", 
"ACTUAL_PROJECT_YEAR", 
"APY_COST_inf2018", 
"Original_COST", 
"Activity_Code", 
"Institute_Code", 
"Acronym_institute_name", 
"full_institute_name", 
"Compressed Names", 
"Project_Count", 
"Grant_Type_Name", 
"Data_Type_by_Drug", 
"Data_Type_by_APY", 
"Data_Type_by_PMID", 
"UQ_COMBO_APY", 
"UQ_COMBO_PMID", 
"TIME_CUT(Download_Used)", 
"TIME_CUT", 
"TIME_CUT_FIC", 
"Years_from_app", 
"Years_from_app_FIC_TARGET", 
"FY_START_Years_from_app_FIC_TARGET", 
"Pre_Approval_APY", 
"Pre_Approval_APY_FIC", 
"0_12_Year_Rule_FIC"
]
Grant_code_output = Grant_code_output[cols_order]


# In[ ]:





# In[54]:


Grant_code_output.to_csv("IRA_10Drugs_Main_Output.csv",index=None)


# In[55]:


Grant_code['CC_ID']= Grant_code['Brand_Name']+"_"+Grant_code['Data_Type_by_Drug']+"_"+Grant_code['ACTUAL_PROJECT_YEAR']


# In[56]:


Grant_code_ID=Grant_code[['CC_ID','Data_Type_by_Drug','Brand_Name']]
Grant_code_ID=Grant_code_ID.drop_duplicates()


# In[57]:


def Applied_Basic_research_split(row):
    if row['Data_Type_by_Drug'] == "Drug":
        if isinstance(row['Search_ID'], (int, float)) and not pd.isna(row['Search_ID']):
            return False  
    return True

Grant_code = Grant_code[Grant_code.apply(Applied_Basic_research_split, axis=1)]


# In[58]:


############ USE FOR FIC CC analysis only
Grant_code_CC=Grant_code.copy()

Grant_code_CC=Grant_code_CC[['CC_ID','APY','Original_COST','TIME_CUT_FIC','Data_Type_by_Drug']]
Grant_code_CC=Grant_code_CC.drop_duplicates()


Grant_code_CC['Years_from_app'] = Grant_code_CC['APY'].astype(int) - Grant_code_CC['TIME_CUT_FIC'].astype(int)
Grant_code_CC['Pre_Approval_APY'] = Grant_code_CC['Years_from_app'] <= 0
Grant_code_CC = Grant_code_CC[Grant_code_CC['Pre_Approval_APY'] == 1]
Grant_code_CC = Grant_code_CC[Grant_code_CC['Years_from_app'] >= -12]

Grant_code_CC=Grant_code_CC[['CC_ID','Years_from_app','Original_COST']]
Grant_code_CC=Grant_code_CC.drop_duplicates()


# In[ ]:





# In[59]:


Grant_code_CC=pd.pivot_table(Grant_code_CC,index='CC_ID',columns='Years_from_app',values='Original_COST', aggfunc='sum')
Grant_code_CC=Grant_code_CC.fillna(0)


# In[60]:


Grant_code_CC.to_csv('CC_CHECK.csv')


# In[61]:


Grant_code_CC=pd.read_csv('CC_CHECK.csv')


# In[62]:


Grant_code_CC=Grant_code_CC.merge(Grant_code_ID,how='inner')


# In[63]:


Grant_code_CC.to_csv('CC_DEBUG.csv',index=None)


# In[64]:


Grant_code_CC_Grouped = Grant_code_CC.groupby(['Brand_Name', 'Data_Type_by_Drug']).sum().reset_index()


# In[65]:


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


year_cols = [col for col in Grant_code_CC_Grouped.columns if is_number(col)]
sorted_cols = sorted(year_cols, key=float, reverse=True)  # Sorts 0, -1, -2...

for col in sorted_cols:
    Grant_code_CC_Grouped[col] = Grant_code_CC_Grouped[col].astype(int)
Grant_code_CC_Grouped['Discount_Type']= 'Original'
final_cols = ['Brand_Name', 'Data_Type_by_Drug','Discount_Type'] + sorted_cols
Grant_code_CC_Grouped = Grant_code_CC_Grouped[final_cols]

#Grant_code_CC_Grouped.to_csv('CC_10Drugs_Raw.csv',index=None)
Grant_code_CC_Grouped_Raw=Grant_code_CC_Grouped.copy()


# In[ ]:


def modify_rightmost_column(df):
    # Identify the rightmost column
    rightmost_col = df.columns[-1]
    
    # Multiply all values in the rightmost column by 1 (Hand analysis for 10.5% or 3% as needed)
    df[rightmost_col] *= 1
    
    return df

df_modified = modify_rightmost_column(Grant_code_CC_Grouped)


# In[67]:


import pandas as pd

def modify_adjacent_column(df, col_name):
    # Identify the specified column and the one to its left
    adjacent_col = df.columns[df.columns.get_loc(col_name) - 1]
    
    # Check if both columns are numeric
    if pd.api.types.is_numeric_dtype(df[col_name]) and pd.api.types.is_numeric_dtype(df[adjacent_col]):
        # For each row, add the value from the specified column to the value in the adjacent column, then multiply by 1.105
        df[adjacent_col] = (df[adjacent_col] + df[col_name]) * 1
    
    return df

def modify_until_discount(df):
    # Begin with the column left of the rightmost column and move leftward until "Discount_Type"
    col_index = len(df.columns) - 2  # Starting one column left of the rightmost column
    while col_index > 0:  # Safety to ensure we don't go out of bounds
        if df.columns[col_index] == "Discount_Type":
            break
        df = modify_adjacent_column(df, df.columns[col_index])
        col_index -= 1

    return df

df_modified = modify_until_discount(df_modified)


# In[68]:


df_modified['Discount_Type'] ="0%"


# In[69]:


df_modified=pd.concat([df_modified, Grant_code_CC_Grouped_Raw])


# In[70]:


df_modified.sort_values(by='Brand_Name')
df_modified.sort_values(by=['Discount_Type', 'Brand_Name'])
df_modified.to_csv('CC_10Drugs_0%_HandCheck.csv',index=None)


# In[72]:


Grant_code_FIC_RULE = Grant_code.loc[Grant_code['0_12_Year_Rule_FIC'] == 1]


# In[73]:


Grant_code_unique_COST = Grant_code_FIC_RULE.drop_duplicates(subset=['Brand_Name', 'Data_Type_by_Drug', 'ACTUAL_PROJECT_YEAR', 'APY_COST_inf2018'])

result_COST = Grant_code_unique_COST.groupby(['Brand_Name', 'Data_Type_by_Drug'])['APY_COST_inf2018'].sum().reset_index()

result_COST = result_COST.rename(columns={'APY_COST_inf2018':'COST_SUM'})


# In[74]:


Grant_code_unique_APY=Grant_code_FIC_RULE[['Brand_Name', 'Data_Type_by_Drug', 'ACTUAL_PROJECT_YEAR']]
Grant_code_unique_APY=Grant_code_unique_APY.drop_duplicates()

result_APY = Grant_code_unique_APY.groupby(['Brand_Name', 'Data_Type_by_Drug'])['ACTUAL_PROJECT_YEAR'].nunique().reset_index()
result_APY.rename(columns={'ACTUAL_PROJECT_YEAR': 'Count_of_APY(0-12Year)'}, inplace=True)


# In[75]:


Grant_code_unique_PMID=Grant_code_FIC_RULE[['Brand_Name', 'Data_Type_by_PMID', 'PMID']]
Grant_code_unique_PMID=Grant_code_unique_PMID.drop_duplicates()


Grant_code_unique_PMID.rename(columns={'Data_Type_by_PMID': 'Data_Type_by_Drug'}, inplace=True)
result_PMID = Grant_code_unique_PMID.groupby(['Brand_Name', 'Data_Type_by_Drug'])['PMID'].nunique().reset_index()
result_PMID.rename(columns={'PMID': 'NIH_PMID(0_12Year)'}, inplace=True)


result_Items = result_PMID.merge(result_APY, on=['Brand_Name', 'Data_Type_by_Drug'], how='left')
Item_Cost_Overview_12_YEARS = result_Items.merge(result_COST, on=['Brand_Name', 'Data_Type_by_Drug'], how='left')

Item_Cost_Overview_12_YEARS['COST_SUM_Mills']=Item_Cost_Overview_12_YEARS['COST_SUM']/1000000
Item_Cost_Overview_12_YEARS.drop('COST_SUM', axis=1, inplace=True)


# In[76]:


Grant_code_unique_COST = Grant_code.drop_duplicates(subset=['Brand_Name', 'Data_Type_by_Drug', 'ACTUAL_PROJECT_YEAR', 'APY_COST_inf2018'])

result_COST = Grant_code_unique_COST.groupby(['Brand_Name', 'Data_Type_by_Drug'])['APY_COST_inf2018'].sum().reset_index()

result_COST = result_COST.rename(columns={'APY_COST_inf2018':'COST_SUM'})

Grant_code_unique_APY=Grant_code[['Brand_Name', 'Data_Type_by_Drug', 'ACTUAL_PROJECT_YEAR']]
Grant_code_unique_APY=Grant_code_unique_APY.drop_duplicates()

#Grant_code_unique_APY.rename(columns={'Data_Type_by_APY': 'Data_Type_by_Drug'}, inplace=True)

result_APY = Grant_code_unique_APY.groupby(['Brand_Name', 'Data_Type_by_Drug'])['ACTUAL_PROJECT_YEAR'].nunique().reset_index()
result_APY.rename(columns={'ACTUAL_PROJECT_YEAR': 'Count_of_APY(ALL_YEAR)'}, inplace=True)

Grant_code_unique_PMID=Grant_code[['Brand_Name', 'Data_Type_by_PMID', 'PMID']]
Grant_code_unique_PMID=Grant_code_unique_PMID.drop_duplicates()


Grant_code_unique_PMID.rename(columns={'Data_Type_by_PMID': 'Data_Type_by_Drug'}, inplace=True)
result_PMID = Grant_code_unique_PMID.groupby(['Brand_Name', 'Data_Type_by_Drug'])['PMID'].nunique().reset_index()
result_PMID.rename(columns={'PMID': 'NIH_PMID(ALL_YEAR)'}, inplace=True)


result_Items = result_PMID.merge(result_APY, on=['Brand_Name', 'Data_Type_by_Drug'], how='left')
Item_Cost_Overview_ALL_YEARS = result_Items.merge(result_COST, on=['Brand_Name', 'Data_Type_by_Drug'], how='left')

Item_Cost_Overview_ALL_YEARS['COST_SUM_Mills(ALL_YEAR)']=Item_Cost_Overview_ALL_YEARS['COST_SUM']/1000000
Item_Cost_Overview_ALL_YEARS.drop('COST_SUM', axis=1, inplace=True)


# In[77]:


Item_Cost_Overview_ALL = Item_Cost_Overview_12_YEARS.merge(Item_Cost_Overview_ALL_YEARS, on=['Brand_Name', 'Data_Type_by_Drug'], how='inner')
Item_Cost_Overview_ALL = PMID_DOWNLOADED_Analysis.merge(Item_Cost_Overview_ALL, on=['Brand_Name', 'Data_Type_by_Drug'], how='inner')


# In[79]:


Item_Cost_Overview_ALL['%_PMID_Funded(0_12)']=Item_Cost_Overview_ALL['NIH_PMID(0_12Year)']/Item_Cost_Overview_ALL['PMID_Downloaded']
Item_Cost_Overview_ALL['%_PMID_Funded(ALL_YEAR)']=Item_Cost_Overview_ALL['NIH_PMID(ALL_YEAR)']/Item_Cost_Overview_ALL['PMID_Downloaded']


# In[81]:


Item_Cost_Overview_ALL.to_csv('Item_Cost_Overview.csv',index=None)
os.remove('CC_CHECK.csv')

