# NIH-funding-for-first-round-inflation-reduction-act-drugs-approvals
Federal financial support for discovery and development of drugs selected for price negotiation in year one of the Inflation Reduction Act: Python analysis code
See README file in function_data for RePORTER input reference link table creation if needed.  

Last Update 11-09-2023| by Edward W. Zhou


###########################################################################################################

This code takes search terms for basic (target / Biological target or Mechanism Of Action or MOA from: SearchTerm_TargetOnly) and applied (drug names from :SearchTerm_Drug) to get total pubmed PMID list
and NIH funding from RePORTER (Last Exporter csv download 9/18/23 - Accepted data range 1985 - 2020).

Step_1_PMID_APY_DOWNLOAD__10_31_23: Takes search_term_inputs csv and output PMID list on terms and time range (TIME_CUT)
Step_2_FUNDING_Analysis_10_31_23  : Format output to per drug funding totals (see paper for method section)

	Basic: Any Search_ID with only numbers

	Applied: Any Search_ID with drug###

contact ezhou@bentley.edu for any questions/assistance on code.  
See README file in function_data for RePORTER input file creation if needed.  



#NOTE: 
    #Folders must be on the same folder layer as code

    
    # Make sure TIME_CUT_TRUE = real approval year / TIME_CUT = Real year +1 (if using same method as from paper)
    # Code downloads one year at a time to bypass 10k per request PubMed Download limit (input term can fail if over 10,000 paper in a single year - no error will be shown for this so be careful)
    
    #function_data = RePORTER local link for grant data on PMID 1985-2020 (some inconsistency exists between link csv and NIH RePORTER website outputs, likely due to inconsistent upkeep on csv archives)

    #Make sure search_term_inputs\BrandName_Linker + SearchTerm_Drug + SearchTerm_TargetOnly line up if new inputs are needed
    #search_term_inputs >> BrandName_Linker.csv = your SearchID and Brand_Name cols (KEEP AS 2 COL, So Brand_Name col may not always be unique)
    #search_term_inputs >> SearchTerm_Drug and SearchTerm_TargetOnly.csv = Input / SearchID / TIME_CUT
        #Input = PubMed Search term of interest (DO NOT ADD "NOT REVIEW[pt] or TIME RANGEs", code adds them later if needed)
        #SearchID = This should be the Drug## or target number of the inputs (or whatever you want them to be, DO NOT LEAVE BLANK)
            # for step 2 to work, Drug ID must start with 'drug' , and target ID must be all numbers (234 for example)
            # if step 2 is not needed, any unique ID will work


#NOTE: csv / os / time / pandas / numpy should be a part of most python interface tools pre-load list
#pip install biopython
#pip install beautifulsoup4
#pip install Bio
