# /usr/bin/python

from selenium import webdriver
from string import ascii_lowercase
import os
from time import sleep
import sys
import pandas

browser = webdriver.Chrome()
time_delay=4
download_dir='/Users/dhruvachandramohan/Downloads/'
#outdir='/Users/dhruvachandramohan/GitHub/multiDiff/p2_output/selected_enricher/'

if len(sys.argv)<4:
	outdir='/Users/dhruvachandramohan/GitHub/multiDiff/p2_output/selected_enricher_3_dmc_prom/'
else:
	outdir=sys.argv[3]
# run_ont-True
# run_pathways=True

## Read in input file, each line of which is a path to file
## with gene names
# print(sys.argv,'\n')
input_file=open(sys.argv[2], 'r')
file_list=input_file.read().split('\n')
input_file.close()


base_url='http://amp.pharm.mssm.edu/Enrichr/'

input_file=open('/Users/dhruvachandramohan/GitHub/multiDiff/annotations/enricher/category_list.txt', 'r')
category_list=input_file.read().split('\n')
input_file.close()

#analysis_df=pandas.read_table('/Users/dhruvachandramohan/GitHub/multiDiff/annotations/enricher/enricher_annotation_list.txt')
analysis_df=pandas.read_table(sys.argv[1])

#pathways_list=[ pathways_list[ i ] for i in select_pways ]

#pathways_list=pathways_list[ range(0,2 )]

def download_table(x):
	buttonElem =browser.find_element_by_id(x+'-link')
	webdriver.ActionChains(browser).move_to_element_with_offset(buttonElem, 0, 0).click().perform()
	
	sleep(2)
	
	buttonElem =browser.find_element_by_id(x+'-Table-link')
	webdriver.ActionChains(browser).move_to_element_with_offset(buttonElem, 5, 3).click().perform()
	
	sleep(2)

	tableElem =browser.find_element_by_id(x+'-TableExport-link')
	webdriver.ActionChains(browser).move_to_element_with_offset(tableElem , 0 , 0).click().perform()
	

#Upload loop
def download_selected_enricher_from_abs_to_genelist(f):
	
	file_lab=os.path.basename(f).split('.')[0]
	print ("Processing "+file_lab)
	gene_list=open(f, 'r').read()

	if len(gene_list)>0:
	
		browser.get(base_url)
		
		##Enter gene list
		textElem =browser.find_element_by_id('text-area')
		textElem.send_keys(gene_list)
		
		##Enter description
		textElem =browser.find_element_by_name('description')
		textElem.send_keys(file_lab)

		buttonElem =browser.find_element_by_id('proceed-button')
		buttonElem.click()
		sleep(time_delay)

		#Process each category
		for category in category_list:
			print ("..."+category)
			
			analysis_list=analysis_df.loc[analysis_df['Category'] ==category][ 'Annotation' ]
			#print(analysis_list.shape)
			if (analysis_list.shape[0]==0):
				continue	

			if category=="Disease/Drugs":
				cat_string="Disease_or_Drugs"
				fixed_cat=category
			elif category=="Cell Types":
				fixed_cat='Cell_Types'
				cat_string=category
			else:
				cat_string=category
				fixed_cat=category

			buttonElem =browser.find_element_by_id(fixed_cat+'-link')
			buttonElem.click()
			sleep(time_delay)

			#analysis_list=analysis_df.loc[analysis_df['Category'] ==category][ 'Annotation' ]
			#Process each analysis
			for analysis in analysis_list:
				print ("......"+analysis)
				
				old_file=download_dir+analysis +'_table.txt'
				new_file='%s/%s/%s_%s.txt' % (outdir, cat_string,file_lab,analysis)
				
				if not os.path.isfile(new_file):
					
					download_table(analysis)
					
					while not os.path.isfile(old_file):
					    sleep(1)
					os.rename( old_file, new_file)
	
		print ("...Done")

for i in file_list:
	if i:
		download_selected_enricher_from_abs_to_genelist(i)

browser.quit()