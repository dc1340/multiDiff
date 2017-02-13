# /usr/bin/python

from selenium import webdriver
from string import ascii_lowercase
import os
from time import sleep
import sys
from selenium import webdriver
from selenium.webdriver.common.keys import Keys


browser = webdriver.Chrome()
time_delay=1

input_file=open('/Users/dhruvachandramohan/GitHub/multiDiff/p2_output/cd_files',
 'r')
file_list=input_file.read().split('\n')
input_file.close()

base_url='http://amp.pharm.mssm.edu/Enrichr/'

#Upload loop
for f in file_list:
	browser.find_element_by_tag_name('body').send_keys(Keys.COMMAND + 't') 
	file_lab=os.path.basename(f).split('.')[0]
	print ("Processing "+file_lab)
	gene_list=open(f, 'r').read()


	old_file='/Users/dhruvachandramohan/Downloads/GO_Biological_Process_2015_table.txt'
	new_file='/Users/dhruvachandramohan/GitHub/multiDiff/p2_output/enricher_pathways/'+file_lab+'.txt'

	if (len(gene_list)>0 & !os.path.isfile(new_file) ):
		browser.get(base_url)
		
		##Enter gene list
		textElem =browser.find_element_by_id('text-area')
		textElem.send_keys(gene_list)
		
		##Enter description
		textElem =browser.find_element_by_name('description')
		textElem.send_keys(file_lab)

		buttonElem =browser.find_element_by_id('proceed-button')
		buttonElem.click()
		#sleep(time_delay)

		

		buttonElem =browser.find_element_by_id('Ontologies-link')
		buttonElem.click()
		#sleep(time_delay)

		buttonElem =browser.find_element_by_id('GO_Biological_Process_2015-link')
		buttonElem.click()
		#sleep(time_delay)

		buttonElem =browser.find_element_by_id('GO_Biological_Process_2015-Table-link')
		webdriver.ActionChains(browser).move_to_element_with_offset(buttonElem, 2, 2).click().perform()
		#sleep(time_delay)

		#tableElem =browser.find_element_by_id('GO_Biological_Process_2015-TableExport-link')
		#tableElem.click()
		#sleep(time_delay)

		#old_file='/Users/dhruvachandramohan/Downloads/GO_Biological_Process_2015_table.txt'
		#new_file='/Users/dhruvachandramohan/GitHub/multiDiff/p2_output/enricher_pathways/'+file_lab+'.txt'

		# while True:
		#         if os.path.isfile(old_file+'.part'):
		#             sleep(time_delay)
		#         elif os.path.isfile(old_file+'.part'):
		#             break
		#         else:
		#             sleep(time_delay)

		#os.rename( old_file, new_file)
		browser.find_element_by_tag_name('body').send_keys(Keys.COMMAND + 't') 
		print ("...Done")


