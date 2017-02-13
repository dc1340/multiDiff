# /usr/bin/python

from selenium import webdriver
from string import ascii_lowercase
import os
from time import sleep

browser = webdriver.Chrome()
time_delay=10	


base_url='http://amp.pharm.mssm.edu/Enrichr/'

#Upload loop
for f in file_list:
	file_lab=os.path.basename(f).split('.')[0]
	gene_list=[line.strip() for line in open(f, 'r')]

	browser.get(base_url)

	textElem =browser.find_element_by_id('black')
	textElem.
	
	print ("Processing "+file_lab)

	buttonElem =browser.find_element_by_id('Ontologies-link')
	buttonElem.click()
	sleep(time_delay)

	buttonElem =browser.find_element_by_id('GO_Biological_Process_2015-link')
	buttonElem.click()
	sleep(time_delay)

	buttonElem =browser.find_element_by_id('GO_Biological_Process_2015-Table-link')
	buttonElem.click()
	sleep(time_delay)

	tableElem =browser.find_element_by_id('GO_Biological_Process_2015-TableExport-link')
	tableElem.click()
	sleep(time_delay)

	old_file='/Users/dhruvachandramohan/Downloads/GO_Biological_Process_2015_table.txt'
	new_file='/Users/dhruvachandramohan/Dropbox/P2_ERRBS/enricher_pathways/'+file_lab+'.txt'
	os.rename( old_file, new_file)


