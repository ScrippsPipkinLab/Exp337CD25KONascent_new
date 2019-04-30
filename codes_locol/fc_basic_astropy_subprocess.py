#####----- Updated 2017-04-06 -----#####

#import tkinter    # For opening file user interface
#from tkFileDialog import askopenfilename    # For opening file user interface
from astropy.io import ascii    # For using ascii table to open csv
from astropy.table import Table, Column    # For using astropy table functions
import string    # For creating 96 well plate list
import xlsxwriter    # For writing in Excel
import re    # For using the char filter
import numpy as np  # For use Astropy table
import seaborn as sns # For plotting
import matplotlib.pyplot as plt # For plotting
import statistics # For statistics calculations
import xlsxwriter

# A function to call bash commands
def run_sh(command):                                     
  process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  retcode = process.wait()
  if retcode != 0:
    raise Exception, "Problem running command: " + command
  stdout, stderr = process.communicate()
  return stdout

# A function to strip the path of the file
def Getfilename(name_x):
	return name_x.split('/')[-1]

# A function to get the path of the file
def Getfilepath(name_x):
	file_path_list=name_x.split("/")[:-1]
	file_path=("/".join(file_path_list))
	return file_path

# A function to get astropy data type name from type function call output
def astropy_typename(str_x):
	str_x=str_x.replace("<","")
	str_x=str_x.replace(">","")
	str_x=str_x.replace("_","")
	str_x_list=str_x.split(" ")
	str_x=str_x_list[1]
	str_x=str_x.replace("'","")
	str_x_list=str_x.split(".")
	str_x=str_x_list[-1]
	return str_x

# A function to remove empty columns and empty rows in table
def RMempty(table_x):
	table_col_len=len(table_x[0])
	table_row_len=len(table_x)
	remove_column_list=[]
	remove_row_list=[]
	for x in xrange(0,table_col_len):    # Loop all columns
		type_x=astropy_typename(str(type(table_x[1][x])))
		if "Masked" in type_x:
			remove_column_list.append(table_x.colnames[x])
	table_x.remove_columns(remove_column_list)
	for y in xrange(0,table_row_len):    # Loop all rows
		type_y=astropy_typename(str(type(table_x[y][1])))
		if "Masked" in type_y:
			remove_row_list.append(y)
	table_x.remove_rows(remove_row_list)
	return table_x

# A function to set colnames if colnames can not be read
def setcolnames(table_x):
	old_colnames=table_x.colnames
	if old_colnames[0] == 'col1':
		new_colnames=list(table_x[0])
		table_x.remove_row(0)
		for x in xrange(0, len(table_x[0])):
			table_x[old_colnames[x]].name=new_colnames[x]
	return table_x

# A function to generate list from string seperated by comma
def strlist(str_x):
	str_list = str_x.split(',')
	return str_list

# A function to generate list from string seperated by |	
def strlist_1(str_x):
	str_list = str_x.split('|')
	return str_list

# A function to generate integer list from string seperated by comma
def intlist(str_x):
	str_list = str_x.split(',')
	return map(int, str_list)

# A function to generate float list from string seperated by comma
def floatlist(str_x):
	str_list = str_x.split(',')
	return map(float, str_list)

# A function to get filename without file formats
def filenamenoformat(file_name_full):
	file_name_full_list=file_name_full.split('.')
	type_len = len(file_name_full_list[-1])
	full_len = len(file_name_full)
	real_len = full_len - type_len - 1
	return file_name_full[0:real_len]

# A function to decide if an element is NaN or not
def isNaN(num):
	return num != num

# A function to write astropy table into xlsx worksheet
def table_write(tab, sheet):
	tab_col_no = len(tab[0])
	tab_row_no = len(tab)
	# Write colnames
	for loop_col in xrange(0, tab_col_no):
		sheet.write(0, loop_col, tab.colnames[loop_col])
	# Write data
	for loop_row in xrange(0, tab_row_no):
		for loop_col in xrange(0, tab_col_no):
			sheet.write(loop_row+1, loop_col, tab[loop_row][loop_col])







