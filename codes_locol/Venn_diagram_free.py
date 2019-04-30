#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 00:00:34 2017
@author: yolandatiao

Updated on Tue Apr 30 15:54:00 2019
@author: yolandatiao
"""

#####------------------ Import START ------------------#####
import os # For changing directory
import subprocess # For calling bash 
from astropy.io import ascii # For using ascii table to open csv
from astropy.table import Table, Column, join   # For using astropy table functions
import matplotlib.pyplot as plt # For plotting
from matplotlib_venn import *
import Tkinter,tkFileDialog
import matplotlib.pyplot as plt
#####------------------ Import END ------------------#####



#####------------------ Config START ------------------#####

#####------------------ Config END ------------------#####



#####------------------ Self Defined functions START ------------------#####
def strlist(listx):
    outlist=[]
    for i in listx:
        outlist.append(str(i))
    return outlist

# A function to get filename without file formats
def filenamenoformat(file_name_full):
	file_name_full_list=file_name_full.split('.')
	type_len = len(file_name_full_list[-1])
	full_len = len(file_name_full)
	real_len = full_len - type_len - 1
	return file_name_full[0:real_len]

# A function to strip the path of the file
def Getfilename(name_x):
	return name_x.split('/')[-1]

# A function to get the path of the file
def Getfilepath(name_x):
	file_path_list=name_x.split("/")[:-1]
	file_path=("/".join(file_path_list))
	return file_path

# A function to set colnames if colnames can not be read
def setcolnames(table_x):
	old_colnames=table_x.colnames
	if old_colnames[0] == 'col1':
		new_colnames=list(table_x[0])
		table_x.remove_row(0)
		for x in xrange(0, len(table_x[0])):
			table_x[old_colnames[x]].name=new_colnames[x]
	return table_x
#####------------------ Self Defined functions END ------------------######


#####------------------ Main START ------------------######
#--- Open input files
root = Tkinter.Tk()
root.withdraw()
root.update()
print 'Please choose files for making venn diagram: '
in_filenames = list(tkFileDialog.askopenfilenames(parent=root,title='Choose a file'))
in_filenames_len = len(in_filenames) #Get number of opened files
#--- Check numbers and names of files
print "You are opening %d files: " % (in_filenames_len) 
for x in xrange (0, in_filenames_len):
	print in_filenames[x].split('/')[-1] #Get rid of directory
print '\n'

#--- Name files
print 'Please name the data from the files, seperate by comma: '
for x in xrange(0, in_filenames_len):
    print filenamenoformat(Getfilename(in_filenames[x]))
plot_names=raw_input()
plot_nameslist=plot_names.split(",")
print '\n'

#--- Output name
o_n_list=[]
for i in in_filenames:
    o_n_list.append(filenamenoformat(Getfilename(i)))
outname="--".join(o_n_list)
outname="%s.pdf"%outname
    

#--- change dir
wkdir=Getfilepath(in_filenames[0])
os.chdir(wkdir)

#--- Getting data lists from files
if len(in_filenames) == 3:
    data_1=ascii.read(in_filenames[0])
    data_2=ascii.read(in_filenames[1])
    data_3=ascii.read(in_filenames[2])
    
    data_1=setcolnames(data_1)
    data_2=setcolnames(data_2)
    data_3=setcolnames(data_3)
    
    set_1=set(list(data_1.columns[0]))
    set_2=set(list(data_2.columns[0]))
    set_3=set(list(data_3.columns[0]))
    
    setlist=[set_1, set_2, set_3]
    
    for x in xrange(0, len(in_filenames)):
        print "%s %s" %(x, plot_nameslist[x])
    print "1st circle:"
    r_n=int(raw_input())
    print '\n'
    print "2nd circle:"
    g_n=int(raw_input())
    print '\n'
    print "3rd circle:"
    p_n=int(raw_input())
    print '\n'    
    keylist=[r_n,g_n,p_n]
    
    for x in xrange(0, len(in_filenames)):
        print "%s %s" %(x, plot_nameslist[keylist[x]])
    print "Solid outline circle: (non:99)"
    sc=int(raw_input())
      
    v3=venn3([setlist[r_n],setlist[g_n], setlist[p_n]],(plot_nameslist[r_n],plot_nameslist[g_n],plot_nameslist[p_n]))
    if sc == 99:
        plt.savefig(outname)
    else:
        idlist=[]
        for x in xrange(0,3):
            if x == sc:
                idlist.append("1")
            else:
                idlist.append("0")
        idxyz=str("".join(strlist(idlist)))
        v3.get_patch_by_id(idxyz).set_alpha(0)
        v3c=venn3_circles([setlist[r_n],setlist[g_n], setlist[p_n]],linestyle='solid')
        for x in xrange(0,len(in_filenames)):
            if x != sc:
                v3c[x].set_lw(0)
        plt.savefig(outname)
    
    
elif len(in_filenames) == 2:
    data_1=ascii.read(in_filenames[0])
    data_2=ascii.read(in_filenames[1])
    
    data_1=setcolnames(data_1)
    data_2=setcolnames(data_2)
    
    set_1=set(list(data_1.columns[0]))
    set_2=set(list(data_2.columns[0]))
    
    setlist=[set_1, set_2]
    
    for x in xrange(0, len(in_filenames)):
        print "%s %s" %(x, plot_nameslist[x])
    print "1st circle:"
    r_n=int(raw_input())
    print '\n'
    print "2nd circle:"
    g_n=int(raw_input())
    print '\n'
    keylist=[r_n,g_n]
    
    for x in xrange(0, len(in_filenames)):
        print "%s %s" %(x, plot_nameslist[keylist[x]])
    print "Solid outline circle: (non:99)"
    sc=int(raw_input())
    
    v2=venn2([setlist[r_n],setlist[g_n]],(plot_nameslist[r_n],plot_nameslist[g_n]))
    if sc == 99:
        plt.savefig(outname)
    else:
        idlist=[]
        for x in xrange(0,2):
            if x == sc:
                idlist.append("1")
            else:
                idlist.append("0")
        idxy=str("".join(strlist(idlist)))
        v2.get_patch_by_id(idxy).set_alpha(0)
        v2c=venn2_circles([setlist[r_n],setlist[g_n]],linestyle='solid')
        for x in xrange(0,len(in_filenames)):
            if x != sc:
                v2c[x].set_lw(0)
        plt.savefig(outname)
#####------------------ Main END ------------------######































