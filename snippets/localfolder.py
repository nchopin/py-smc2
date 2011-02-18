#!/usr/local/bin/python
import pat  # touch pat.py in the cwd and try to import the empty file under Linux
import string 
import os,sys

# This function is only used to get the path to
# the folder where these files are, it has
# nothing to do with block IMH.

################## os.getcwd()##############
def get_path():
    PAT = str(pat).replace("<module 'snippets.pat' from ", "")[1:-9]
    PAT = PAT.replace("/snippets", "")
    #PAT=str(pat).split()[3][1:-9] # PATH extracted..
    sig=None
    try:
        sig=os.remove(os.path.join(PAT,'pat.pyc'))# get_rid...
    except OSError:
        pass
#        PAT=PAT +'/'
#        sig=os.remove(PAT + 'pat.pyc')# Fix for mutiple calls..  
    return PAT
###############################
#LOCATE=get_path()
#print LOCATE # 

