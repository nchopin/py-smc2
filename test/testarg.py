from snippets.localfolder import get_path
import os, imp, sys 

nargs = len(sys.argv)
if nargs == 1:
    print "using default userfile ('userfile.py' in the root folder)..."
    userfilepath = "userfile.py"
if nargs == 2:
    userfilepath = sys.argv[1]
    print "working with user file %s" % userfilepath 
if nargs > 2:
    raise ValueError, "stop: too many arguments; you should only specify a path to a user file"

THISPATH = get_path()
userfilefolder =  os.path.join(THISPATH, os.path.dirname(userfilepath))
userfilebasename = os.path.basename(userfilepath)
userfilebasename = userfilebasename.replace(".py", "")
sys.path.append(userfilefolder)
f, filename, description = imp.find_module(userfilebasename)
userfile = imp.load_module("userfile", f, filename, description)
