import numpy
import sys
import os
import math
import scipy.misc
import ctypes
import revolve

def init(argv):
    argc=len(sys.argv)
    if(argc != 3):
        sys.exit("Usage: [timesteps] [memory checkpoints available]")
    files=sys.argv
    files.pop(0) 
    return files

args=init(sys.argv)

fine=int(args[0])  # total number of timesteps
snaps_in=int(args[1]) # number of memory checkpoints available

check=int(-1)
capo=int(0)
info=int(2)
ret=revolve.whatodo(check,capo,fine,snaps_in,info)
oldcapo=capo
whatodo=ret[0]
check=ret[1]
capo=ret[2]
fine=ret[3]
info=ret[4]
if(whatodo==1): 
    print "Advance from" + str(oldcapo) + " to " + str(capo) + "."
if(whatodo==2): 
    print "Store in checkpoint number " + str(check) + "."
if(whatodo==3): 
    print "First turn: Initialize adjoints and reverse first step."
if(whatodo==4): 
    print "Forward and reverse one step."
if(whatodo==5): 
    print "Restore in checkpoint number" + check + "."
if(whatodo==-1): 
    print "Error!"
while (whatodo != 6):
      ret=revolve.whatodo(check,capo,fine,snaps_in,info)
      oldcapo=capo
      whatodo=ret[0]
      check=ret[1]
      capo=ret[2]
      fine=ret[3]
      info=ret[4]
      if(whatodo==1): 
          print "Advance from " + str(oldcapo) + " to " + str(capo) + "."
      if(whatodo==2): 
        print "Store in checkpoint number " + str(check) + "."
      if(whatodo==3): 
        print "First turn: Initialize adjoints and reverse first step."
      if(whatodo==4): 
        print "Forward and reverse one step."
      if(whatodo==5): 
        print "Restore in checkpoint number " + str(check) + "."
      if(whatodo==-1): 
        print "Error!"
