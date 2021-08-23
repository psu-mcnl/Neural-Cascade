import os
import sys
import inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
projectdir = os.path.dirname(currentdir)
sys.path.insert(0, projectdir) 
