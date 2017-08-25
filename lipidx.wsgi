import os
import sys
import site
import logging

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), 'lipid_analysis.log'))

logging.basicConfig(filename=os.path.join(os.path.dirname(os.path.realpath(__file__)), 'lipid_analysis.log'), level=logging.DEBUG)

from lipidx import app as application
