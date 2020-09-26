#!/usr/bin/env python

#############################################################################
##  this file is part of the FastRoot package
##  see LICENSE for terms and conditions of usage.
#############################################################################

PROGRAM_NAME = "FastRoot"
PROGRAM_AUTHOR = ["Uyen Mai","Merve Kilic","Erfan Sayyari","Siavash Mirarab"]
PROGRAM_LICENSE = "MIT License"
PROGRAM_VERSION = "1.4"
PROGRAM_YEAR = "2017"
PROGRAM_INSTITUTE = "University of California at San Diego"

import logging
from sys import stdout
    
def new_logger(myName,myLevel=logging.INFO,myStream=stdout):
    logger = logging.getLogger(myName)
    logger.setLevel(myLevel)
    handler = logging.StreamHandler(myStream)
    formatter = logging.Formatter('%(levelname)s:%(name)s:%(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.propagate = False
    
    return logger
