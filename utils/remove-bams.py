#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Remove obsolete bam (index) files
"""

import sys
import os

args = sys.argv

remove_list = args[1:len(args)]

for file in remove_list:
    os.remove(file)
    
    