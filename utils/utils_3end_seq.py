#!/usr/bin/env python3

import glob
import os
import sys
import subprocess
import shutil



script_dir = os.path.abspath(os.path.dirname(__file__))
script_dir = os.path.dirname(script_dir)
sys.path.append(os.path.join(script_dir, "utils"))
import utils_general as utils