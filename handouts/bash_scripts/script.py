#!/usr/bin/env python3 

import subprocess

cmd = "for file in {a..c}.txt; do echo $file; cat $file; done"

subprocess.run(cmd , shell=True)