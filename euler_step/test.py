#! /usr/bin/env python
import os
import sys

def run_cmd():
    NX = [1, 2, 3]
    NY = [1, 2, 3]
    for i in NX:
        for j in NY:
            cmd = 'echo %s %s' %(i, j)
            os.system(cmd)

def get_data():
    fout = fopen("1.log", "r")

def main(argv=None):
    run_cmd()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
