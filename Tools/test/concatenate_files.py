"""
A macro to concatenate PANDAS data frames
"""

import argparse
import glob
import os

PARSER = argparse.ArgumentParser(description=__doc__)

PARSER.add_argument("folderInput",
                    help="Input folder with data-frame csv files \
                          to be concatenated.")

ARGS = PARSER.parse_args()

file_names = []

if os.path.isdir(ARGS.folderInput):
    file_names = glob.glob(os.path.join(ARGS.folderInput,"df*.csv"))

WRITE_HEADER = True
with open(ARGS.folderInput.strip("/") + "_total.csv", 'w') as out_file:

    print "[" + __file__ + "] writing : {} ".format(out_file)

    for file_name in file_names:

        with open(file_name) as in_file:

            print "[" + __file__ + "] processing : {} ".format(file_name)

            in_file_line = 0
            for line in in_file:
                in_file_line += 1

                if in_file_line == 1:

                    if WRITE_HEADER:
                        WRITE_HEADER = False
                    else:
                        continue

                out_file.write(line)
                
