"""
Add a header line beginning with ">" to all files in
a file listing files prefixes. I tried to implement
a linux principle of havinga tool to do a well a
specific job. This is for adding headers to fasta files
only.
"""

import os
import argparse

# Parse the command line options
parser = argparse.ArgumentParser(description='Insert a line at the beginning of each file.')
parser.add_argument('-s', '--suffix', help='Suffix for the new file', required=True)
parser.add_argument('-n', '--names_file', help='File containing the list of names', required=True)
args = parser.parse_args()
suffix = args.suffix
names_file = args.names_file

# Get the current working directory
folder_path = os.getcwd()

# Read the names from the provided file
with open(names_file, 'r') as file:
    names = file.read().splitlines()

# Iterate through the files in the folder
for name in names:
    print(f"Adding header to {name}{suffix}")
    file_path = os.path.join(folder_path, f"{name}{suffix}")
    if os.path.exists(file_path):
        with open(file_path, 'r+') as file:
            content = file.read()
            file.seek(0, 0)
            file.write(f">{name}\n{content}\n")
