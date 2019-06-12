##################################################
# Update metadata
##################################################

import json
import fileinput
import os
import sys

base_path = '.\\'
if len(sys.argv) == 2:
    base_path = sys.argv[1].rstrip('\\') + '\\'

path_exports = base_path + 'exports\\'
path_metadata = base_path + 'metadata\\'

def load_json(path):
    with open(path) as file:
        return json.load(file)

def change_emsa_xtilt(path, xtilt):
    with fileinput.FileInput(path, inplace=True) as file:
        for line in file:
            if line.startswith("#XTILTSTGE"):
                line = "#XTILTSTGE   : " + xtilt + "\n"
            print(line, end='')

def get_index(string):
    s = ".".join(string.split('.')[:-1]).split('-')[0]
    start = -1
    for i in range(len(s)-1, -1, -1):
        if not s[i].isdigit():
            start = i
            break
    return s[start+1:]

def main():
    for filename in os.listdir(path_exports):
        if filename.endswith('.emsa') or filename.endswith('.msa'):
            emsa_path = path_exports + filename
            index = get_index(filename)
            metadata_path = path_metadata + 'meta' + index + '.json' # Replace any file ending with '.json'
            try:
                j = load_json(metadata_path)
                xtilt = j['tilt_x'].split(' ')[0]
                change_emsa_xtilt(emsa_path, xtilt)
                print('Updated tilt_x of file ' + filename)
            except FileNotFoundError:
                print('File `' + metadata_path + '` not found')

main()
