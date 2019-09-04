import os
import glob
 
# Get a list of all the file paths
moncs = glob.glob('MONC_pp.o*')
texts = glob.glob('*.txt')
arrays = glob.glob('array.*')

all_files = moncs + texts + arrays

for filePath in all_files:
    try:
        os.remove(filePath)
    except:
        print("Error while deleting file : ", filePath)
