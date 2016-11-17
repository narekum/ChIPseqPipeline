import os
import sys
import csv

class ChromosomeEnds(dict):
    def __init__(self,  build):
        infile = os.path.expanduser("~/mount/publicdata/"+build+"/chrmSizes."+build)
        
        reader = csv.reader(open(infile), delimiter="\t")
        
        for row in reader:
            self[row[0]] = int(row[1])
