###################################################################
# This script is to get protein id from Panther API
# Change 'curl' to the corresponding download command of your platform.

import os
import time

# convert human protein to mouse protein
# read input protein list
diff_pro = []
with open('diff_pro_names_1.5.txt','r') as f:
    for line in f:
        line = line.strip()
        diff_pro.append(line)

# retrieve ID from PANTHER, change organism and targetOrganism accordingly 
for i, q in enumerate(diff_pro):
    os.system('curl -X POST "http://pantherdb.org/services/oai/pantherdb/ortholog/matchortho?geneInputList=%s&organism=10090&targetOrganism=9606&orthologType=all" -H "accept: application/json" > ./json/diff_pro_%s.json'%(q, i))
    time.sleep(15)
