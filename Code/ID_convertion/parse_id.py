#################################################################
# This script is for parsing ID retrieved from PANTHER

import json

d = {}

for i in range(108):
    with open('json/diff_pro_%s.json'%(i)) as f:
        data = json.load(f)
        mapped = data['search']['mapping']['mapped']
        if type(mapped)!=type(d):
            gene = mapped[0]['gene'].split("|")[2].split('=')[1]
            d[gene]=[]
            for m in mapped:
                target_gene = m['target_gene'].split("|")[2].split('=')[1]
                d[gene].append(target_gene)
        else:
            try:
            # if multiple gene exists
                target_gene = data['search']['mapping']['mapped']['target_gene'].split("|")[2].split('=')[1]
                gene = data['search']['mapping']['mapped']['gene'].split("|")[2].split('=')[1]
                d[gene] = [target_gene]
            except:
                continue

fw = open('diff_map.tsv','w') 
for k in d.keys():
    fw.write(k + '\t' + '\t'.join(d[k]) + '\n')

fw.close()


fw = open('diff_human_1.5.txt','w')
for k in d.keys():
    for e in d[k]:
        fw.write(e+'\n')
fw.close()
