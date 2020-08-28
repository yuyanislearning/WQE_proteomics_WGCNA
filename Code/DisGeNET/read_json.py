################################################################
# This script is for reading retrieved disgenet json file

import json

def read_and_write(data, fw):
    for i in range(len(data)):
        pro = str(data[i]['uniprotid'])
        gene_dsi = str(data[i]['gene_dsi'])
        gene_dpi = str(data[i]['gene_dpi'])
        gene_pli = str(data[i]['gene_pli'])
        disease_name = str(data[i]['disease_name'])
        score = str(data[i]['score'])
        ei = str(data[i]['ei'])
        el = str(data[i]['el'])
        source = str(data[i]['source'])
        fw.write('\t'.join([pro,gene_dsi,gene_dpi,gene_pli,disease_name,score,ei,el,source]) + '\n')
    return None



with open('case_animal_models_1.json') as f:
    data = json.load(f)

fw = open('case.tsv','w')
fw.write('\t'.join(['pro','gene_dsi','gene_dpi','gene_pli','disease_name','score','ei','el','source'])+ '\n')

read_and_write(data, fw)


with open('case_Curated_1.json') as f:
    data = json.load(f)
read_and_write(data, fw)


fw.close()
