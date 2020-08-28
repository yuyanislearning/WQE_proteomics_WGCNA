##############################################################
# This script is for getting api from disgenet

import os
import time



pros = []
with open('/Users/yuyan/Desktop/Projects/WQE/CaseOLAP/top_caseolap_pro_id.txt') as f:
    for line in f:
        pros.append(line.strip())

pros = '%2C'.join(pros)


os.system('curl -X GET "https://www.disgenet.org/api/gda/gene/uniprot/%s?source=CURATED&disease_class=C14" -H "accept: */*" > case_Curated_1.json'%(pros))
time.sleep(15)
os.system('curl -X GET "https://www.disgenet.org/api/gda/gene/uniprot/%s?source=ANIMAL_MODELS&disease_class=C14" -H "accept: */*" > case_animal_models_1.json'%(pros))
