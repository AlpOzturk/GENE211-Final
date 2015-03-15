
from dis_gen_net_parser import parseDisGenNet

ALL_FILE_NAME = "all_gene_disease_associations.txt"
OUT_FILE_NAME = "diseases.txt"
DISEASE_INDEX = 1

def updateDiseaseList():
    associations = parseDisGenNet(ALL_FILE_NAME)
    diseases = set([association[DISEASE_INDEX] for association in associations])
    diseases = list(diseases)
    diseases.sort()
    outFile = open(OUT_FILE_NAME, "w")
    outFile.write("\n".join(diseases))
    outFile.close()

updateDiseaseList()