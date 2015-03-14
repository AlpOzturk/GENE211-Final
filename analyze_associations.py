
import sys
from dis_gen_net_parser import parseDisGenNet

ALL_FILE_NAME = "all_gene_disease_associations.txt"

def run():
    print "Parsing input data."
    associations = parseDisGenNet(ALL_FILE_NAME)

run()