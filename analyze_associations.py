
import sys
from dis_gen_net_parser import parseDisGenNet
from regex_tester import regexMatch

ALL_FILE_NAME = "all_gene_disease_associations.txt"
TARGET_FILE_NAME = "targets.txt"
DELIMITER = "\t"

ASSOC_GENE_INDEX = 0
ASSOC_DISEASE_INDEX = 1
ASSOC_SCORE_INDEX = 2

SCORE_CUTOFF = 0.1

def run():
    print "Getting target disorders and regexes..."
    regexMap = getTargetDisorderMap()
    #targetDisorders = disorderMap.values()
    print "Parsing input data..."
    associations = parseDisGenNet(ALL_FILE_NAME)
    disorderMap, scoreMap = processAssociations(associations, regexMap)
    print "Pruning gene list"
    pruneGenes(disorderMap)
    for disorder in disorderMap:
        print disorder + ": " + str(len(disorderMap[disorder]))

def getTargetDisorderMap():
    inputFile = open(TARGET_FILE_NAME, "r")
    regexMap= dict()
    for line in inputFile:
        disease, regexStr = line.strip().split(DELIMITER)
        regexMap[regexStr] = disease
    inputFile.close()
    return regexMap

def processAssociations(associations, regexMap):
    disorderMap = dict()
    scoreMap = dict() 
    for association in associations:
        gene, disease, score = association
        if score > SCORE_CUTOFF:
            for regexStr in regexMap:
                if regexMatch(regexStr, disease):
                    actualDisorder = regexMap[regexStr]
                    genesOfDisorder = disorderMap.get(actualDisorder, set())
                    genesOfDisorder.add(gene)
                    disorderMap[actualDisorder] = genesOfDisorder
                    scoreMapKey = (actualDisorder, gene)
                    scoreMap[scoreMapKey] = max(scoreMap.get(scoreMapKey, 0.0), score)
                    break
    return disorderMap, scoreMap

def pruneGenes(geneMap):
    pass

run()