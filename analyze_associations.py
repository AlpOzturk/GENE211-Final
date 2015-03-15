
import sys
import numpy
import pylab
import scipy
import scipy.cluster.hierarchy as hier
import scipy.spatial.distance as dist
from dis_gen_net_parser import parseDisGenNet
from regex_tester import regexMatch



ALL_FILE_NAME = "all_gene_disease_associations.txt"
TARGET_FILE_NAME = "targets.txt"
DELIMITER = "\t"

ASSOC_GENE_INDEX = 0
ASSOC_DISEASE_INDEX = 1
ASSOC_SCORE_INDEX = 2

CUTOFF_LIST = [0.0, 0.1]

OUTPUT_IMAGE_NAMES = ["diseaseDentrogram.png", "diseaseDentrogramCutoff.png",]
LEAF_FONT_SIZE = 7.5

def run():
    print "Getting target disorders and regexes..."
    regexMap = getTargetDisorderMap()
    print "Parsing input data..."
    associations = parseDisGenNet(ALL_FILE_NAME)
    for i in range(len(CUTOFF_LIST)): 
        cutoff = CUTOFF_LIST[i]
        print "Processing for cutoff: " + str(cutoff)
        disorderMap, scoreMap = processAssociations(associations, regexMap, cutoff)
        geneList = set()
        print "Genes per disorder:"
        for disorder in disorderMap:
            print disorder + ": " + str(len(disorderMap[disorder]))
            geneList = geneList.union(disorderMap[disorder])
        geneList = list(geneList)
        geneList.sort()
        print "Performing Heirarchical Clustering ..."
        performHeirarchicalClustering(geneList, disorderMap, scoreMap, i)
        print "Conducting genetic analysis ..."
        conductGeneAnalysis(geneList, disorderMap)

def conductGeneAnalysis(geneList, disorderMap):
    resultList = list()
    for gene in geneList:
        disorders = [disorder for disorder in disorderMap if gene in disorderMap[disorder]]
        resultList.append((gene, len(disorders)))
    resultList = sorted(resultList,key=lambda x: x[1], reverse=True)
    print resultList
    print len(resultList)

def getTargetDisorderMap():
    inputFile = open(TARGET_FILE_NAME, "r")
    regexMap= dict()
    for line in inputFile:
        disease, regexStr = line.strip().split(DELIMITER)
        regexMap[regexStr] = disease
    inputFile.close()
    return regexMap

def processAssociations(associations, regexMap, scoreCutoff):
    disorderMap = dict()
    scoreMap = dict() 
    for association in associations:
        gene, disease, score = association
        if score > scoreCutoff:
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

def performHeirarchicalClustering(geneList, disorderMap, scoreMap, imageNameIndex):
    dataMatrix = list() 
    disorderList = list(disorderMap.keys())
    disorderList.sort()
    for disorder in disorderList:
        dataMatrix.append([scoreMap.get((disorder, gene), 0.0) for gene in geneList])
    dataMatrix = numpy.array(dataMatrix)
    distanceMatrix = dist.squareform(dist.pdist(dataMatrix))
    linkageMatrix = hier.linkage(distanceMatrix)
    dendrogram = hier.dendrogram(linkageMatrix, labels=disorderList, leaf_font_size = LEAF_FONT_SIZE)
    #leaves = dendrogram["leaves"]
    #reorderedData = dataMatrix[leaves,:]
    pylab.savefig(OUTPUT_IMAGE_NAMES[imageNameIndex])
    pylab.cla()

run()
