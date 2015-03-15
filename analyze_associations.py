
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

SCORE_CUTOFF = 0.1

LEAF_FONT_SIZE = 7.5

def run():
    print "Getting target disorders and regexes..."
    regexMap = getTargetDisorderMap()
    print "Parsing input data..."
    associations = parseDisGenNet(ALL_FILE_NAME)
    disorderMap, scoreMap = processAssociations(associations, regexMap)
    print "Pruning gene list"
    geneList = set()
    print "Genes per disorder:"
    for disorder in disorderMap:
        print disorder + ": " + str(len(disorderMap[disorder]))
        geneList = geneList.union(disorderMap[disorder])
    geneList = list(geneList)
    geneList.sort()
    performHeirarchicalClustering(geneList, disorderMap, scoreMap)


def performHeirarchicalClustering(geneList, disorderMap, scoreMap):
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
    pylab.savefig( "diseaseDentrogram.png" )
    pylab.cla()

#get your data into a 2d array where rows are genes, and columns 
#are conditions
#data = numpy.array(dataMatrix)

#calculate a distance matrix
#distMatrix = dist.pdist(data)

#convert the distance matrix to square form. The distance matrix 
#calculated above contains no redundancies, you want a square form 
#matrix where the data is mirrored across the diagonal.
# distSquareMatrix = dist.squareform(distMatrix)

#calculate the linkage matrix 
#linkageMatrix = hier.linkage(distSquareMatrix)

#dendro = hier.dendrogram(linkageMatrix)

#get the order of rows according to the dendrogram 
#leaves = dendro['leaves'] 

#reorder the original data according to the order of the 
#dendrogram. Note that this slice notation is numpy specific.
#It just means for every row specified in the 'leaves' array,
#get all the columns. So it effectively remakes the data matrix
#using the order from the 'leaves' array.
#transformedData = data[leaves,:]


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

run()
