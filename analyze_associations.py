
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
OUTPUT_FILE_NAMES = ["gene_out.csv", "gene_out_cutoff.csv"]
LEAF_FONT_SIZE = 7.5
BLANK_SPACE = "-"
PROB_PRECISION = 3


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
        conductGeneAnalysis(geneList, disorderMap, OUTPUT_FILE_NAMES[i])

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
    pylab.savefig(OUTPUT_IMAGE_NAMES[imageNameIndex])
    pylab.cla()

def conductGeneAnalysis(geneList, disorderMap, outFileName):
    resultList = list()
    disorderList = sorted(disorderMap.keys())
    for gene in geneList:
        disorders = [disorder for disorder in disorderMap if gene in disorderMap[disorder]]
        resultList.append((gene, disorders))
    resultList = sorted(resultList,key=lambda x: len(x[1]), reverse=True)
    maxOverlap = len(resultList[0][1])
    overlapMap = dict()
    pooledOverlap = getCounterList(maxOverlap + 1)
    for disorder in disorderList:
        # Remember, this will be 1-indexed
        overlapMap[disorder] = getCounterList(maxOverlap + 1)
    for result in resultList:
        gene, disorders = result
        for disorder in disorders:
            overlapMap[disorder][len(disorders)] += 1
        pooledOverlap[len(disorders)] += 1
    outFile = open(outFileName, 'w')
    outFile.write(DELIMITER.join(disorderList) + "\n \n")
    overlapLists = list()
    for disorder in disorderList:
        toPrint = [str(i) for i in overlapMap[disorder]]
        outFile.write(DELIMITER.join([str(i) for i in overlapMap[disorder]]) + "\n")
        overlapLists.append(overlapMap[disorder])
    outFile.write(DELIMITER.join([str(i) for i in pooledOverlap]) + "\n")
    overlapLists.append(pooledOverlap)
    for overlapList in overlapLists:
        outFile.write(DELIMITER.join(getProbList(overlapList)) + "\n")
    outFile.write("\n")
    makeConnectionMatrix(disorderList, disorderMap, outFile)
    outFile.close()

def getCounterList(length):
    return [0 for i in range(length)]

def getProbList(countList):
    total = sum(countList)
    return [str(round(float(count) / total, PROB_PRECISION)) for count in countList]

def makeConnectionMatrix(disorderList, disorderMap, outFile):
    disorderList.insert(0, " ")
    numDisorders = len(disorderList)
    for x in range(numDisorders):
        for y in range(numDisorders):
            if y == 0:
                outFile.write(disorderList[x])
            elif x == 0:
                outFile.write(disorderList[y])
            elif y >= x:                 
                outFile.write(BLANK_SPACE)
            else:
                genesX = disorderMap[disorderList[x]]
                genesY = disorderMap[disorderList[y]]
                genesCommon = genesX.intersection(genesY)
                outFile.write(str(len(genesCommon)))
            outFile.write(DELIMITER)
        outFile.write("\n")
run()
