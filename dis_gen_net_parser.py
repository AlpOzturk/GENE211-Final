
DELIMITER = "\t" 
GENE_INDEX = 1
DISEASE_INDEX = 4
SCORE_INDEX = 5

def parseDisGenNet(inputFileName):
    inputFile = open(inputFileName, "r")
    associations = set()
    isFirst = True
    for line in inputFile:
        if isFirst:
            isFirst = False
        else:
            splitLine = line.strip().split(DELIMITER)
            association = (splitLine[GENE_INDEX], splitLine[DISEASE_INDEX], float(splitLine[SCORE_INDEX]))
            associations.add(association)
    inputFile.close()
    return associations
