
# The purpose of this script is to generate the list of target disorders and associated regex's.
# This could be done by hand, but this reduces the possibility of typos.

OUTPUT_FILE = "targets.txt"
DELIMITER = "\t"
SUB_DELIMITER = "|"

#def writeToFile(outFile, disorderName, disorderRegexStrs):

outFile = open(OUTPUT_FILE, "w")
writeToFile(outFile, "Autistic Disorder", [".*autis.*"])
outFile.close()

#DISORDER_LIST = ["Autistic Disorder", "Dementia", "Alzheimers", "Bipolar Disorder"]

