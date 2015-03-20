
import re
import sys

DISEASES_LIST = "diseases.txt"
OUTFILE_NAME = "temp.txt"

def run():
    diseasesList = list()
    diseasesFile = open(DISEASES_LIST, "r")
    for line in diseasesFile:
        diseasesList.append(line.strip())
    diseasesFile.close()
    outFile = open(OUTFILE_NAME, "w")
    while(True):
        regexStr = raw_input("Query disease list with regex (q to quit): ")
        if regexStr == "q":
            break
        for disease in diseasesList:
            if regexMatch(regexStr, disease):
                outFile.write(disease + "\n")
                print disease
    outFile.close()

def regexMatch(regexStr, disease):
    regex = re.compile(regexStr.lower())
    return regex.match(disease.lower())



# Only run when not called via import
if __name__ == "__main__":
    if regexMatch(sys.argv[1], sys.argv[2]):
        print "yes"
    print ""
