
import re

DISEASES_LIST = "diseases.txt"

def run():
    diseasesList = list()
    diseasesFile = open(DISEASES_LIST, "r")
    for line in diseasesFile:
        diseasesList.append(line.strip())
    diseasesFile.close()
    while(True):
        regexStr =  raw_input("Query disease list with regex (q to quit): ")
        if regexStr == "q":
            break
        for disease in diseasesList:
            if regexMatch(regexStr, disease):
                print disease


def regexMatch(regexStr, disease):
    regex = re.compile(regexStr.lower())
    return regex.match(disease.lower())



# Only run when not called via import
if __name__ == "__main__":
    run()