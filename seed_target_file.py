
# The purpose of this script is to generate the list of target disorders and associated regex's.
# This could be done by hand, but this reduces the possibility of typos.

TESTING_DISORDER_LIST = ["Autistic Disorder", "Dementia", "Depression", "ADHD", "Bipolar disorder", "Epilepsy", "OCD", "Schizophrenia", "Tourette's Disorder", "Alzheimer's Disease"]

OUTPUT_FILE = "targets.txt"
DELIMITER = "\t"
SUB_DELIMITER = "|"


AUTISM_REGEX = ".*autis.*"
DEMENTIA_REGEX = "(frontotemporal |Mild |Mixed |moderate |familial |severe |progressive |semantic |(pre)?senile.* |.*vascular |)dementia.*"
DEPRESSION_REGEX = "(chronic |atypical |endogenous |.*major.*|recurrent |severe |clinical |)depress.*"
ADHD_REGEX = ".*attention.*defi.*"
BIPOLAR_REGEX = ".*bipolar.*(type|disorder)"
EPILEPSY_REGEX = "(.*infantile |.*generalized.*|.*temporal |myoclonic.*|nocturnal |photogenic |reflex |refractory.*|.*benign.*|rolandic |severe.*|idiopathic |.*onset.*|.*dependent.*|)epilep.*"
OCD_REGEX = ".*compulsive.*"
SCHIZOPHRENIA_REGEX = ".*schizophrenia.*"
TOURETTES_REGEX = ".*tourette.*"
ALZHEIMERS_REGEX=".*alzheimer.*"


outFile = open(OUTPUT_FILE, "w")
outFile.write(DELIMITER.join(["Autistic Disorder", AUTISM_REGEX + "\n"]))
outFile.write(DELIMITER.join(["Dementia", DEMENTIA_REGEX + "\n"]))
outFile.write(DELIMITER.join(["Depression", DEPRESSION_REGEX + "\n"]))
outFile.write(DELIMITER.join(["ADHD", ADHD_REGEX + "\n"]))
outFile.write(DELIMITER.join(["Bipolar disorder", BIPOLAR_REGEX + "\n"]))
outFile.write(DELIMITER.join(["Epilepsy", EPILEPSY_REGEX + "\n"]))
outFile.write(DELIMITER.join(["OCD", OCD_REGEX + "\n"]))
outFile.write(DELIMITER.join(["Schizophrenia", SCHIZOPHRENIA_REGEX + "\n"]))
outFile.write(DELIMITER.join(["Tourette's Disorder", TOURETTES_REGEX + "\n"]))
outFile.write(DELIMITER.join(["Alzheimer's Disease", ALZHEIMERS_REGEX]))
outFile.close()
