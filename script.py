from parseBRENDA import BRENDA

workDir = 'C:/Users/tinta/OneDrive/Documents/Projects/BRENDA'
dataFile = workDir + '/brenda_download.txt'

brenda = BRENDA(dataFile)
r = brenda.getReactions('1.1.1.1')[0]

r.getPHData()

lines = r._getDataLines('PR')
for line in lines[:10]:
    # print(line)
    print(r._extractDataLineInfo(line, numeric_value=False))
r.getSpecies()
KMs = r.getKMvalues()
[KM['value'] for KM in KMs['D-ribose 5-phosphate']]

s = "#102# 8.47 {(S)-N-benzyl-3-pyrrolidinol}  (#102# pH 8.0 <185>) <185>"
r._extractDataLineInfo(line, numeric_value=False)

r._getDictOfEnzymeProperties('KM')

r.getKMvalues()
