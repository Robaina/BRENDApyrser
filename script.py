from parseBRENDA import BRENDA

workDir = 'C:/Users/tinta/OneDrive/Documents/Projects/BRENDA'
dataFile = workDir + '/brenda_download.txt'

brenda = BRENDA(dataFile)
r = brenda.getReactions('2.2.1.1')[0]

lines = r._getDataLines('SA')
i = 4
print(lines[i])
r._extractDataLineInfo(lines[i])

KMs = r.getKMvalues()
[KM['value'] for KM in KMs['D-ribose 5-phosphate']]
