from parseBRENDA import BRENDA

workDir = 'C:/Users/tinta/OneDrive/Documents/Projects/BRENDA'
dataFile = workDir + '/brenda_download.txt'

brenda = BRENDA(dataFile)
r = brenda.reactions.get_by_id('2.7.1.40')
r
human_enzymes = brenda.reactions.filter_by_organism('Homo sapiens neanderthalensis')
len(human_enzymes)
human_enzymes[1]
human_enzymes[1].getOrganisms()

lines = r._getDataLines('PR')
for line in lines[:10]:
    print(r._extractDataLineInfo(line, numeric_value=False))
r.getSpecies()
KMs = r.getKMvalues()
[KM['value'] for KM in KMs['D-ribose 5-phosphate']]

KMs = r.getKMvalues()
KMs
KMs.filter_by_organism('Sus scrofa')

r.getKMvalues().filter_by_organism('Homo sapiens').keys()
r.getKMvalues()['pyruvate']
