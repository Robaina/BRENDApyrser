import re
import numpy as np


class reaction:

    def __init__(self, rxn_data):
        self.data = rxn_data



def findRxnData(data, ec_number):
    idx_1 = data.find(f'ID\t{ec_number}')
    idx_2 = data[idx_1:].find('///')
    return data[idx_1:idx_2]



class BRENDA:
    """
    Provides methods to parse the BRENDA database (https://www.brenda-enzymes.org/)
    """
    def __init__(self, path_to_database):
        def getECnumberIndices() -> dict:
            EC_pattern = '(?<=ID\\t)(.*)(?=\\n)'

#             EC_pattern = '(?<=ID\\t)(.*)(?=///\\nID\\t)'
#             EC_pattern = 'ID\t(.+?)///\nID\t'
#             return {self.data[m.start():m.end()]: (m.start(), m.end())
#                     for m in re.finditer(EC_pattern, self.data)}

            return {self.data[m.start():m.end()]: (z.start(), z.end())
                    for m, z in zip(re.finditer(EC_pattern, self.data),
                                   re.finditer(rxn_data_pattern, self.data)
                                   )
                   }


        with open(path_to_database, encoding="iso-8859-1") as file:
            self.data = file.read()
        self.ECIndices = getECnumberIndices()


    def getReactionDataIndices(self) -> dict:
        """
        Returns a dict with keys equal to EC numbers and values corresponding
        to the start and end indices of the data corresponding to that reaction.
        """
        '///\nID\t1.1.1.1'

    def getProteins(self, ec_number) -> dict:
        """
        Returns a dict listing all proteins for given EC number
        """
        def getPRlines(search_indices):
            return [p.group(1)
                    for p in re.finditer("PR\t(.+?)\nPR",
                                         self.data[search_indices[0]:search_indices[1]]
                                        )
                   ]


    def getReferences(self, ec_number: str) -> dict:
        """
        Returns a dict listing all references for given EC number
        """

    def getReactionName(self, ec_number: str) -> str:
        """
        Returns a the systematic name of the reaction
        """
        return self.data[self.ECIndices[ec_number][0]:].find('(?<=SN\\t)(.*)(?=\\n)')



    def getNextEnzymeIdx(self, ec_number: str) -> int:
        try:
            key_idx = list(self.ECIndices).index(ec_number)
            if key_idx == len(self.ECIndices):
                return len(self.data)
            else:
                return self.ECIndices[list(self.ECIndices.keys())[key_idx + 1]][0]
        except:
            raise ValueError("Ec number not in database")

    def getKMvalues(self, ec_number: str, substrate: str=None) -> dict:
        """
        Returns a dictionary with all KM values of the enzyme with
        given EC number. If a substrate is given, then results are
        restricted to that substrate.
        """
        search_indices = (self.ECIndices[ec_number][1], self.getNextEnzymeIdx(ec_number))

        def getEnzymeSubstrates(KM_lines):
            substrates = []
            for line in KM_lines:
                sub = extractKMInfo(line)['substrate']
                if sub not in substrates:
                    substrates.append(sub)
            return substrates

        def extractKMInfo(KM_line):
            res = {}
            try:
                species = re.search('#(.+?)#', KM_line).group(1).split(',')
                res['species'] = species
            except:
                res['species'] = ''
            try:
                KM_value = re.search('# (.+?) {', KM_line).group(1)
                res['KM'] = KM_value
            except:
                res['KM'] = np.nan
            try:
                substrate = re.search('{(.+?)}', KM_line).group(1)
                res['substrate'] = substrate
            except:
                res['substrate'] = ''
            try:
                meta = re.search('\((.+?)\)', KM_line).group(1)
                res['meta'] = meta
            except:
                res['meta'] = ''
            try:
                references = re.search('<(.+?)>', KM_line).group(1)
                res['references'] = references
            except:
                res['references'] = ''
            return res

        def getKMlines(search_indices):
            return [p.group(1) for p in re.finditer("KM\t(.+?)\nKM",
                                                    self.data[search_indices[0]:search_indices[1]])]

        KM_lines = getKMlines(search_indices)
        enzyme_substrates = getEnzymeSubstrates(KM_lines)
        if substrate is None:
            res = {s: [] for s in enzyme_substrates}
            for line in KM_lines:
                KM_info = extractKMInfo(line)
                res[KM_info['substrate']].append(float(KM_info['KM']))
            return res
        else:
            res = {substrate: []}
            for line in KM_lines:
                KM_info = extractKMInfo(line)
                if KM_info['substrate'] == substrate:
                    res[substrate].append(float(KM_info['KM']))
            return res











































'''
Provides classes and functions to parse the BRENDA data base
'''
# import re
# import numpy as np
# WorkDir = 'C:/Users/robaina/OneDrive/Documents/BRENDA'
# DataFile = WorkDir + '/brenda_download.txt'

# with open(DataFile, encoding="iso-8859-1") as file:
#     data = file.read()

# ECidx = [( m.start(), m.end() )
#          for m in re.finditer('(?<=ID\\t)(.*)(?=\\n)', data)] # Enzyme number

# MEidx = [( m.start(), m.end() )
#       for m in re.finditer('(?<=ME\\t)(.*)(?=\\n)', data)] # Metals/ions

# CFidx = [( m.start(), m.end() )
#     for m in re.finditer('(?<=CF\\t)(.*)(?=\\n)', data)] # Cofactors

# i=2000;data[ECidx[i][0]:ECidx[i][1]]
# data[MEidx[i][0]:MEidx[i][1]]
# data[MEidx[i][0]:CFidx[i][1]]
