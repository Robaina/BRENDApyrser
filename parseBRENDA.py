import re
import numpy as np

fields = {
    'AC': 'activating compound',
    'AP': 'application',
    'CF': 'cofactor',
    'CL': 'cloned',
    'CR': 'crystallization',
    'EN': 'engineering',
    'EXP': 'expression',
    'GI': 'general information on enzyme',
    'GS': 'general stability',
    'IC50': 'IC-50 Value',
    'ID': 'EC-class',
    'IN': 'inhibitors',
    'KKM': 'Kcat/KM-Value substrate in {...}',
    'KI': 'Ki-value, inhibitor in {...}',
    'KM': 'KM-value, substrate in {...}',
    'LO': 'localization',
    'ME': 'metals/ions',
    'MW': 'molecular weight',
    'NSP': 'natural substrates/products	reversibilty information in {...}',
    'OS': 'oxygen stability',
    'OSS': 'organic solvent stability',
    'PHO': 'pH-optimum',
    'PHR': 'pH-range',
    'PHS': 'pH stability',
    'PI': 'isoelectric point',
    'PM': 'posttranslation modification',
    'PR': 'protein',
    'PU': 'purification',
    'RE': 'reaction catalyzed',
    'RF': 'references',
    'REN': 'renatured',
    'RN': 'accepted name (IUPAC)',
    'RT': 'reaction type',
    'SA': 'specific activity',
    'SN': 'synonyms',
    'SP': 'substrates/products, reversibilty information in {...}',
    'SS': 'storage stability',
    'ST': 'source/tissue',
    'SU': 'subunits',
    'SY': 'systematic name',
    'TN': 'turnover number, substrate in {...}',
    'TO': 'temperature optimum',
    'TR': 'temperature range',
    'TS': 'temperature stability'
}


class BRENDA:
    """
    Provides methods to parse the BRENDA database (https://www.brenda-enzymes.org/)
    """
    def __init__(self, path_to_database):

        with open(path_to_database, encoding="iso-8859-1") as file:
            self.data = file.read()
        self.__ec_numbers = [ec.group(1)
                             for ec in re.finditer('(?<=ID\\t)(.*)(?=\\n)', self.data)]

    def __getRxnData(self, ec_number):
        idx_1 = self.data.find(f'ID\t{ec_number}')
        idx_2 = self.data[idx_1:].find('///')
        return self.data[idx_1:(idx_1+idx_2)]

    def getReactions(self, ec_number: str = None) -> list:
        """
        Get reaction object corresponding to the given EC number or all reaction objects
        in the database if left as None.
        """
        if ec_number is None:
            return [Reaction(self.__getRxnData(ec)) for ec in self.__ec_numbers]
        else:
            if ec_number in self.__ec_numbers:
                return [Reaction(self.__getRxnData(ec_number))]
            else:
                raise ValueError(f'Enzyme {ec_number} not found in database')

    def getOrganisms(self) -> list:
        """
        Get list of represented species in BRENDA
        """

    def getMetaInformation(self) -> dict:
        """
        Get BRENDA version, licence and meta information
        """


class Reaction:
    def __init__(self, reaction_data):
        self.__reaction_data = reaction_data
        self.ec_number = re.search('(?<=ID\t)(.*)(?=\n)', self.__reaction_data).group(1)
        self.systematic_name = re.search(
            '(?<=SN\t)(.*)(?=\n)', self.__reaction_data).group(1)
        self.name = re.search('(?<=RN\t)(.*)(?=\n)', self.__reaction_data).group(1)
        self.mechanism = re.search('(?<=RE\t)(.*)(?=\n)', self.__reaction_data).group(1)
        self.reaction_type = re.search('(?<=RT\t)(.*)(?=\n)', self.__reaction_data).group(1)

    def _getDataLines(self, pattern: str):
        try:
            search_pattern = f'{pattern}\t(.+?)\n(?!\t)'
            return [p.group(1)
                    for p in re.finditer(
                        search_pattern, self.__reaction_data, flags=re.DOTALL)]
        except Exception:
            return []

    @staticmethod
    def __removeTabs(line):
        return line.replace('\n', '').replace('\t', '').strip()

    @staticmethod
    def __extractDataField(line, regex_tags: tuple):
        try:
            l, r = regex_tags
            searched_s = re.search(f'{l}(.+?){r}', line)
            span = searched_s.span()
            matched_s = line[span[0] + 1:span[1] - 1].strip()
            line = line.replace(f'{searched_s.group()}', '')
            return (line, matched_s)
        except Exception:
            return (line, '')

    def _extractDataLineInfo(self, line: str, numeric_value=False):
        """
        Extracts data fields in each data line according to the tags used by BRENDA
        and described in the REAMDE.txt file. What remains after extracting all tags
        is the value of that particular data field, e.g., KM value.
        """
        line = self.__removeTabs(line)
        line, meta = self.__extractDataField(line, ('\(', '.*\)'))
        line, refs = self.__extractDataField(line, ('<', '>'))
        line, species = self.__extractDataField(line, ('#', '#'))
        line, specific_info = self.__extractDataField(line, ('{', '.*}'))
        if numeric_value:
            value = float(line.strip())
        else:
            value = line.strip()

        return {'value': value, 'species': species.split(','),
                'meta': meta, 'refs': refs.split(','),
                'specific_info': specific_info}

    def _getDictOfEnzymeActuators(self, pattern: str) -> dict:
        res = {}
        lines = self._getDataLines(pattern)
        for line in lines:
            data = self._extractDataLineInfo(line)
            if data['value'] != 'more':
                res[data['value']] = {'species': data['species'],
                                      'meta': data['meta'],
                                      'refs': data['refs']}
        return res

    def _getDictOfEnzymeProperties(self, pattern):
        res = {}
        lines = self._getDataLines(pattern)
        for line in lines:
            data = self._extractDataLineInfo(line, numeric_value=True)
            substrate = data['specific_info']
            if substrate != 'more':
                if substrate not in res.keys():
                    res[substrate] = []
                res[substrate].append({'value': data['value'],
                                       'species': data['species'],
                                       'meta': data['meta'],
                                       'refs': data['refs']})
        return res

    def getCofactors(self):
        return self._getDictOfEnzymeActuators('CF')

    def getMetals(self):
        return self._getDictOfEnzymeActuators('ME')

    def getInhibitors(self):
        return self._getDictOfEnzymeActuators('IN')

    def getActivators(self):
        return self._getDictOfEnzymeActuators('AC')

    def getKMvalues(self):
        return self._getDictOfEnzymeProperties('KM')

    def getKIvalues(self):
        return self._getDictOfEnzymeProperties('KI')

    def getKKMvalues(self):
        return self._getDictOfEnzymeProperties('KKM')

    def getKcatvalues(self):
        return self._getDictOfEnzymeProperties('TN')

    def getSpecificActivities(self):
        lines = self._getDataLines('SA')
        return [self._extractDataLineInfo(line, numeric_value=True) for line in lines]

    def getSubstratesAndProducts(self) -> list:
        """
        Returns list of dicts with evaluated "natural" substrates and products
        of the enzyme across organisms.
        """
        substrates, products, res = [], [], []
        lines = self._getDataLines('NSP')
        for line in lines:
            data = self._extractDataLineInfo(line)
            is_full_rxn = '=' in data['value']
            rxn = data['value'].replace(
                '{}', '').replace('?', '').replace('more', '').strip()
            if is_full_rxn:
                subs, prods = rxn.split('=')
                subs = [s.strip() for s in subs.split('+') if s.strip() != '']
                prods = [s.strip() for s in prods.split('+') if s.strip() != '']
                subs.sort()
                prods.sort()
                if (subs not in substrates and len(subs) > 0 and len(prods) > 0):
                    substrates.append(subs)
                    products.append(prods)
                    res.append({'substrates': subs, 'products': prods})
        return res

    def getTemperatureData(self):

        def getData(data_type):
            values = []
            lines = self._getDataLines(data_type)
            if data_type != 'TR':
                eval_temp = lambda t: float(t) if not re.search('\d-\d', t) else np.mean(
                    [float(s) for s in t.split('-')])
            else:
                eval_temp = lambda t: [float(s) for s in t.split('-')]

            for line in lines:
                res = self._extractDataLineInfo(line)
                values.append({'value': eval_temp(res['value']),
                               'species': res['species'],
                               'meta': res['meta'], 'refs': res['refs']})
            return values

        return {'optimum': getData('TO'),
                'range': getData('TR'),
                'stability': getData('TS')}

    def getSpecies(self) -> dict:
        """
        Returns a dict listing all proteins for given EC number
        """
        species = {}
        lines = self._getDataLines('PR')
        for line in lines:
            res = self._extractDataLineInfo(line)
            species[res['species'][0]] = {'name': res['value'], 'refs': res['refs']}
        return species

    def getReferences(self):
        """
        Returns a dict listing the bibliography cited for the given EC number
        """
        references = {}
        lines = self._getDataLines('RF')
        for line in lines:
            line = self.__removeTabs(line)
            line, refs = self.__extractDataField(line, ('<', '>'))
            references[refs[0]] = line
        return references
