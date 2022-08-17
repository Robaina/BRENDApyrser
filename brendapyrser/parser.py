from __future__ import annotations
import re
import numpy as np
import pandas as pd

__version__ = '0.0.1'
__author__ = 'Semidán Robaina Estévez, 2020'


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

units = {
    'KM': 'mM',
    'KI': 'mM',
    'TN': '$s^{-1}$',
    'SA': '$µmol.min^{-1}.mg^{-1}$',
    'KKM': '$mM^{-1}.s^{-1}$',
    'TO': '${}^oC$',
    'TR': '${}^oC$',
    'TS': '${}^oC$',
    'MW': 'Da'
}


class BRENDA:
    """
    Provides methods to parse the BRENDA database (https://www.brenda-enzymes.org/)
    """
    def __init__(self, path_to_database):

        with open(path_to_database, encoding="iso-8859-1") as file:
            self.__data = file.read()
        self.__ec_numbers = [ec.group(1)
                             for ec in re.finditer('(?<=ID\\t)(.*)(?=\\n)', self.__data)]
        self.__reactions = self.__initializeReactionObjects()
        self.__copyright = ("""Copyrighted by Dietmar Schomburg, Techn. University
        Braunschweig, GERMANY. Distributed under the License as stated
        at http:/www.brenda-enzymes.org""")
        self.__fields = fields
        self.__units = units

    def _repr_html_(self):
        """This method is executed automatically by Jupyter to print html!"""
        return """
        <table>
            <tr>
                <td><strong>Number of Enzymes</strong></td><td>{n_ec}</td>
            </tr><tr>
                <td><strong>BRENDA copyright</strong></td><td>{cr}</td>
            </tr><tr>
                <td><strong>Parser version</strong></td><td>{parser}</td>
            </tr><tr>
                <td><strong>Author</strong></td><td>{author}</td>
            </tr>
        </table>
        """.format(n_ec=len(self.__reactions),
                   cr=self.__copyright,
                   parser=__version__,
                   author=__author__)

    def __getRxnData(self):
        rxn_data = [r.group(0)
                    for r in re.finditer('ID\\t(.+?)///', self.__data, flags=re.DOTALL)]
        del self.__data
        return rxn_data

    def __initializeReactionObjects(self):
        return [Reaction(datum) for datum in self.__getRxnData()]

    @property
    def fields(self):
        return self.__fields

    @property
    def units(self):
        return self.__units

    @property
    def reactions(self):
        return ReactionList(self.__reactions)

    @property
    def copyright(self):
        return self.__copyright

    def getOrganisms(self) -> list:
        """
        Get list of all represented species in BRENDA
        """
        species = set()
        for rxn in self.__reactions:
            species.update([s['name'] for s in rxn.proteins.values()])
        species.remove('')
        species = list(set([s for s in species if 'no activity' not in s]))
        return species

    def getKMcompounds(self) -> list:
        """
        Get list of all substrates in BRENDA with KM data
        """
        cpds = set()
        for rxn in self.__reactions:
            cpds.update([s for s in rxn.KMvalues.keys()])
        try:
            cpds.remove('')
        except Exception:
            pass
        return list(cpds)


class ReactionList(list):
    # Make ReactionList slicing return ReactionList object
    def __init__(self, seq=None):
        super(self.__class__, self).__init__(seq)

    def __getslice__(self, start, stop):
        return self.__class__(super(self.__class__, self).__getslice__(start, stop))

    def __getitem__(self, key):
        if isinstance(key, slice):
            return self.__class__(super(self.__class__, self).__getitem__(key))
        else:
            return super(self.__class__, self).__getitem__(key)

    def get_by_id(self, id: str):
        try:
            return [rxn for rxn in self if rxn.ec_number == id][0]
        except Exception:
            raise ValueError(f'Enzyme with EC {id} not found in database')

    def get_by_name(self, name: str):
        try:
            return [rxn for rxn in self if rxn.name.lower() == name.lower()][0]
        except Exception:
            raise ValueError(f'Enzyme {name} not found in database')

    def filter_by_substrate(self, substrate: str) -> list[Reaction]:
        """
        Filter reactions by a specific substrate
        """
        return [
            rxn
            for rxn in self
            if any(
                [substrate in mets["substrates"]
                for mets in rxn.substratesAndProducts]
            )
            ]

    def filter_by_product(self, product: str) -> list[Reaction]:
        """
        Filter reactions by a specific product
        """
        return [
            rxn
            for rxn in self
            if any(
                [product in mets["products"]
                for mets in rxn.substratesAndProducts]
            )
            ]

    def filter_by_compound(self, compound: str) -> list[Reaction]:
        """
        Filter reactions by a substrate or product
        """
        return [
            rxn
            for rxn in self
            if any(
                [(compound in mets["substrates"] 
                or compound in mets["products"])
                for mets in rxn.substratesAndProducts]
            )
            ]

    def filter_by_organism(self, species: str):
        def is_contained(p, S): return any([p in s.lower() for s in S])
        return self.__class__(
                            [rxn for rxn in self if is_contained(species.lower(), rxn.organisms)]
                            )
        


class EnzymeDict(dict):
    def filter_by_organism(self, species: str):
        filtered_dict = {}
        def is_contained(p, S): return any([p in s for s in S])
        for k in self.keys():
            filtered_values = [v for v in self[k] if is_contained(species, v['species'])]
            if len(filtered_values) > 0:
                filtered_dict[k] = filtered_values
        return self.__class__(filtered_dict)

    def get_values(self):
        return [v['value'] for k in self.keys() for v in self[k]]


class EnzymePropertyDict(EnzymeDict):
    def filter_by_compound(self, compound: str):
        try:
            return self.__class__({compound: self[compound]})
        except Exception:
            return self.__class__({compound: []})
            # raise KeyError(f'Invalid compound, valid compounds are: {", ".join(list(self.keys()))}')


class EnzymeConditionDict(EnzymeDict):
    def filter_by_condition(self, condition: str):
        try:
            return self.__class__({condition: self[condition]})
        except Exception:
            raise KeyError(f'Invalid condition, valid conditions are: {", ".join(list(self.keys()))}')


class Reaction:
    def __init__(self, reaction_data):
        self.__reaction_data = reaction_data
        self.__ec_number = self.__extractRegexPattern('(?<=ID\t)(.*)(?=\n)')
        self.__systematic_name = self.__extractRegexPattern('(?<=SN\t)(.*)(?=\n)')
        self.__name = self.__extractRegexPattern('(?<=RN\t)(.*)(?=\n)').capitalize()
        self.__mechanism_str = (self.__extractRegexPattern('(?<=RE\t)(.*)(?=\n\nREACTION_)',
                                                       dotall=True).replace('=', '<=>')
                               .replace('\n\t', ''))
        self.__reaction_type = self.__extractRegexPattern('(?<=RT\t)(.*)(?=\n)').capitalize()
        self.__proteins = self.getSpeciesDict()
        self.__references = self.getReferencesDict()

    def getSpeciesDict(self) -> dict:
        """
        Returns a dict listing all proteins for given EC number
        """
        species = {}
        lines = self.__getDataLines('PR')
        for line in lines:
            res = self.extractDataLineInfo(line)
            species_name, protein_ID = self.__splitSpeciesFromProteinID(res['value'])
            species[res['species'][0]] = {'name': species_name,
                                          'proteinID': protein_ID,
                                          'refs': res['refs']}
        return species

    def getReferencesDict(self):
        """
        Returns a dict listing the bibliography cited for the given EC number
        """
        references = {}
        lines = self.__getDataLines('RF')
        for line in lines:
            line = self.__removeTabs(line)
            line, refs = self.__extractDataField(line, ('<', '>'))
            references[refs[0]] = line
        return references

    def printReactionSummary(self):
        data = {'EC number': self.__ec_number,
                'Name': self.__name,
                'Systematic name': self.__systematic_name,
                'Reaction type': self.__reaction_type,
                'Mechanism': self.__mechanism}
        return pd.DataFrame.from_dict(data, orient='index', columns=[''])

    def _repr_html_(self):
        """This method is executed automatically by Jupyter to print html!"""
        return """
        <table>
            <tr>
                <td><strong>Enzyme identifier</strong></td><td>{ec}</td>
            </tr><tr>
                <td><strong>Name</strong></td><td>{name}</td>
            </tr><tr>
                <td><strong>Systematic name</strong></td><td>{sys_name}</td>
            </tr><tr>
                <td><strong>Reaction type</strong></td><td>{rxn_type}</td>
            </tr><tr>
                <td><strong>Reaction</strong></td><td>{rxn_str}</td>
            </tr>
        </table>
        """.format(ec=self.__ec_number,
                   name=self.__name,
                   sys_name=self.__systematic_name,
                   rxn_type=self.__reaction_type,
                   rxn_str=self.reaction_str)

    def __extractRegexPattern(self, pattern, dotall=False):
        if dotall:
            flag = re.DOTALL
        else:
            flag = 0
        try:
            return re.search(pattern, self.__reaction_data, flags=flag).group(1)
        except Exception:
            return ''

    def __getDataLines(self, pattern: str):
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

    @staticmethod
    def __eval_range_value(v):
        try:
            if not re.search('\d-\d', v):
                return float(v)
            else:
                return np.mean([float(s) for s in v.split('-')])
        except Exception:
            return -999

    @staticmethod
    def __splitSpeciesFromProteinID(line):
        try:
            idx = re.search('[A-Z]{1}[0-9]{1}', line).start()
            return (line[:idx].strip(), line[idx:].strip())
        except Exception:
            return (line.strip(), '')

    def extractDataLineInfo(self, line: str, numeric_value=False):
        """
        Extracts data fields in each data line according to the tags used by BRENDA
        and described in the REAMDE.txt file. What remains after extracting all tags
        is the value of that particular data field, e.g., KM value.
        """
        line = self.__removeTabs(line)
        line, specific_info = self.__extractDataField(line, ('{', '.*}'))
        line, meta = self.__extractDataField(line, ('\(', '.*\)'))
        line, refs = self.__extractDataField(line, ('<', '>'))
        line, species = self.__extractDataField(line, ('#', '#'))
        if numeric_value:
            value = self.__eval_range_value(line.strip())
        else:
            value = line.strip()
        return {'value': value, 'species': species.split(','),
                'meta': meta, 'refs': refs.split(','),
                'specific_info': specific_info}

    def __extractReactionMechanismInfo(self, line: str):
        """
        Extracts reaction string and mechanism info
        """
        line = self.__removeTabs(line)
        line, meta = self.__extractDataField(line, ('\(', '.*\)'))
        rxn_str = line.strip()
        meta_list = []
        for meta_line in meta.split(';'):
            meta_line, refs = self.__extractDataField(meta_line, ('<', '>'))
            meta_line, species = self.__extractDataField(meta_line, ('#', '#'))
            meta_list.append({'species': species.split(','),
                              'refs': refs.split(','),
                              'meta': meta_line.strip()})
        return (rxn_str, meta_list)

    def __getBinomialNames(self, species_list: list) -> list:
        """
        Returns a list with binomial names mapped to the species codes
        employed by BRENDA to attach species to protein entries
        """
        species_dict = self.__proteins
        return list(set([species_dict[s]['name'] for s in species_list
                         if s in species_dict.keys()]))

    def __getFullReferences(self, refs_list: list) -> list:
        """
        Returns a list with full reference mapped to the refs codes
        employed by BRENDA in each entry
        """
        refs_dict = self.__references
        return [refs_dict[s] for s in refs_list if s in refs_dict.keys()]

    def __getDictOfEnzymeActuators(self, pattern: str) -> dict:
        res = {}
        lines = self.__getDataLines(pattern)
        for line in lines:
            data = self.extractDataLineInfo(line)
            if data['value'] != 'more':
                res[data['value']] = {'species': self.__getBinomialNames(data['species']),
                                      'meta': data['meta'],
                                      #'refs': data['refs']}
                                      'refs': self.__getFullReferences(data['refs'])}
        return EnzymePropertyDict(res)

    def __getDictOfEnzymeProperties(self, pattern: str) -> dict:
        res = {}
        lines = self.__getDataLines(pattern)
        for line in lines:
            data = self.extractDataLineInfo(line, numeric_value=True)
            substrate = data['specific_info']
            if substrate != 'more':
                if substrate not in res.keys():
                    res[substrate] = []
                res[substrate].append({'value': data['value'],
                                       'species': self.__getBinomialNames(data['species']),
                                       'meta': data['meta'],
                                       #'refs': data['refs']})
                                       'refs': self.__getFullReferences(data['refs'])})
        return EnzymePropertyDict(res)

    def __extractTempOrPHData(self, data_type: str) -> list:
        values = []
        lines = self.__getDataLines(data_type)
        if 'R' not in data_type:
            eval_value = self.__eval_range_value
        else:
            def eval_value(v):
                try:
                    return [float(s) for s in v.split('-')]
                except Exception:
                    return [-999, -999]

        for line in lines:
            data = self.extractDataLineInfo(line)
            values.append({'value': eval_value(data['value']),
                           'species': self.__getBinomialNames(data['species']),
                           'meta': data['meta'],
                           'refs': data['refs']})
        return values

    @property
    def summary(self):
        return self.printReactionSummary()

    @property
    def ec_number(self):
        return self.__ec_number

    @property
    def name(self):
        return self.__name

    @property
    def systematic_name(self):
        return self.__systematic_name

    @property
    def reaction_str(self):
        return self.__extractReactionMechanismInfo(self.__mechanism_str)[0]

    @property
    def mechanism(self):
        return self.__extractReactionMechanismInfo(self.__mechanism_str)[1]

    @property
    def reaction_type(self):
        return self.__reaction_type

    @property
    def cofactors(self):
        return self.__getDictOfEnzymeActuators('CF')

    @property
    def metals(self):
        return self.__getDictOfEnzymeActuators('ME')

    @property
    def inhibitors(self):
        return self.__getDictOfEnzymeActuators('IN')

    @property
    def activators(self):
        return self.__getDictOfEnzymeActuators('AC')

    @property
    def KMvalues(self):
        return self.__getDictOfEnzymeProperties('KM')

    @property
    def KIvalues(self):
        return self.__getDictOfEnzymeProperties('KI')

    @property
    def KKMvalues(self):
        return self.__getDictOfEnzymeProperties('KKM')

    @property
    def Kcatvalues(self):
        return self.__getDictOfEnzymeProperties('TN')

    @property
    def specificActivities(self):
        lines = self.__getDataLines('SA')
        return [self.extractDataLineInfo(line, numeric_value=True) for line in lines]

    @property
    def substratesAndProducts(self) -> list:
        """
        Returns list of dicts with evaluated "natural" substrates and products
        of the enzyme across organisms.
        """
        substrates, products, res = [], [], []
        lines = self.__getDataLines('NSP')
        for line in lines:
            data = self.extractDataLineInfo(line)
            rxn = data['value'].replace(
                '{}', '').replace('?', '').replace('more', '').strip()
            try:
                subs, prods = rxn.split('=')
                subs = [s.strip() for s in subs.split('+') if s.strip() != '']
                prods = [s.strip() for s in prods.split('+') if s.strip() != '']
                subs.sort()
                prods.sort()
                if (subs not in substrates and len(subs) > 0 and len(prods) > 0):
                    substrates.append(subs)
                    products.append(prods)
                    res.append({'substrates': subs, 'products': prods})
            except Exception:
                pass
        return res

    @property
    def temperature(self):
        return EnzymeConditionDict({'optimum': self.__extractTempOrPHData('TO'),
                                    'range': self.__extractTempOrPHData('TR'),
                                    'stability': self.__extractTempOrPHData('TS')})

    @property
    def PH(self):
        return EnzymeConditionDict({'optimum': self.__extractTempOrPHData('PHO'),
                                    'range': self.__extractTempOrPHData('PHR'),
                                    'stability': self.__extractTempOrPHData('PHS')})

    @property
    def proteins(self) -> dict:
        return self.__proteins

    @property
    def organisms(self) -> list:
        """
        Returns a list containing all represented species in the database for this reaction
        """
        organisms = list(set([s['name'] for s in self.proteins.values()]))
        organisms.sort()
        return organisms

    @property
    def references(self):
        return self.__references
