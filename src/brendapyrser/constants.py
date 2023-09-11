#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
BRENDA data fields and units of measurement.
"""

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