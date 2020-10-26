# A python parser for the BRENDA database

This project provides python classes and functions to parse the text file containing the entire BRENDA enzyme database (https://www.brenda-enzymes.org)

This is an ongoing project!


```python
from parseBRENDA import BRENDA
workDir = 'C:/Users/tinta/OneDrive/Documents/Projects/BRENDA'
dataFile = workDir + '/brenda_download.txt'
```


```python
# Let's load the database
brenda = BRENDA(dataFile)
brenda
```





<table>
    <tr>
        <td><strong>Number of Enzymes</strong></td><td>7558</td>
    </tr><tr>
        <td><strong>BRENDA copyright</strong></td><td>Copyrighted by Dietmar Schomburg, Techn. University
Braunschweig, GERMANY. Distributed under the License as stated
at http:/www.brenda-enzymes.org</td>
    </tr><tr>
        <td><strong>Parser version</strong></td><td>0.0.1</td>
    </tr><tr>
        <td><strong>Author</strong></td><td>Semidán Robaina Estévez, 2020</td>
    </tr>
</table>





```python
# We can retrieve an enzyme entry by its EC number like this
r = brenda.reactions.get_by_id('2.7.1.40')
r
```





<table>
    <tr>
        <td><strong>Enzyme identifier</strong></td><td>2.7.1.40</td>
    </tr><tr>
        <td><strong>Name</strong></td><td>pyruvate kinase</td>
    </tr><tr>
        <td><strong>Systematic name</strong></td><td>ATP:pyruvate 2-O-phosphotransferase</td>
    </tr><tr>
        <td><strong>Reaction type</strong></td><td>phospho group transfer</td>
    </tr><tr>
        <td><strong>Mechanism</strong></td><td>ATP + pyruvate <=> ADP + phosphoenolpyruvate (#6,23,32,49,69,70,92#mechanism <5,92>; #69# compulsory-ordered tri-bi mechanism <100>; #103#model for allosteric regulation <91>; #120,121# allosteric enzyme <98>;#16,18,36# allosteric enzyme: homotropic <73,74,75>; #12# hyperbolickinetics <96>; #118# sigmoidal kinetics with respect tophosphoenolpyruvate <93>; #52# catalyzes the addition of a proton andthe loss of a phosphoryl group which is transferred to ADP <48>; #98#sigmoidal saturation curves with substrate and metal ions <89>; #140#modeling of the catalytic mechanism, overview <248>; #149# the kineticmechanism is random order with a rapid equilibrium <255>)</td>
    </tr>
</table>





```python
# Here are all the KM values associated with this enzyme
r.KMvalues.filter_by_organism('Bos taurus').filter_by_compound('phosphoenolpyruvate').get_values()
```




    [0.051500000000000004, 0.18]




```python
# We can also get information about operational temperatures of the enzyme
temp = r.temperature
```


```python
print([l['value'] for l in temp['optimum']])
```

    [30.0, 37.0, 95.0, 22.0, 23.0, 44.0, 28.0, 40.0, 55.0, 45.0, -999.0, 60.0, 15.0, 20.5, 85.0, 50.0, 80.0, 25.0, 98.0]
    


```python
print([l['value'] for l in temp['range']][0])
```

    [35.0, 53.0]
    


```python
r.Kcatvalues.filter_by_compound('phosphoenolpyruvate').get_values()
```




    [3.2,
     161.0,
     1182.0,
     66.0,
     12.1,
     13.9,
     12.2,
     21.5,
     232.0,
     1.97,
     226.0,
     41.3,
     58.4,
     1.0,
     0.38,
     39.67,
     50.12,
     73.04,
     3204.0,
     1736.0]




```python
r.substratesAndProducts   
```




    [{'substrates': ['AKT1S1', 'ATP'], 'products': ['ADP', 'phospho-AKT1S1']},
     {'substrates': ['TDP', 'phosphoenolpyruvate'],
      'products': ['TTP', 'pyruvate | 95% yield |']},
     {'substrates': ['ATP', 'pyruvate'],
      'products': ['ADP', 'phosphoenolpyruvate']},
     {'substrates': ['ADP', 'phosphoenolpyruvate'],
      'products': ['ATP', 'pyruvate']},
     {'substrates': ['ATP', 'prothymosin alpha'],
      'products': ['ADP', 'phospho-prothymosin alpha']}]


