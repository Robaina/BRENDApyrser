# A python parser for the BRENDA database

This project provides python classes and functions to parse the text file containing the entire BRENDA enzyme database (https://www.brenda-enzymes.org)

Due to BRENDA's license, BRENDA's database cannot be downloaded directly by the parser, instead, the user is asked to download the database as a text file after accepting usage conditions [here](https://www.brenda-enzymes.org/download_brenda_without_registration.php).

This is an ongoing project!

## Installation
1. ```pip install brendapyrser```

or

2. Git clone project to local directory.

   In terminal navigate to directory and enter: ```python setup.py install```


```python
import numpy as np
from matplotlib import pyplot as plt
from brendapyrser import BRENDA

dataFile = 'data/brenda_download.txt'
```

## 1. Parsing BRENDA


```python
# Let's load the database
brenda = BRENDA(dataFile)
brenda
```





<table>
    <tr>
        <td><strong>Number of Enzymes</strong></td><td>7609</td>
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
# Plot all Km values in the database
BRENDA_KMs = np.array([v for r in brenda.reactions 
                       for v in r.KMvalues.get_values()])
values = BRENDA_KMs[(BRENDA_KMs < 1000) & (BRENDA_KMs >= 0)]
plt.hist(values)
plt.title(f'Median KM value: {np.median(values)}')
plt.xlabel('KM (mM)')
plt.show()
print(f'Minimum and maximum values in database: {values.min()} mM, {values.max()} mM')
```


    
![png](README_files/output_5_0.png)
    


    Minimum and maximum values in database: 0.0 mM, 997.0 mM



```python
# Plot all Km values in the database
BRENDA_Kcats = np.array([v for r in brenda.reactions 
                       for v in r.Kcatvalues.get_values()])
values = BRENDA_Kcats[(BRENDA_Kcats < 1000) & (BRENDA_Kcats >= 0)]
plt.hist(values)
plt.title(f'Median Kcat value: {np.median(values)}')
plt.xlabel('Kcat (1/s)')
plt.show()
print(f'Minimum and maximum values in database: {values.min()} 1/s, {values.max()} 1/s')
```


    
![png](README_files/output_6_0.png)
    


    Minimum and maximum values in database: 5.83e-10 1/s, 997.0 1/s



```python
# Plot all enzyme optimal temperature values in the database
BRENDA_TO = np.array([v for r in brenda.reactions 
                       for v in r.temperature.filter_by_condition(
                           'optimum').get_values()])
values = BRENDA_TO[(BRENDA_TO >= 0)]
plt.hist(values)
plt.title(f'Median Optimum Temperature: {np.median(values)}')
plt.xlabel('TO (${}^oC$)')
plt.show()
print(f'Minimum and maximum values in database: {values.min()} °C, {values.max()} °C')
```


    
![png](README_files/output_7_0.png)
    


    Minimum and maximum values in database: 0.0 °C, 125.0 °C


We see that the median optimal temperature for all enzymes in the BRENDA database is 37 °C! That's interesting... perhaps all organisms have agreed to prefer that temperature over other ones... or, more likely, it could be that BRENDA database is biased towards mammals and microorganisms that live within mammals... such as human pathogens.

Let's filter results for a particular species, let's try with a hyperthermophylic baterial genus, _Thermotoga_


```python
# Plot all enzyme optimal temperature values in the database
species = 'Thermotoga'
BRENDA_TO = np.array([v for r in brenda.reactions.filter_by_organism(species)
                       for v in r.temperature.filter_by_condition('optimum').filter_by_organism(species).get_values()])
values = BRENDA_TO[(BRENDA_TO >= 0)]
plt.hist(values)
plt.title(f'Median Optimum Temperature: {np.median(values)}')
plt.xlabel('TO (${}^oC$)')
plt.show()
print(f'Minimum and maximum values in database: {values.min()} °C, {values.max()} °C')
```


    
![png](README_files/output_9_0.png)
    


    Minimum and maximum values in database: 20.0 °C, 105.0 °C


We can see that the median optimal temperature among all enzymes in the genus, 80°C, is much higher than in the case of the entire database. That's consistent with the fact that _Thermotoga_ are hyperthermophylic... alright!

## 2. Extracting data for _Pyruvate kinase_


```python
# We can retrieve an enzyme entry by its EC number like this
r = brenda.reactions.get_by_id('2.7.1.40')
r
```





<table>
    <tr>
        <td><strong>Enzyme identifier</strong></td><td>2.7.1.40</td>
    </tr><tr>
        <td><strong>Name</strong></td><td>Pyruvate kinase</td>
    </tr><tr>
        <td><strong>Systematic name</strong></td><td>ATP:pyruvate 2-O-phosphotransferase</td>
    </tr><tr>
        <td><strong>Reaction type</strong></td><td>Phospho group transfer</td>
    </tr><tr>
        <td><strong>Reaction</strong></td><td>ATP + pyruvate <=> ADP + phosphoenolpyruvate</td>
    </tr>
</table>





```python
# Here are all the KM values for phosphoenolpyruvate associated with this enzyme class
compound = 'phosphoenolpyruvate'
kms = r.KMvalues.filter_by_compound(compound).get_values()
plt.hist(kms)
plt.xlabel('KM (mM)')
plt.title(f'{r.name} ({compound})')
plt.show()
```


    
![png](README_files/output_13_0.png)
    



```python
# Here are all the KM values for phosphoenolpyruvate associated with this enzyme class
compound = 'phosphoenolpyruvate'
KMs = r.KMvalues.filter_by_compound(compound).get_values()
plt.hist(KMs)
plt.xlabel('KM (mM)')
plt.title(f'{r.name} ({compound})')
plt.show()
```


    
![png](README_files/output_14_0.png)
    



```python
# And further filtered by organism
r.KMvalues.filter_by_organism('Bos taurus').filter_by_compound('phosphoenolpyruvate').get_values()
```




    [0.051500000000000004, 0.18]




```python
# Here are all the Kcat values for phosphoenolpyruvate associated with this enzyme class
compound = 'phosphoenolpyruvate'
kcats = r.Kcatvalues.filter_by_compound(compound).get_values()
plt.hist(kcats)
plt.xlabel('Kcat ($s^{-1}$)')
plt.title(f'{r.name} ({compound})')
plt.show()
```


    
![png](README_files/output_16_0.png)
    


## 3 Finding all KM values for a given substrate and organism
Next, we will retrieve KM values associated to a particular substrate for all enzymes in a given species. Will t he KM values distribute around a narrow or wider concentration range? Since substrate concentration in cytoplasma is the same for all enzymes it makes sense that all cytoplasmi enzymes utilizing that substrate have similar KM values. Let's test this idea with _Escherichia coli_ and some common substrates participating in the central carbon metabolism.


```python
species, compound = 'Escherichia coli', 'NADH'
KMs = np.array([v for r in brenda.reactions.filter_by_organism(species)
                for v in r.KMvalues.filter_by_compound(compound).filter_by_organism(species).get_values()])

if len(KMs) > 0:
    plt.hist(KMs)
    plt.xlabel('KM (mM)')
    plt.title(f'{species} KMs ({compound}), median = {np.median((KMs))}')
    plt.show()
else:
    print('No KM values for compound')
```


    
![png](README_files/output_18_0.png)
    


That's interesting! typical NADH concentrations are low in _Escherichia coli_, e.g., from [BioNumbers](http://book.bionumbers.org/what-are-the-concentrations-of-free-metabolites-in-cells/) we get a value of 0.083 mM. The median KM value for NADH among all enzymes binding it is lower as we see in the plot above! Hence, it looks like most enzymes are (nearly) saturated for NADH and thus fluxes are sort of independent of NADH concentration.

# 4 Filtering reactions by specific compound

We can also filter reactions in BRENDA by a specific compound: substrate, product or either of the two. Let's filter reactions containg _geraniol_ as a substrate, product or both to exemplify this feature


```python
substrate_rxns = brenda.reactions.filter_by_substrate("phosphoenolpyruvate")
substrate_rxns[2]
```





<table>
    <tr>
        <td><strong>Enzyme identifier</strong></td><td>2.5.1.19</td>
    </tr><tr>
        <td><strong>Name</strong></td><td>3-phosphoshikimate 1-carboxyvinyltransferase</td>
    </tr><tr>
        <td><strong>Systematic name</strong></td><td>phosphoenolpyruvate:3-phosphoshikimate 5-O-(1-carboxyvinyl)-transferase</td>
    </tr><tr>
        <td><strong>Reaction type</strong></td><td>Enolpyruvate group transfer (#3,52,55# induced-fit mechanism, formation</td>
    </tr><tr>
        <td><strong>Reaction</strong></td><td>phosphoenolpyruvate + 3-phosphoshikimate <=> phosphate +5-O-</td>
    </tr>
</table>





```python
compound_rxns = brenda.reactions.filter_by_compound("phosphoenolpyruvate")
compound_rxns[7]
```





<table>
    <tr>
        <td><strong>Enzyme identifier</strong></td><td>2.5.1.7</td>
    </tr><tr>
        <td><strong>Name</strong></td><td>Udp-n-acetylglucosamine 1-carboxyvinyltransferase</td>
    </tr><tr>
        <td><strong>Systematic name</strong></td><td>phosphoenolpyruvate:UDP-N-acetyl-D-glucosamine</td>
    </tr><tr>
        <td><strong>Reaction type</strong></td><td>Carboxyvinyl group transfer</td>
    </tr><tr>
        <td><strong>Reaction</strong></td><td>phosphoenolpyruvate + UDP-N-acetyl-alpha-D-glucosamine <=> phosphate +UDP-N-acetyl-3-O-</td>
    </tr>
</table>



