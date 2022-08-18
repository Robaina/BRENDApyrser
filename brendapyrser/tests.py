#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unit tests for brendapyrser.
"""

import unittest
from brendapyrser import Reaction, ReactionList


rxn_data = """ID	1.1.1.304
********************************************************************************
*                                                                              *
* Copyrighted by Dietmar Schomburg, Techn. University Braunschweig, GERMANY    *
* Distributed under the License as stated at http:/www.brenda-enzymes.org      *
*                                                                              *
********************************************************************************

PROTEIN
PR	#1# Staphylococcus aureus   <4>
PR	#2# Geobacillus stearothermophilus   <5,8>
PR	#3# Klebsiella aerogenes   <1>
PR	#4# Bacillus licheniformis   <10>
PR	#5# Rhodococcus erythropolis   <7>
PR	#6# Columba livia   <2,3>
PR	#7# Klebsiella pneumoniae Q48436 UniProt <6>
PR	#8# Paenibacillus polymyxa KC505218 GenBank <9>
PR	#9# Mycobacterium sp. W8VSK8 UniProt <11>

RECOMMENDED_NAME
RN	diacetyl reductase [(S)-acetoin forming]


SYSTEMATIC_NAME
SN	(S)-acetoin:NAD+ oxidoreductase


SYNONYMS
SY	#2# S-stereospecific diacetyl reductase <8>
SY	#2# BSDR <8>
SY	#4# budC <10>
SY	#5# AdR <7>
SY	#5# acetoin(diacetyl) reductase <7>
SY	#5# ReADR <7>
SY	#8# diacetyl reductase <9>
SY	#8# DAR <9>
SY	#9# ADS1 <11>

REACTION
RE	(S)-acetoin + NAD+ = diacetyl + NADH + H+ (#6# Theorell-Chance
	mechanism with NADH as the leading substrate <3>)
RE	(S)-acetoin + NAD+ = diacetyl + NADH + H+ (#6# Theorell-Chance
	mechanism with NADH as the leading substrate <3>)

SOURCE_TISSUE
ST	#6# liver <2,3>

NATURAL_SUBSTRATE_PRODUCT
NSP	#2,5# diacetyl + NADH + H+ = (S)-acetoin + NAD+ <7,8>
NSP	#2,5# diacetyl + NADH + H+ = (S)-acetoin + NAD+ {ir} <7,8>
NSP	#2,5# diacetyl + NADH + H+ = (S)-acetoin + NAD+ {r} <7,8>
NSP	#2,5,8# more = ? (#5# Rhodococcus erythropolis WZ010 is capable of
	producing optically pure (2S,3S)-2,3-butanediol in alcoholic
	fermentation.  <7>; #8# the enzyme also catalyzes the stereospcific
	reaction of (S)-acetoin reduction to butanediol, EC 1.1.1.76 <9>; #2#
	the enzyme shows an S-enantioselectivity in the reversible reduction of
	acetoin so it might be responsible of the meso-butenediol formation
	from R-acetoin. It acts on racemic acetoin and (S)-acetoin to form
	(2S,3S)-butane-2,3-diol, EC 1.1.1.76, but also on the
	(2R,3R)-butane-2,3-diol isomer in the reverse reaction, EC 1.1.1.4 <8>)
	<7,8,9>

SUBSTRATE_PRODUCT
SP	#1# 2,3-pentanedione + beta-NADH + H+ = L-3-hydroxy-2-pentanone +
	beta-NAD+ {ir} <4>
SP	#1# diacetyl + beta-NADH + H+ = (S)-acetoin + beta-NAD+ (#1# 86.9% of
	the activity with pentane-2,3-dione <4>) {ir} <4>
SP	#1# ethyl pyruvate + beta-NADH + H+ = ? + beta-NAD+ (#1# 38.4% of the
	activity with pentane-2,3-dione <4>) {ir} <4>
SP	#1# methyl pyruvate + beta-NADH + H+ = ? + beta-NAD+ (#1# 22.8% of the
	activity with pentane-2,3-dione <4>) {ir} <4>
SP	#1,2,4,5,8,9# more = ? (#1# no activity with alpha-NADH or NADPH <4>;
	#5# Rhodococcus erythropolis WZ010 is capable of producing optically
	pure (2S,3S)-2,3-butanediol in alcoholic fermentation.  <7>; #8# the
	enzyme also catalyzes the stereospcific reaction of (S)-acetoin
	reduction to butanediol, EC 1.1.1.76 <9>; #2# the enzyme shows an
	S-enantioselectivity in the reversible reduction of acetoin so it might
	be responsible of the meso-butenediol formation from R-acetoin. It acts
	on racemic acetoin and (S)-acetoin to form (2S,3S)-butane-2,3-diol, EC
	1.1.1.76, but also on the (2R,3R)-butane-2,3-diol isomer in the reverse
	reaction, EC 1.1.1.4 <8>; #5# the enzyme displays absolute
	stereospecificity in the reduction of diacetyl to
	(2S,3S)-2,3-butanediol via (S)-acetoin. The enzyme shows higher
	catalytic efficiency for (S)-1-phenylethanol oxidation than that for
	acetophenone reduction. ReADR-catalyzed asymmetric reduction of
	diacetyl is coupled with stereoselective oxidation of 1-phenylethanol,
	which simultaneously forms both (2S,3S)-2,3-butanediol and
	(R)-1-phenylethanol in great conversions and enantiomeric excess
	values.The enzyme accepts a broad range of substrates including
	aliphatic and aryl alcohols, aldehydes, and ketones <7>; #9# enzyme
	shows activity as a reductase specific for (S)-acetoin, EC 1.1.1.76,
	and both diacetyl reductase (EC 1.1.1.304) and NAD+-dependent alcohol
	dehydrogenase (EC 1.1.1.1) activities <11>; #4# enzyme shows oxidative
	activity to racemic 2,3-butanediol but no activity toward racemic
	acetoin in the presence of NAD+ <10>) <4,7,8,9,10,11>
SP	#2,3,5,6,7,9# diacetyl + NADH + H+ = (S)-acetoin + NAD+ (#7# 87% of the
	(R)-2,3-butanediol dehydrogenase activity with substrate acetoin <6>)
	<1,2,3,5,6,7,8,11>
SP	#2,3,5,6,7,9# diacetyl + NADH + H+ = (S)-acetoin + NAD+ (#7# 87% of the
	(R)-2,3-butanediol dehydrogenase activity with substrate acetoin <6>)
	{ir} <1,2,3,5,6,7,8,11>
SP	#2,3,5,6,7,9# diacetyl + NADH + H+ = (S)-acetoin + NAD+ (#7# 87% of the
	(R)-2,3-butanediol dehydrogenase activity with substrate acetoin <6>)
	{r} <1,2,3,5,6,7,8,11>
SP	#3# ethyl pyruvate + NADH + H+ = ? + NAD+ (#3# 57.7% of the activity
	with diacetyl <1>) <1>
SP	#3# methyl glyoxal + NADH + H+ = ? + NAD+ (#3# 11% of the activity with
	diacetyl <1>) <1>
SP	#3# methyl pyruvate + NADH + H+ = ? + NAD+ (#3# 49% of the activity
	with diacetyl <1>) <1>
SP	#3,7# 2,3-pentanedione + NADH + H+ = 3-hydroxy-2-pentanone + NAD+ (#7#
	77% of the (R)-2,3-butanediol dehydrogenase activity with substrate
	acetoin <6>; #3# 85.6% of the activity with diacetyl <1>) <1,6>
SP	#4# (2S,3S)-2,3-butanediol + NAD+ = (3S)-acetoin + NADH + H+ {r} <10>
SP	#4# (3S)-acetoin + NADH + H+ = (2S,3S)-2,3-butanediol + NAD+ (#4# 97%
	of the activity with diacetyl <10>) |#4# 97.3% enantiomeric excess and
	96.5% diastereomeric excess <10>| {r} <10>
SP	#4# 1,2-propanediol + NAD+ = ? + NADH + H+ |#4# 0.5% of the activity
	with 2,3-butanediol <10>| {r} <10>
SP	#4# diacetyl + NADH + H+ = (3S)-acetoin + NAD+ |#4# 97.3% enantiomeric
	excess <10>| {ir} <10>
SP	#4# 2,3-pentanedione + NADH + H+ = ? + NAD+ |#4# 69% of the activity
	with diacetyl <10>| <10>
SP	#4# 2,3-hexanedione + NADH + H+ = ? + NAD+ |#4# 66% of the activity
	with diacetyl <10>| <10>
SP	#4# 3,4-hexanedione + NADH + H+ = ? + NAD+ |#4# 10% of the activity
	with diacetyl <10>| <10>
SP	#6,8# diacetyl + NADPH + H+ = (S)-acetoin + NADP+ {ir} <3,9>

TURNOVER_NUMBER
TN	#4# 748 {NAD+}  (#4# pH 10.0, 30°C <10>) <10>
TN	#4# 202 {(3S)-acetoin}  (#4# pH 6.0, 30°C <10>) <10>
TN	#4# 591 {(2S,3S)-2,3-butanediol}  (#4# pH 10.0, 30°C <10>) <10>
TN	#4# 1222 {diacetyl}  (#4# pH 6.0, 30°C <10>) <10>
TN	#4# 1274 {NADH}  (#4# pH 6.0, 30°C <10>) <10>
TN	#9# 110 {NADH}  (#9# pH 7.0, 30°C <11>) <11>
TN	#9# 163 {diacetyl}  (#9# pH 7.0, 30°C <11>) <11>

KM_VALUE
KM	#1# 0.045 {NADH}  (#1# 25°C, pH 6.0, cosubstrate diacetyl <4>) <4>
KM	#1# 0.095 {NADH}  (#1# 25°C, pH 6.0, cosubstrate methyl pyruvate <4>)
	<4>
KM	#1# 0.025 {NADH}  (#1# 25°C, pH 6.0, cosubstrate 2,3-pentanedione <4>)
	<4>
KM	#1# 0.11 {NADH}  (#1# 25°C, pH 6.0, cosubstrate ethyl pyruvate <4>) <4>
KM	#1# 24 {ethyl pyruvate}  (#1# 25°C, pH 6.0 <4>) <4>
KM	#1# 16 {Methyl pyruvate}  (#1# 25°C, pH 6.0 <4>) <4>
KM	#1# 6 {2,3-Pentanedione}  (#1# 25°C, pH 6.0 <4>) <4>
KM	#1# 15 {diacetyl}  (#1# 25°C, pH 6.0 <4>) <4>
KM	#2# 19 {diacetyl}  (#2# pH 7.5, 25°C <5>) <5>
KM	#3# 0.005 {NADH}  (#3# cosubstrate acetoin, pH 7.0, 25°C <1>) <1>
KM	#3# 20 {ethyl pyruvate}  (#3# pH 7.0, 25°C <1>) <1>
KM	#3# 0.007 {NADH}  (#3# cosubstrate diacetyl, pH 7.0, 25°C <1>) <1>
KM	#3# 1.6 {diacetyl}  (#3# pH 7.0, 25°C <1>) <1>
KM	#3# 6 {pentane-2,3-dione}  (#3# pH 7.0, 25°C <1>) <1>
KM	#3# 18 {Methyl pyruvate}  (#3# pH 7.0, 25°C <1>) <1>
KM	#3# 75 {methyl glyoxal}  (#3# pH 7.0, 25°C <1>) <1>
KM	#4# 0.25 {NADH}  (#4# pH 6.0, 30°C <10>) <10>
KM	#4# 0.34 {NAD+}  (#4# pH 10.0, 30°C <10>) <10>
KM	#4# 0.47 {(3S)-acetoin}  (#4# pH 6.0, 30°C <10>) <10>
KM	#4# 7.25 {(2S,3S)-2,3-butanediol}  (#4# pH 10.0, 30°C <10>) <10>
KM	#4# 72.4 {diacetyl}  (#4# pH 6.0, 30°C <10>) <10>
KM	#5# -999 {more}  (#5# Michaelis-Menten-type kinetics <7>) <7>
KM	#6# 0.1 {NADH}  (#6# pH 6.1, 25°C <3>) <3>
KM	#6# 3.1 {diacetyl}  (#6# pH 6.1, 25°C <3>) <3>
KM	#6# 3 {diacetyl}  (#6# 25°C, pH 6.1 <2>) <2>
KM	#6# 0.087 {NADH}  (#6# 25°C, pH 5.9 <2>) <2>
KM	#6# 0.116 {NADH}  (#6# 25°C, pH 6.1 <2>) <2>
KM	#6# 2.64 {diacetyl}  (#6# 25°C, pH 6.7 <2>) <2>
KM	#6# 0.135 {NADH}  (#6# 25°C, pH 6.7 <2>) <2>
KM	#6# 2.81 {diacetyl}  (#6# 25°C, pH 5.9 <2>) <2>
KM	#9# 0.05 {NADH}  (#9# pH 7.0, 30°C <11>) <11>
KM	#9# 4.47 {diacetyl}  (#9# pH 7.0, 30°C <11>) <11>

PH_OPTIMUM
PHO	#1,8# 6 (#8# assay at <9>) <4,9>
PHO	#2# 6.5 (#2# assay at <8>) <8>
PHO	#4# 5 (#4# reduction of diacetyl <10>) <10>
PHO	#4# 10 (#4# oxidation of butanediol <10>) <10>
PHO	#5# 7 (#5# diacetyl reduction <7>) <7>
PHO	#6# 6.1 <3>

PH_RANGE
PHR	#4# 5-8 (#4# reduction of diacetyl <10>) <10>
PHR	#6# 5 (#6# 5 min, 30% loss of activity <2>) <2>
PHR	#6# 5.1 (#6# 5 min, 20% loss of activity <2>) <2>
PHR	#6# 5.4-7.6 (#6# stable within <2>) <2>
PHR	#6# 4.8 (#6# 5 min, 60% loss of activity <2>) <2>

SPECIFIC_ACTIVITY
SA	#2# 71.4 (#2# pH 7.5, 25°C <5>) <5>
SA	#4# 120.0 (#4# substrate diacetyl, 30°C, pH 6.0 <10>) <10>
SA	#8# 72.6 (#8# purified recombinant enzyme, NADPH, pH 6.0, 30°C <9>) <9>

TEMPERATURE_OPTIMUM
TO	#2# 50 <5>
TO	#5,8# 30 (#8# assay at <9>; #5# diacetyl reduction <7>) <7,9>

COFACTOR
CF	#1,2,3,4,5,6# NADH (#2# dependent on <8>; #1# beta-NADH <4>; #3#
	specific for beta-NADH <1>) <1,3,4,5,7,8,10>
CF	#1,4,8# more (#4# no activity with NADPH <10>; #8# inactive with NADH
	<9>; #1# no activity with alpha-NADH or NADPH <4>) <4,9,10>
CF	#4,5# NAD+ <7,10>
CF	#8# NADPH (#8# dependent on <9>) <9>

ACTIVATING_COMPOUND
AC	#5# DMSO (#5# DMSO at a final concentration of 30% v/v added into the
	assay mixture, increases the activity up to 120% of the control enzyme
	activity <7>) <7>

INHIBITORS
IN	#1# diacetyl (#1# substrate inhibition at concentrations above 80-90 mM
	<4>) <4>
IN	#1# ethyl pyruvate (#1# substrate inhibition at concentrations above
	80-90 mM <4>) <4>
IN	#1# Methyl pyruvate (#1# substrate inhibition at concentrations above
	80-90 mM <4>) <4>
IN	#4# Cu2+ (#4# 1 mM, 0.2% of initial activity with substrate diacetyl,
	1% with substrate 2,3-butanediol, respectively <10>) <10>
IN	#4# Ag+ (#4# 1 mM, 0.5% of initial activity with substrate diacetyl,
	0.5% with substrate 2,3-butanediol, respectively <10>) <10>
IN	#4# Fe3+ (#4# 1 mM, 2% of initial activity with substrate diacetyl,
	3.5% with substrate 2,3-butanediol, respectively <10>) <10>
IN	#4# EDTA (#4# 1 mM, 91% of initial activity with substrate diacetyl,
	90% with substrate 2,3-butanediol, respectively <10>) <10>
IN	#4,5# Al3+ (#4# 1 mM, 4% of initial activity with substrate diacetyl,
	6% with substrate 2,3-butanediol, respectively <10>) <7,10>
IN	#4,5# Zn2+ (#5# inhibits 94.6% at 2 mM <7>; #4# 1 mM, 75% of initial
	activity with substrate diacetyl, 80% with substrate 2,3-butanediol,
	respectively <10>) <7,10>
IN	#5# Fe2+ (#5# inhibits 91.6% at 2 mM <7>) <7>
IN	#6# 2-oxoglutarate (#6# noncompetitive <3>) <3>
IN	#6# hexane-2,5-dione (#6# noncompetitive <3>) <3>
IN	#6# NAD+ (#6# competitive, product inhibition <3>) <3>
IN	#6# acetoin (#6# noncompetitive, product inhibition <3>) <3>
IN	#6# acetone (#6# competitive for diacetyl, uncompetitive for NADH <3>)
	<3>
IN	#6# Pentane-3-one (#6# competitive for diacetyl, uncompetitive for NADH
	<3>) <3>

METALS_IONS
ME	#5# more (#5# addition of EDTA or the cations at 1 mM, such as Na+, K+,
	Mn2+, Mg2+, and Ca2+, have no significant effect on the activity of
	ReADR <7>) <7>
ME	#5# Mn2+ (#5# activates by 201.6 to 265.6% at 2 mM <7>) <7>
ME	#5# K+ (#5# activates by 201.6 to 265.6% at 2 mM <7>) <7>
ME	#5# Na+ (#5# activates by 201.6 to 265.6% at 2 mM <7>) <7>

MOLECULAR_WEIGHT
MW	#1# 68000 (#1# gel filtration <4>) <4>
MW	#2# 26000 (#2# 2 * 26000, SDS-PAGE <5>) <5>
MW	#2# 49000 (#2# gel filtration <5>) <5>
MW	#3# 61000 (#3# gel filtration <1>) <1>
MW	#3# 28000 (#3# 2 * 28000, SDS-PAGE <1>) <1>
MW	#4# 125000 (#4# gel filtration <10>) <10>
MW	#4# 30000 (#4# 4 * 30000, SDS-PAGE <10>) <10>
MW	#5# 26864 (#5# 2 * 26864, sequence calculation <7>) <7>
MW	#7# 96000 (#7# gel filtration <6>) <6>
MW	#7# 26591 (#7# 4 * 26591, calculated <6>) <6>
MW	#8# 118000 (#8# recombinant enzyme, gel filtration <9>) <9>
MW	#8# 28500 (#8# 4 * 28500, recombinant enzyme, SDS-PAGE <9>) <9>
MW	#9# 150000 (#9# gel filtration <11>) <11>
MW	#9# 36000 <11>

SUBUNITS
SU	#1# monomer (#1# 1 * 68000, SDS-PAGE <4>) <4>
SU	#2,3# dimer (#3# 2 * 28000, SDS-PAGE <1>; #2# 2 * 26000, SDS-PAGE <5>)
	<1,5>
SU	#4,7,8,9# tetramer (#4# 4 * 30000, SDS-PAGE <10>; #7# 4 * 26591,
	calculated <6>; #8# 4 * 28500, recombinant enzyme, SDS-PAGE <9>; #9# 4
	* 36000, SDS-PAGE, 4 * 35971, calculated <11>) <6,9,10,11>
SU	#5# homodimer (#5# 2 * 26864, sequence calculation <7>) <7>

PI_VALUE
PI	#3# 6.8 (#3# isoelectric focusing <1>) <1>
PI	#7# 5.9-7.2 (#7# isoelectric focusing <6>) <6>
PI	#9# 4.8 (#9# isoelectric focusing <11>) <11>
PI	#9# 5.1 (#9# calculated <11>) <11>

APPLICATION
AP	#5,8# synthesis (#5# acetoin(diacetyl) reductase, i.e. 2,3-butanediol
	dehydrogenase, is one of the key enzymes in the microbial production of
	2,3-butanediol, a platform with extensive industrial applications in
	the production of plastics, printing inks, perfumes, fumigants,
	spandex, moistening and softening agents, plasticizers, and
	pharmaceutical carrier <7>; #8# the enzyme is used for production of
	S-acetoin with higher than 99.9% optical purity from diacetyl using
	whole cells of engineered Escherichia coli <9>) <7,9>

CLONED
CL	#4# (expression in Escherichia coli) <10>
CL	#5# (gene adr, DNA and amino acid sequence determination and analysis,
	sequence comparison, expression of His-tagged enzyme in Escherichia
	coli) <7>
CL	#7# (expression in Escherichia coli) <6>
CL	#8# (gene dar, fucntional expression in Escherichia coli strain Rosetta
	(DE3), resulting in production of S-acetoin with higher than 99.9%
	optical purity from diacetyl) <9>

PURIFICATION
PU	#1# <4>
PU	#2# <5>
PU	#2# (native enzyme by adsorption on diethylaminoethylSepharose and
	hydrophobic interaction chromatography) <8>
PU	#3# <1>
PU	#5# (recombinant His-tagged enzyme from Escherichia coli by nickel
	affinity chromatography) <7>
PU	#8# (recombinant enzyme from Escherichia coli strain Rosetta (DE3)) <9>

GENERAL_STABILITY
GS	#2# (unstable to dilution, kept diluted at 0°C for ca. 60 min it will
	lose 62% of activity. This inactivation is almost completely reversed
	by the addition of NAD+) <5>

ORGANIC_SOLVENT_STABILITY
OSS	#5# DMSO (#5# the enzyme retains 53.6% of the initial activity after 4
	h incubation with 30% v/v DMSO at 4°C <7>) <7>

PH_STABILITY
PHS	#7# 7-8 <6>

STORAGE_STABILITY
SS	#2# (storage at 0°C in the presence of 20% glycerol, 0.1 mM EDTA, 5 mM
	2-mercaptoethanol and 0.6 mM NAD+ in TEA buffer, pH 7.5, half-life of
	one month) <5>

REFERENCE
RF	<1> Carballo, J.; Martin, R.; Bernardo, A.; Gonzalez, J.: Purification,
	characterization and some properties of diacetyl(acetoin) reductase
	from Enterobacter aerogenes. Eur. J. Biochem. (1991) 198, 327-332.
	{Pubmed:2040298} (c)
RF	<2> Martin, R.; Diez, V.; Burgos, J.: Pigeon liver diacetyl reductase.
	Effects of pH on the kinetic parameters of the reaction. Biochim.
	Biophys. Acta (1976) 429, 293-300. {Pubmed:4124}
RF	<3> Burgos, J.; Martin, R.; Diez, V.: Pigeon liver diacetyl reductase.
	Kinetic and thermodynamic studies with NADH as coenzyme. Biochim.
	Biophys. Acta (1974) 364, 9-16. {Pubmed:4373071}
RF	<4> Vidal, I.; Gonzalez, J.; Bernardo, A.; Martin, R.: Purification and
	classification of diacetyl-reducing enzymes from Staphylococcus aureus.
	Biochem. J. (1988) 251, 461-466. {Pubmed:3041963}
RF	<5> Giovannini, P.P.; Medici, A; Bergamini, C.M.; Rippa, M.: Properties
	of diacetyl (acetoin) reductase from Bacillus stearothermophilus.
	Bioorg. Med. Chem. (1996) 4, 1197-1201. {Pubmed:8879540}
RF	<6> Ui, S.; Okajima, Y.; Mimura, A.; Kanai, H.; kobayashi, T.; Kudo,
	T.: Sequence analysis of the gene for and characterization of D-acetoin
	forming meso-2,3-butanediol dehydrogenase of Klebsiella pneumoniae
	expressed in Escherichia coli. J. Ferment. Bioeng. (1997) 83, 32-37.
	{Pubmed:}
RF	<7> Wang, Z.; Song, Q.; Yu, M.; Wang, Y.; Xiong, B.; Zhang, Y.; Zheng,
	J.; Ying, X.: Characterization of a stereospecific acetoin(diacetyl)
	reductase from Rhodococcus erythropolis WZ010 and its application for
	the synthesis of (2S,3S)-2,3-butanediol. Appl. Microbiol. Biotechnol.
	(2014) 98, 641-650. {Pubmed:23568047}
RF	<8> Giovannini, P.; Mantovani, M.; Grandini, A.; Medici, A.; Pedrini,
	P.: New acetoin reductases from Bacillus stearothermophilus: meso- and
	2R,3R-butanediol as fermentation products. J. Mol. Catal. B (2011) 69,
	15-20. {Pubmed:}
RF	<9> Gao, J.; Xu, Y.; Li, F.; Ding, G.: Production of S-acetoin from
	diacetyl by Escherichia coli transformant cells that express the
	diacetyl reductase gene of Paenibacillus polymyxa ZJ-9. Lett. Appl.
	Microbiol. (2013) 57, 274-281. {Pubmed:23701367}
RF	<10> Xu, G.C.; Bian, Y.Q.; Han, R.Z.; Dong, J.J.; Ni, Y.: Cloning,
	expression, and characterization of budC gene encoding
	meso-2,3-butanediol dehydrogenase from Bacillus licheniformis. Appl.
	Biochem. Biotechnol. (2016) 178, 604-617. {Pubmed:26494135}
RF	<11> Takeda, M.; Anamizu, S.; Motomatsu, S.; Chen, X.; Thapa Chhetri,
	R.: Identification and characterization of a mycobacterial
	NAD+-dependent alcohol dehydrogenase with superior reduction of
	diacetyl to (S)-acetoin. Biosci. Biotechnol. Biochem. (2014) 78,
	1879-1886. {Pubmed:25082080}

KI_VALUE
KI	#1# 300 {diacetyl}  (#1# 25°C, pH 6.0 <4>) <4>
KI	#1# 150 {ethyl pyruvate}  (#1# 25°C, pH 6.0 <4>) <4>
KI	#1# 150 {Methyl pyruvate}  (#1# 25°C, pH 6.0 <4>) <4>

KCAT_KM_VALUE
KKM	#4# 16.9 {diacetyl}  (#4# pH 6.0, 30°C <10>) <10>
KKM	#4# 81.5 {(2S,3S)-2,3-butanediol}  (#4# pH 10.0, 30°C <10>) <10>
KKM	#4# 432 {(3S)-acetoin}  (#4# pH 6.0, 30°C <10>) <10>
KKM	#4# 2192 {NAD+}  (#4# pH 10.0, 30°C <10>) <10>
KKM	#4# 5072 {NADH}  (#4# pH 6.0, 30°C <10>) <10>
KKM	#5# 8.519 {NADH}  (#5# pH 7.0, 30°C, recombinant enzyme <7>) <7>
KKM	#9# 210 {NADH}  (#9# pH 7.0, 30°C <11>) <11>
KKM	#9# 36.4 {diacetyl}  (#9# pH 7.0, 30°C <11>) <11>

GENERAL_INFORMATION
GI	#2,5# metabolism (#5# acetoin(diacetyl) reductase, also known as
	2,3-butanediol dehydrogenase, is one of the key enzymes in the
	microbial production of 2,3-butanediol <7>; #2# the enzyme is involved
	in the butanediol cycle, overview <8>) <7,8>
GI	#5# evolution (#5# the enzyme belongs to the family of the short-chain
	dehydrogenase/reductases <7>) <7>"""



class TestReaction(unittest.TestCase):
    def test_ec_number(self):
        rxn = Reaction(rxn_data)
        self.assertEqual(
            rxn.ec_number, "1.1.1.304",
            "Failed to correctly retrieve EC number"
            )
    def test_name(self):
        rxn = Reaction(rxn_data)
        self.assertEqual(
            rxn.name, "Diacetyl reductase [(s)-acetoin forming]",
            "Failed to correctly retrieve reaction name"
            )
    def test_sysname(self):
        rxn = Reaction(rxn_data)
        self.assertEqual(
            rxn.systematic_name, "(S)-acetoin:NAD+ oxidoreductase",
            "Failed to correctly retrieve systematic reaction name"
            )
    def test_KMvalues(self):
        rxn = Reaction(rxn_data)
        self.assertEqual(
            rxn.KMvalues.get_values()[:4], [0.045, 0.095, 0.025, 0.11],
            "Failed to correctly retrieve KM values"
            )
    def test_KKMvalues(self):
        rxn = Reaction(rxn_data)
        self.assertEqual(
            rxn.KKMvalues.get_values()[:4], [16.9, 36.4, 81.5, 432.0],
            "Failed to correctly retrieve KKM values"
            )
    def test_Kcatvalues(self):
        rxn = Reaction(rxn_data)
        self.assertEqual(
            rxn.Kcatvalues.get_values()[:4], [748.0, 202.0, 591.0, 1222.0],
            "Failed to correctly retrieve Kcat values"
            )
    def test_temperature(self):
        rxn = Reaction(rxn_data)
        self.assertEqual(
            rxn.temperature["optimum"][0]["value"], 50.0,
            "Failed to correctly retrieve temperature values"
            )

 


if __name__ == '__main__':
    unittest.main()


