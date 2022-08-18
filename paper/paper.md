---
title: 'Brendapyrser: a Python package to parse and manipulate the BRENDA database'
tags:
  - Python
  - bioinformatics
  - BRENDA
  - metabolism
  - enzymes
  - biochemistry
authors:
  - name: Semidán Robaina Estévez
    orcid: 0000-0003-0781-1677
    email: srobaina@ull.edu.es
    affiliation: 1
affiliations:
 - name: Department of Microbiology. Universidad de La Laguna.
   index: 1
date: 18 August 2022
bibliography: paper.bib
---

# Summary

Cellular metabolism, the collection of all (bio)chemical reactions that occur within the cell, is a complex system which is driven by multiple actuators, including enzymes and chemical compounds, i.e., metabolites. BRENDA [@brenda] is a database of biochemical reactions and their associated properties.

The forces on stars, galaxies, and dark matter under external gravitational
fields lead to the dynamical evolution of structures in the universe. The orbits
of these bodies are therefore key to understanding the formation, history, and
future state of galaxies. The field of "galactic dynamics," which aims to model
the gravitating components of galaxies to study their structure and evolution,
is now well-established, commonly taught, and frequently used in astronomy.
Aside from toy problems and demonstrations, the majority of problems require
efficient numerical tools, many of which require the same base code (e.g., for
performing numerical orbit integration).

# Statement of need

`Brendapyrser` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package [@astropy] (`astropy.units` and
`astropy.coordinates`).

`Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
<!-- ![Caption for example figure.\label{fig:example}](figure.png) -->
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
<!-- ![Caption for example figure.](figure.png){ width=20% } -->

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References