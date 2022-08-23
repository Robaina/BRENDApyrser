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
 - name: Department of Microbiology. University of La Laguna. Spain.
   index: 1
date: 18 August 2022
bibliography: paper.bib
---

# Summary

Enzymes &mdash; proteins with specialized, catalytic functions &mdash; constitute the workforce of cellular metabolism, catalyzing thousands of biochemical reactions within cells. Enzymes present different physicochemical properties that affect their function. Knowledge about these enzyme functional properties is fundamental to understanding how biochemical reactions operate and are controlled by the cell. The BRENDA [@brenda] database is a widely-used, publicly available collection of enzyme functional information obtained from the primary literature. The development of computational tools to parse and query BRENDA would facilitate its integration in analyses of cellular metabolism.

# Statement of need

Users can access the BRENDA database directly on the website [@brenda-web], which provides searching and filtering capabilities. BRENDA also offers an API to access the database programmatically within several programming languages. However, obtaining specific data through the browser or the API turns inefficient in certain applications, for instance, when extracting certain data fields for thousands of enzymes to conduct statistical analyses. Here, we present `Brendapyrser`, a Python package to parse and manipulate the BRENDA database. Instead of accessing BRENDA via its API, `Brendapyrser` provides a collection of objects and methods to parse BRENDA locally as a text file &mdash; currently sized under 300 MB &mdash;, thus extracting data fields more quickly. Moreover, its syntax and object-oriented organization are well-suited for exploratory analyses within the Python ecosystem.

# Acknowledgements

We acknowledge constructive feedback from Brendapyrser users that has helped improve the package. This work has been conducted without any financial or commercial support.

# References