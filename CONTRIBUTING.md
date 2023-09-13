# Contributing to BRENDApyrser

First of all, thanks for taking the time to contribute! :tada::+1:

Here you will find a set of guidelines for contributing to BRENDApyrser. Feel free to propose changes to this document in a pull request.

## Code of conduct

This project and everyone participating in it is governed by the [Contributor Covenant, v2.0](.github/CODE_OF_CONDUCT.md) code of conduct. By participating, you are expected to uphold this code.

## I have a question!

If you only have a question about all things related to BRENDApyrser, the best course of actions for you is to open a new [discussion](https://github.com/Robaina/BRENDApyrser/discussions).

## How can I contribute?

### 1. Reporting bugs

We all make mistakes, and the developers behind BRENDApyrser are no exception... So, if you find a bug in the source code, please open an [issue](https://github.com/Robaina/BRENDApyrser/issues) and report it. Please, first search for similar issues that are currrently open.

### 2. Suggesting enhancements

Are you missing some feature that would like BRENDApyrser to have? No problem! You can contribute by suggesting an enhancement, just open a new issue and tag it with the [```enhancement```](https://github.com/Robaina/BRENDApyrser/labels/enhancement) label. Please, first search for similar issues that are currrently open.

### 3. Improving the documentation

Help is always needed at improving the [documentation](https://robaina.github.io/BRENDApyrser/). Either adding more detailed docstrings, usage explanations or new examples.

## First contribution

Unsure where to begin contributing to BRENDApyrser? You can start by looking for issues with the label [```good first issue```](https://github.com/Robaina/BRENDApyrser/labels/good%20first%20issue). If you are unsure about how to set a developer environment for BRENDApyrser, do take a look at the section below. Thanks!

## Setting up a local developer environment

To setup up a developer environment for BRENDApyrser:

1. Fork and download repo, cd to downloaded directory. You should create a new branch to work on your issue.

2. Create conda environment with required dependencies:

The file `envs/BRENDApyrser-dev.yml` contains all dependencies required to use BRENDApyrser. Conda is very slow solving the environment. It is recommended to use [mamba](https://github.com/mamba-org/mamba) instead:

```bash
mamba env create -n BRENDApyrser-dev -f envs/BRENDApyrser-dev.yml
conda activate BRENDApyrser-dev
```

3. Build package

```bash
(BRENDApyrser-dev) poetry build
```

4. Install BRENDApyrser

```bash
(BRENDApyrser-dev) pip install dist/BRENDApyrser*.whl
```

5. Run tests

```bash
(BRENDApyrser-dev) python -m unittest discover tests
```

## Building the documentation

The documentation is formed by a series of markdown files located in directory [docs](https://github.com/Robaina/BRENDApyrser/tree/main/docs). This repo uses [mkdocs](https://www.mkdocs.org/) to automatically generate documentation pages from markdown files. Also, [MathJax](https://github.com/mathjax/MathJax) syntax is allowed!

This means that, to modify the [API reference](https://robaina.github.io/BRENDApyrser/references/api/), all you need to do is to modify the docstring directly in the source file where the definion/class is located. And, to update the documentation pages, you just have to update the corresponding markdown file in the [docs](https://github.com/Robaina/BRENDApyrser/tree/main/docs) directory. Note that, if you need to change the documentation structure (e.g., add or new pages),you would need to tell mkdocs about this change through its [configuration file](https://github.com/Robaina/BRENDApyrser/blob/main/mkdocs.yml). Or just open an issue and ask for help!

When all the changes are ready to deploy, just open a pull request. After reviewing and merging the changes, the documentation will be automatically deployed.

Run the documentation locally with:

> mkdocs serve

## Tests on push and pull request to main

BRENDApyrser's repo contains a [GitHub Action](https://github.com/features/actions) to perform build and integration tests which is triggered automatically on push and pull request events to the main brach. Currently the tests include building and installing BRENDApyrser in Ubuntu and MacOS and running the [test](tests) suit.
