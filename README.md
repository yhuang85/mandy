# Mandy

This repo contains the source code to build this [Github-page](https://yhuang85.github.io/mandy/) where I write about stuff that I find interesting.

## Build tools
The Github page is served from the **docs** folder, which is automatically generated by [Sphinx](https://www.sphinx-doc.org) - a software primarily designed for documenting Python projects and clearly abused here. Details about the build can be found in the official [Sphinx tutorial](https://www.sphinx-doc.org/en/master/tutorial/index.html). In particular, I heavily use [MathJax 3](https://www.mathjax.org/) to render math formulae. My MathJax configurations can be found in **build-source/source/conf.py**, which actually hosts all Sphinx configurations. On top of Sphinx, I've made the following personal choices:

- [Poetry](https://python-poetry.org/) is used to manage dependencies.

- [Pydata-Sphinx-Theme](https://github.com/pydata/pydata-sphinx-theme) is used to style the website.

- [Inkscape](https://inkscape.org) is used to draw pictures in the format of SVG files.

- Since Github page can only be rendered under certain path (I use **main/docs**), I followed suggestions in this [post](https://www.docslikecode.com/articles/github-pages-python-sphinx/) to organize folders in this repo.

## Work with the repo
- Run `poetry install` from the root folder to install the dependencies.

- Run `poetry run make livehtml` from the folder **build-source** to view the changes live on `localhost:8000`.

- Run `poetry run make github` from the folder **build-source** to create the source files for the GitHub-page and copy them to the **docs** folder so that it's ready for the gh-page.
