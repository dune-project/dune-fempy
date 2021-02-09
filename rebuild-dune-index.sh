#!/bin/bash

rm -rf env
rm -rf package
python3 -m virtualenv -p $(which python3) env
source env/bin/activate
python -m pip install -r requirements.txt

getDune () {
  # Clone the Dune modules that we need
  rm -rf dune
  mkdir dune
  cd dune
  git clone https://gitlab.dune-project.org/core/dune-common.git
  git clone https://gitlab.dune-project.org/core/dune-geometry.git
  git clone https://gitlab.dune-project.org/core/dune-grid.git
  git clone https://gitlab.dune-project.org/core/dune-istl.git
  git clone https://gitlab.dune-project.org/core/dune-localfunctions.git
  cd ..
  # Add packages from these Dune modules
  pip2pi packages ./dune/dune-common
  pip2pi packages ./dune/dune-geometry
  pip2pi packages ./dune/dune-grid
  pip2pi packages ./dune/dune-istl
  pip2pi packages ./dune/dune-localfunctions
}

# getDune

# Add some packages that are needed for Python 3.7, but not 3.8
# This is not really a satisfactory way of solving the problem
# that the index depends on the platform+Python version it was built on.
pip2pi packages importlib-metadata typing-extensions

# Additional packages needed to build the tutorial
# There is an issue with the correct jinja package so it is included here
# with specific lower version
pip2pi packages Pygments==2.7.3
pip2pi packages jupyterlab-pygments==0.1.2
pip2pi packages fenics-ufl # UFL==2017.1.0
pip2pi packages \
  apidoc\
  gmsh\
  ipython\
  jinja2\
  jupyter\
  mayavi\
  nbsphinx\
  portalocker\
  pandas\
  pandoc\
  pygmsh\
  pyqt5\
  sparkmonitor\
  sphinx-autobuild\
  sphinx_rtd_theme\
  sphinxcontrib-bibtex\
  triangle\
  wheel\
  x3d


# Upload all the packages
# export TWINE_USERNAME=__token__
rm -rf ./packages/simple
for filename in packages/*
do
  python -m twine upload --verbose --repository gitlab $filename
  echo upload $filename
done
