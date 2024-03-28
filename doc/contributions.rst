.. _contributing:

#############################
Contributing to this tutorial
#############################


Contributions are very welcomed - these can include but are not restricted to

* pointing out typos or suggesting places were the description should be improved
  (open an MR in https://gitlab.dune-project.org/dune-fem/dune-fempy
  or send us an email)

* providing an improvement:

  To make changes or add something to one of the existing scripts, please
  note that the notebooks are generated from the Python scripts using the
  package `jupytext`. So any changes should be done in the scripts. To check
  the effect of any change in the notebook simply install the Python packages
  `jupytext` (and `jupyter`) and after the change run

  .. code-block:: rst

    jupytext script.py --execute --set-kernel python3 --output script_nb.ipynb

  Then you can have a look at the notebook in the usual way.
  Once you are finished open a MR or send us the modified script - any help
  will certainly be appreciated!

* showcasing what you are doing with this project:
  we would like to showcase what user have used this package for,
  by adding sections to the :ref:`userprojects` chapter:

#################################
How to showcase your own project
#################################

To provide an example of what you have been doing with this package,
send us a file *projectname_descr.rst* with a short
description of your project with the authors and links to the project web
page, journal article and so on.
The following is an example taken from the description file of the
:download:`vem project<vemdemo_descr.rst>`

.. code-block:: rst

   ########################################################
   *dune-vem*: implementation of the virtual element method
   ########################################################

   .. sectionauthor:: Andreas Dedner <a.s.dedner@warwick.ac.uk>, Martin Nolte <nolte.mrtn@gmail.com>

   This module is based on DUNE-FEM
   (https://gitlab.dune-project.org/dune-fem/dune-fem)
   and provides implementation for the Virtual Element Method.
   The code is available in the module
   https://gitlab.dune-project.org/dune-fem/dune-vem.

A python script *myproject.py* showcasing your project can be provided as
well. But it should run with a standard set of Python packages and Dune
modules - and not take ages to run. For longer running example simply link
to results on your project webpage.

The script can contain both *markdown* parts for
additional descriptions between python code blocks similar to
a jupyter notebook
(see e.g. the :download:`script<vemdemo.py>` for the virtual element project).
This file will be translated into a :download:`jupyter notebook<vemdemo_nb.ipynb>`
which will also be made available for download. The python module *jupytext* and the sphinx
extension *nbsphinx* is used to generate the restructured text file used in the
documentation. The syntax for combining markdown cells and python code
is straightforward:

.. code-block::

   # %% [markdown]
   # # Laplace problem
   #
   # We first consider a simple Laplace problem with Dirichlet boundary conditions
   # \begin{align*}
   #   -\Delta u &= f, && \text{in } \Omega, \\
   #           u &= g, && \text{on } \partial\Omega,
   # \end{align*}
   # First some setup code:

   # %%
   import dune.vem
   from dune.grid import cartesianDomain, gridFunction
   from dune.vem import voronoiCells
   from ufl import *
   import dune.ufl

   # %% [markdown]
   # Now we define the model starting with the exact solution:

   # %%
   uflSpace = dune.ufl.Space(2, dimRange=1)
   x = SpatialCoordinate(uflSpace)
   exact = as_vector( [x[0]*x[1] * cos(pi*x[0]*x[1])] )

   # next the bilinear form
   u = TrialFunction(uflSpace)
   v = TestFunction(uflSpace)
   a = (inner(grad(u),grad(v))) * dx

If you have any questions or something in unclear let us know!
