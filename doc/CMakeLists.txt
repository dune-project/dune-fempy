include(DuneSphinxDoc)
add_sphinx_targets(sphinx)
add_subdirectory("_templates")

dune_symlink_to_source_files(FILES
additional.md index.md twophaseflow_descr.md contributions.md
mmesh_descr.md todo.md citing.md
corepy.md installation.md vemdemo_descr.md overview.md
stokes.md furtherexamples.md nextsteps.md furthertopics.md
extensionmodules.md inforesources.md userprojects.md
gridviews.md cpp.md femdg.md developers.md slurmbatchscript.md
changelog.md changelog290.md scheme_api.md
3dexample.py          limit.py              spiral.py
crystal.py            lineplot.py           svg2pdf.py            cppfunctions.py
dune-corepy.py        mcf-algorithm.py      twophaseflow.py       chemical.py
concepts.py           dune-fempy.py         euler.py
elasticity.py         mcf.py
uzawa-scipy.py        monolithicStokes.py   fieldsplitStokes.py
othergrids.py         gitlab-formatting.py  gmsh2dgf.py          pandoc-formatting.py
vemdemo.py            chimpl.py             discontinuousgalerkin.py
interpolation.py      wave.py               laplace-adaptive.py   laplace-dwr.py laplace-dwr-algorithm.py
backuprestore.py      parallelization.py    boundary.py
solversInternal.py    solversExternal.py    evalues_laplace.py    biharmonic_IPC0.py    mmesh.py
crystal_nb.ipynb           mcf_nb.ipynb               chemical_nb.ipynb
dune-corepy_nb.ipynb       spiral_nb.ipynb            euler_nb.ipynb
concepts_nb.ipynb          dune-fempy_nb.ipynb        twophaseflow_nb.ipynb
elasticity_nb.ipynb        uzawa-scipy_nb.ipynb       othergrids_nb.ipynb
laplace-adaptive_nb.ipynb  vemdemo_nb.ipynb           chimpl_nb.ipynb
discontinuousgalerkin_nb.ipynb gmsh2dgf_nb.ipynb
monolithicStokes_nb.ipynb  fieldsplitStokes_nb.ipynb
laplace-dwr_nb.ipynb       laplace-dwr-algorithm_nb.ipynb
lineplot_nb.ipynb          wave_nb.ipynb
mcf-algorithm_nb.ipynb     parallelization_nb.ipynb   backuprestore_nb.ipynb
boundary_nb.ipynb          mmesh_nb.ipynb
solversInternal_nb.ipynb  solversExternal_nb.ipynb
cppfunctions_nb.ipynb      laplace-dwr.hh             utility.hh
uzawa.hh                   evalues_laplace_nb.ipynb   biharmonic_IPC0_nb.ipynb
figures _static dune-fempy.bib .special.md
soap.dgf sphere.dgf unitcube-2d.dgf wave_tank.msh quads.msh
config.opts
levelgridview_nb.ipynb
filteredgridview_nb.ipynb
geoview_nb.ipynb
geoview_nb_files
mcf.gif
dune-python-logo.png
favicon.ico
)
