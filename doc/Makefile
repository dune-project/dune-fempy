SHELL := /bin/bash
PATH := bin:$(PATH)

FIGURES = figures/3dexample.png figures/interpolation_discrete.png figures/interpolation_exact.png figures/interpolation_error.png
RSTS = levelgridview_nb.rst geoview_nb.rst
BASES = spiral elasticity wave uzawa-scipy monolithicStokes fieldsplitStokes cylinder\
      laplace-adaptive crystal mmesh cppfunctions evalues_laplace\
      biharmonic_IPC0 mixed_poisson mcf boundary solversInternal solversExternal\
			mcf-algorithm lineplot laplace-dwr laplace-dwr-algorithm othergrids\
			concepts dune-fempy backuprestore parallelization dune-corepy discontinuousgalerkin\
      twophaseflow euler chemical geoview levelgridview filteredgridview\
      overview_of_advection_solver linear_diffusion_overland_flow_dune parallelrun\
      vemdemo chimpl
NBS     = $(addsuffix _nb.ipynb,$(BASES))
SCRIPTS = $(addsuffix _script.py,$(BASES))
MODULES = $(addsuffix .py.jitmodules,$(BASES))
RUN_TARGETS = $(addprefix run-,$(SCRIPTS))

# .EXPORT_ALL_VARIABLES:

export DUNE_LOGMODULES=8
export DUNE_FORCE_BUILD=1
export DUNE_LOG_LEVEL=ERROR
export ALUGRID_VERBOSITY_LEVEL=0
export IPYTHONDIR=.

# ignore SAVE_BUILD since this will appear in notebooks
export DUNE_SAVE_BUILD=none
# export DUNE_SAVE_BUILD=append

# ignore Python warnings since these will appear in notebooks
export PYTHONWARNINGS=ignore

all: $(NBS) $(FIGURES) $(SCRIPTS) $(RSTS)
modules: $(MODULES)
scripts: $(SCRIPTS)

.PHONY: clean distclean linkcheck scripts run-scripts run-scripts-report

linkcheck:
	@sphinx-build -b linkcheck . html

clean:
	@rm -f crystal*.png wave*.png *.vtu *.pvtu *.p $(FIGURES) $(NBS) $(MODULES) $(SCRIPTS) $(RSTS)
	@rm -f geoview_nb_files/geoview_nb_1_0.png
distclean: clean
	@rm -f *.bbl
	@rm -rf html

.SECONDARY:

%.py.jitmodules: %.py
	DUNEPY_DISABLE_PLOTTING=1 DUNE_LOGMODULES=4 python $<

%_script.py: %.py
	@python extractpy.py $*

%_nb.ipynb: %.py %_script.py
	@jupytext $< --execute --set-kernel python3 --output $*_nb.ipynb

%_nb.ipynb: %.ipynb
	@jupyter nbconvert --execute --to notebook --output $*_nb.ipynb $<

%_nb.rst: %_nb.ipynb
	@jupyter-nbconvert $< --to rst

run-scripts: $(RUN_TARGETS)

run-%_script.py: %_script.py
	DUNEPY_DISABLE_PLOTTING=1 python "$<" || { printf '%s\n' "$<" >> failed-scripts.txt; exit 1; }

run-scripts-report:
	@rm -f failed-scripts.txt
	@$(MAKE) -i -j8 run-scripts
	@if [[ -s failed-scripts.txt ]]; then \
	    echo "Failed scripts:"; \
	    cat failed-scripts.txt; \
	else \
	    echo "All scripts succeeded."; \
	fi

figures/3dexample.png: 3dexample.py dune-fempy_nb.ipynb
	@pvpython 3dexample.py

figures/interpolation_discrete.png figures/interpolation_exact.png figures/interpolation_error.png &: interpolation.py dune-corepy_nb.ipynb
	@pvpython interpolation.py
