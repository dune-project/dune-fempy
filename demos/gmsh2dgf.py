# %% [markdown]
#
# ## Converting gmsh to DGF
#
# A simple gmsh to DGF converter allowing to specify boundary conditions.
# A complete description of the DGF format is found [here](https://dune-project.org/doxygen/master/group__DuneGridFormatParser.html#details).
#
# %%
def gmsh2DGF(points, cells, bndDomain = None, periodic = None, dim = None):
    """
    Parameter:
       points     array(list) of points of length dim
       cells      dict containing element vertex numbers
       bndDomain  dict id -> str containing possible boundary domains
       periodic   string containing periodic boundary transformation
       dim        dimension of grid

    Returns
        String containing the mesh description in DGF format.
    """

    if dim is None:
        if "tetra" in cells or "hexa" in cells:
            dim = 3
        elif "quad" in cells or "triangle" in cells:
            dim = 2
        else:
            dim = len(points[0])

    simplex = "triangle" if "triangle" in cells else None
    if dim == 3 and "tetra" in cells:
        simplex = "tetra"

    dgf="DGF\nVertex\n"
    for p in points:
        for i in range(dim):
            dgf += str(p[i]) + " "
        dgf += "\n"
    dgf += "#\n\n"

    if simplex is not None:
        dgf += "Simplex\n"
        for t in cells[simplex]:
            for v in t:
                dgf += str(v) + " "
            dgf += "\n"
        dgf += "#\n\n"

    if "quad" in cells:
        dgf += "Cube\n"
        # gmsh has a different reference quadrilateral
        vxmap = [0,1,3,2] # flip vertex 2 and 3
        for t in cells["quad"]:
            for i in range(4):
                dgf += str(t[vxmap[i]]) + " "
            dgf += "\n"
        dgf += "#\n\n"

    if bndDomain is not None:
        assert isinstance(bndDomain, dict), "Expecting a dictionary for boundary domain"
        dgf += "BoundaryDomain\n"
        for bndid,bnd in bndDomain.items():
            if bnd == 'default':
                dgf += bnd + " " + str(bndid) +"\n"
            else:
                dgf += str(bndid) + " " + bnd +"\n"

        dgf += "#\n\n"
    if periodic is not None:
        dgf += periodic

    return dgf
