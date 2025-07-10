#ifndef Guard_gf_d55a0c54785f18c75ab4b1047a2df261_88d29816137b4270fb400a51f650f4b8
#define Guard_gf_d55a0c54785f18c75ab4b1047a2df261_88d29816137b4270fb400a51f650f4b8

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/grid/function.hh>
#include <dune/python/pybind11/pybind11.h>

PYBIND11_MODULE( gf_d55a0c54785f18c75ab4b1047a2df261_88d29816137b4270fb400a51f650f4b8, module )
{
  typedef pybind11::function Evaluate;
  Dune::Python::registerGridFunction< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView, Evaluate, 0 >( module, "gf", true );
}
#endif
