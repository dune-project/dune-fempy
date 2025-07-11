#ifndef Guard_gf_873f0dfcf004423db817995079003b90_88d29816137b4270fb400a51f650f4b8
#define Guard_gf_873f0dfcf004423db817995079003b90_88d29816137b4270fb400a51f650f4b8

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/python/grid/hierarchical.hh>
#include <dune/vem/misc/grid.hh>

#include <dune/python/grid/function.hh>
#include <dune/python/pybind11/pybind11.h>

PYBIND11_MODULE( gf_873f0dfcf004423db817995079003b90_88d29816137b4270fb400a51f650f4b8, module )
{
  typedef pybind11::function Evaluate;
  Dune::Python::registerGridFunction< Dune::Vem::Grid<2,2>::LeafGridView, Evaluate, 0 >( module, "gf", true );
}
#endif
