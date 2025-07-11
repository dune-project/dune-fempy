#ifndef Guard_gf_021111853d7c16551f3a1bb36a6c1854_88d29816137b4270fb400a51f650f4b8
#define Guard_gf_021111853d7c16551f3a1bb36a6c1854_88d29816137b4270fb400a51f650f4b8

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/grid/function.hh>
#include <dune/python/pybind11/pybind11.h>

PYBIND11_MODULE( gf_021111853d7c16551f3a1bb36a6c1854_88d29816137b4270fb400a51f650f4b8, module )
{
  typedef pybind11::function Evaluate;
  Dune::Python::registerGridFunction< Dune::ALUGrid< 2, 2, Dune::cube >::LeafGridView, Evaluate, 1 >( module, "gf", false );
}
#endif
