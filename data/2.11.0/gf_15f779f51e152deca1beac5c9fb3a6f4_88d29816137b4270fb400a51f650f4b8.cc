#ifndef Guard_gf_15f779f51e152deca1beac5c9fb3a6f4_88d29816137b4270fb400a51f650f4b8
#define Guard_gf_15f779f51e152deca1beac5c9fb3a6f4_88d29816137b4270fb400a51f650f4b8

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fempy/py/gridview.hh>
#include <dune/python/grid/gridview.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/grid/function.hh>
#include <dune/python/pybind11/pybind11.h>

PYBIND11_MODULE( gf_15f779f51e152deca1beac5c9fb3a6f4_88d29816137b4270fb400a51f650f4b8, module )
{
  typedef pybind11::function Evaluate;
  Dune::Python::registerGridFunction< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > >, Evaluate, 0 >( module, "gf", true );
}
#endif
