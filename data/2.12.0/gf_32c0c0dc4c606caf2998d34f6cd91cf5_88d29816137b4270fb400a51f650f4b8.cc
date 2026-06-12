#ifndef Guard_gf_32c0c0dc4c606caf2998d34f6cd91cf5_88d29816137b4270fb400a51f650f4b8
#define Guard_gf_32c0c0dc4c606caf2998d34f6cd91cf5_88d29816137b4270fb400a51f650f4b8

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/python/alugrid/hierarchical.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/grid/function.hh>
#include <dune/python/pybind11/pybind11.h>

PYBIND11_MODULE( gf_32c0c0dc4c606caf2998d34f6cd91cf5_88d29816137b4270fb400a51f650f4b8, module )
{
  typedef pybind11::function Evaluate;
  Dune::Python::registerGridFunction< Dune::ALUGrid< 2, 2, Dune::cube >::LeafGridView, Evaluate, 0 >( module, "gf", true );
}
#endif
