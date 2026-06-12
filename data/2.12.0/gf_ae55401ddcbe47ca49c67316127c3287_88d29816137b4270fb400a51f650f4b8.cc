#ifndef Guard_gf_ae55401ddcbe47ca49c67316127c3287_88d29816137b4270fb400a51f650f4b8
#define Guard_gf_ae55401ddcbe47ca49c67316127c3287_88d29816137b4270fb400a51f650f4b8

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/python/alugrid/hierarchical.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/grid/function.hh>
#include <dune/python/pybind11/pybind11.h>

PYBIND11_MODULE( gf_ae55401ddcbe47ca49c67316127c3287_88d29816137b4270fb400a51f650f4b8, module )
{
  typedef pybind11::function Evaluate;
  Dune::Python::registerGridFunction< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView, Evaluate, 2 >( module, "gf", false );
}
#endif
