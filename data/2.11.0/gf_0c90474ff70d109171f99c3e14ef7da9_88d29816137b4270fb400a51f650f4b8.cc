#ifndef Guard_gf_0c90474ff70d109171f99c3e14ef7da9_88d29816137b4270fb400a51f650f4b8
#define Guard_gf_0c90474ff70d109171f99c3e14ef7da9_88d29816137b4270fb400a51f650f4b8

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/polygongrid/dgf.hh>
#include <dune/polygongrid/grid.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/grid/function.hh>
#include <dune/python/pybind11/pybind11.h>

PYBIND11_MODULE( gf_0c90474ff70d109171f99c3e14ef7da9_88d29816137b4270fb400a51f650f4b8, module )
{
  typedef pybind11::function Evaluate;
  Dune::Python::registerGridFunction< Dune::PolygonGrid< double >::LeafGridView, Evaluate, 0 >( module, "gf", true );
}
#endif
