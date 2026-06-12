#ifndef Guard_fieldvector_0ff7907781efafab7efd863b51b7533b
#define Guard_fieldvector_0ff7907781efafab7efd863b51b7533b

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/python/common/fvector.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

PYBIND11_MODULE( fieldvector_0ff7907781efafab7efd863b51b7533b, module )
{
  using pybind11::operator""_a;
  pybind11::module cls0 = module;
  {
    using DuneType = Dune::FieldVector< double ,3 >;
    auto cls = Dune::Python::insertClass< DuneType >( cls0, "FieldVector",pybind11::buffer_protocol(), Dune::Python::GenerateTypeName("Dune::FieldVector< double ,3 >"), Dune::Python::IncludeFiles{"dune/python/common/fvector.hh"}).first;
    Dune::Python::registerFieldVector( cls0, cls );
  }
}
#endif
