#ifndef Guard_algorithm_90dd6d2676806a462406f64f84948c26_a3bdc680197c909476ccd65ef123d729
#define Guard_algorithm_90dd6d2676806a462406f64f84948c26_a3bdc680197c909476ccd65ef123d729

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/common/fvector.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/space.hh>
#include <dune/python/alugrid/hierarchical.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <pybind11/pybind11.h>

#include <vector>
#include <utility>
#include <dune/fem/misc/linesegmentsampler.hh>

template <class GF, class DT>
std::pair<std::vector<DT>, std::vector<typename GF::RangeType>>
sample(const GF &gf, DT &start, DT &end, int n)
{
  Dune::Fem::LineSegmentSampler<typename GF::GridPartType> sampler(gf.gridPart(),start,end);
  std::vector<DT> coords(n);
  std::vector<typename GF::RangeType> values(n);
  sampler(gf,values);
  sampler.samplePoints(coords);
  return std::make_pair(coords,values);
}


PYBIND11_MODULE( algorithm_90dd6d2676806a462406f64f84948c26_a3bdc680197c909476ccd65ef123d729, module )
{
  module.def( "run", [] ( Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 3, 1 >, Dune::FemPy::GridPart< Dune::ALUGrid< 3, 3, Dune::simplex >::LeafGridView >, Dune::Fem::CodegenStorage> > & arg0, Dune::FieldVector< double, 3 > & arg1, Dune::FieldVector< double, 3 > & arg2, int arg3 ) {
      return sample( arg0, arg1, arg2, arg3 );
    } );
}
#endif
