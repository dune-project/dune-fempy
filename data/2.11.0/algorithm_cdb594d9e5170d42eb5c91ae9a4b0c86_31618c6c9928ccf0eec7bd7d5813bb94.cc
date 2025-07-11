#ifndef Guard_algorithm_cdb594d9e5170d42eb5c91ae9a4b0c86_31618c6c9928ccf0eec7bd7d5813bb94
#define Guard_algorithm_cdb594d9e5170d42eb5c91ae9a4b0c86_31618c6c9928ccf0eec7bd7d5813bb94

#include <config.h>

#define USING_DUNE_PYTHON 1

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

#include <dune/common/fvector.hh>
#include <python/dune/generated/localfunction_acf405591bc6db0270766f36467b8da0_b385bf0bfa709f15b6f537b3344c6b29.cc>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>

PYBIND11_MODULE( algorithm_cdb594d9e5170d42eb5c91ae9a4b0c86_31618c6c9928ccf0eec7bd7d5813bb94, module )
{
  module.def( "run", [] ( UFLLocalFunctions_acf405591bc6db0270766f36467b8da0::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::ALUGrid< 3, 3, Dune::simplex >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 3, 3, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 1 > > > & arg0, Dune::FieldVector< double, 3 > & arg1, Dune::FieldVector< double, 3 > & arg2, int arg3 ) {
      return sample( arg0, arg1, arg2, arg3 );
    } );
}
#endif
