#ifndef Guard_algorithm_ff8afbaac1ff4fc3071d5868f6a506e1_a3bdc680197c909476ccd65ef123d729
#define Guard_algorithm_ff8afbaac1ff4fc3071d5868f6a506e1_a3bdc680197c909476ccd65ef123d729

#include <config.h>

#define USING_DUNE_PYTHON 1

#include <dune/common/fvector.hh>
#include <python/dune/generated/localfunction_aa96a17dd6cfe2b317348181ff5509f7_b385bf0bfa709f15b6f537b3344c6b29.cc>

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


PYBIND11_MODULE( algorithm_ff8afbaac1ff4fc3071d5868f6a506e1_a3bdc680197c909476ccd65ef123d729, module )
{
  module.def( "run", [] ( UFLLocalFunctions_aa96a17dd6cfe2b317348181ff5509f7::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::ALUGrid< 3, 3, Dune::simplex >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 3, 3, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 1 > > > & arg0, Dune::FieldVector< double, 3 > & arg1, Dune::FieldVector< double, 3 > & arg2, int arg3 ) {
      return sample( arg0, arg1, arg2, arg3 );
    } );
}
#endif
