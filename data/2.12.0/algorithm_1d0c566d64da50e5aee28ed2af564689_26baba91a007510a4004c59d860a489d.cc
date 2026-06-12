#ifndef Guard_algorithm_1d0c566d64da50e5aee28ed2af564689_26baba91a007510a4004c59d860a489d
#define Guard_algorithm_1d0c566d64da50e5aee28ed2af564689_26baba91a007510a4004c59d860a489d

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


            #ifndef FEMPY_UTILITY_HH
            #define FEMPY_UTILITY_HH
            #include <vector>
            #include <utility>
            #include <dune/fem/misc/linesegmentsampler.hh>
            #include <dune/fem/gridpart/common/entitysearch.hh>
            #include <dune/fem/function/localfunction/const.hh>

            template <class GF, class DT>
            std::pair<std::vector<DT>, std::vector<typename GF::RangeType>>
            sample(const GF &gf, DT &start, DT &end, int n, bool allowOutside)
            {
              Dune::Fem::LineSegmentSampler<typename GF::GridPartType> sampler(gf.gridPart(),start,end);
              std::vector<DT> coords(n);
              std::vector<typename GF::RangeType> values(n);
              if (allowOutside)
                sampler(gf,values, std::nothrow);
              else
                sampler(gf,values);
              sampler.samplePoints(coords);
              return std::make_pair(coords,values);
            }
            #endif



PYBIND11_MODULE( algorithm_1d0c566d64da50e5aee28ed2af564689_26baba91a007510a4004c59d860a489d, module )
{
  module.def( "run", [] ( Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 3, 1 >, Dune::FemPy::GridPart< Dune::ALUGrid< 3, 3, Dune::simplex >::LeafGridView >, Dune::Fem::CodegenStorage> > & arg0, Dune::FieldVector< double, 3 > & arg1, Dune::FieldVector< double, 3 > & arg2, int arg3, bool arg4 ) {
      return sample( arg0, arg1, arg2, arg3, arg4 );
    } );
}
#endif
