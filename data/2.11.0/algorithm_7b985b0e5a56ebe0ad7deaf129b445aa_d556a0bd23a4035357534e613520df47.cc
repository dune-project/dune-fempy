#ifndef Guard_algorithm_7b985b0e5a56ebe0ad7deaf129b445aa_d556a0bd23a4035357534e613520df47
#define Guard_algorithm_7b985b0e5a56ebe0ad7deaf129b445aa_d556a0bd23a4035357534e613520df47

#include <config.h>

#define USING_DUNE_PYTHON 1


                #ifndef FEMPY_UTILITY_HH
                #define FEMPY_UTILITY_HH
                #include <tuple>
                #include <dune/fem/gridpart/common/entitysearch.hh>
                #include <dune/fempy/py/grid/gridpart.hh>

                template <class GV, class DT>
                auto sample(const GV &gv, DT &point)
                {
                  const auto& gp = Dune::FemPy::gridPart(gv);
                  Dune::Fem::EntitySearch<Dune::FemPy::GridPart<GV>>
                    search(gp);
                  const auto &entity = search(point);
                  const auto localPoint = entity.geometry().local(point);
                  return std::make_tuple(entity,localPoint);
                }
                #endif // FEMPY_UTILITY_HH


#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/common/fvector.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fempy/py/gridview.hh>
#include <dune/python/grid/gridview.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>

PYBIND11_MODULE( algorithm_7b985b0e5a56ebe0ad7deaf129b445aa_d556a0bd23a4035357534e613520df47, module )
{
  module.def( "run", [] ( Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > > & arg0, Dune::FieldVector< double, 2 > & arg1 ) {
      return sample( arg0, arg1 );
    } );
}
#endif
