#ifndef GuardIntegrands_dfb7d77952f419a7b0d1b99691ff3a05v1_3
#define GuardIntegrands_dfb7d77952f419a7b0d1b99691ff3a05v1_3
#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/gridpart/geometrygridpart.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/gridview.hh>
#include <dune/fempy/py/space.hh>
#include <dune/python/grid/gridview.hh>
#include <dune/python/grid/hierarchical.hh>
#include <tuple>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/extensions.h>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/geometry/edgelength.hh>
#include <dune/fempy/py/integrands.hh>

namespace Integrands_dfb7d77952f419a7b0d1b99691ff3a05v1_3
{




  // Integrands
  // ----------

  template< class GridPart >
  struct Integrands
  {
    typedef GridPart GridPartType;
    typedef typename GridPartType::GridViewType GridView;
    typedef typename GridView::ctype ctype;
    typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
    typedef typename GridPartType::IntersectionType IntersectionType;
    typedef typename EntityType::Geometry Geometry;
    typedef typename Geometry::GlobalCoordinate GlobalCoordinateType;
    typedef Dune::Fem::IntersectionSide Side;
    typedef double Condt;
    typedef std::tuple< std::shared_ptr< Condt > > ConstantTupleType;
    template< std::size_t i >
    using ConstantsRangeType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;
    typedef std::tuple<> CoefficientTupleType;
    static constexpr bool gridPartValid = true;
    template< std::size_t i >
    using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
    template< std::size_t i >
    using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;

    Integrands ( const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    {
      std::get< 0 >( constants_ ) = std::make_shared< Condt >( (Condt(0)) );
    }

    bool init ( const EntityType &entity )
    {
      entity_ = entity;
      geometry_.emplace( this->entity().geometry() );
      return true;
    }

    void unbind ()
    {}

    bool init ( const IntersectionType &intersection )
    {
      intersection_ = intersection;
      return (intersection.boundary() && init( intersection.inside() ));
    }
    typedef std::tuple< Dune::FieldVector< double, 3 >, Dune::FieldMatrix< double, 3, 3 > > DomainValueType;
    typedef std::tuple< Dune::FieldVector< double, 3 >, Dune::FieldMatrix< double, 3, 3 > > RangeValueType;
    static constexpr bool _nonlinear = false;

    bool nonlinear () const
    {
      return _nonlinear;
    }

    template< class Point >
    RangeValueType interior ( const Point &x, const DomainValueType &u ) const
    {
      GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
      const auto tmp1 = -1 * tmp0[ 0 ];
      const auto tmp2 = std::get< 0 >( u );
      const auto tmp3 = tmp2[ 0 ] + tmp1;
      const auto tmp4 = -1 * tmp0[ 1 ];
      const auto tmp5 = tmp2[ 1 ] + tmp4;
      const auto tmp6 = -1 * tmp0[ 2 ];
      const auto tmp7 = tmp2[ 2 ] + tmp6;
      const auto tmp8 = std::get< 1 >( u );
      const auto tmp9 = 0.5 * (tmp8[ 0 ])[ 0 ];
      const auto tmp10 = 0.5 + tmp9;
      double tmp11 = constant< 0 >();
      const auto tmp12 = tmp11 * tmp10;
      const auto tmp13 = 0.5 * (tmp8[ 0 ])[ 1 ];
      const auto tmp14 = tmp11 * tmp13;
      const auto tmp15 = 0.5 * (tmp8[ 0 ])[ 2 ];
      const auto tmp16 = tmp11 * tmp15;
      const auto tmp17 = 0.5 * (tmp8[ 1 ])[ 0 ];
      const auto tmp18 = tmp11 * tmp17;
      const auto tmp19 = 0.5 * (tmp8[ 1 ])[ 1 ];
      const auto tmp20 = 0.5 + tmp19;
      const auto tmp21 = tmp11 * tmp20;
      const auto tmp22 = 0.5 * (tmp8[ 1 ])[ 2 ];
      const auto tmp23 = tmp11 * tmp22;
      const auto tmp24 = 0.5 * (tmp8[ 2 ])[ 0 ];
      const auto tmp25 = tmp11 * tmp24;
      const auto tmp26 = 0.5 * (tmp8[ 2 ])[ 1 ];
      const auto tmp27 = tmp11 * tmp26;
      const auto tmp28 = 0.5 * (tmp8[ 2 ])[ 2 ];
      const auto tmp29 = 0.5 + tmp28;
      const auto tmp30 = tmp11 * tmp29;
      return RangeValueType{ { tmp3, tmp5, tmp7 }, { { tmp12, tmp14, tmp16 }, { tmp18, tmp21, tmp23 }, { tmp25, tmp27, tmp30 } } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      double tmp0 = constant< 0 >();
      const auto tmp1 = 0.5 * tmp0;
      return [ tmp1 ] ( const DomainValueType &phi ) {
          return RangeValueType{ { (std::get< 0 >( phi ))[ 0 ], (std::get< 0 >( phi ))[ 1 ], (std::get< 0 >( phi ))[ 2 ] }, { { tmp1 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ], tmp1 * ((std::get< 1 >( phi ))[ 0 ])[ 1 ], tmp1 * ((std::get< 1 >( phi ))[ 0 ])[ 2 ] }, { tmp1 * ((std::get< 1 >( phi ))[ 1 ])[ 0 ], tmp1 * ((std::get< 1 >( phi ))[ 1 ])[ 1 ], tmp1 * ((std::get< 1 >( phi ))[ 1 ])[ 2 ] }, { tmp1 * ((std::get< 1 >( phi ))[ 2 ])[ 0 ], tmp1 * ((std::get< 1 >( phi ))[ 2 ])[ 1 ], tmp1 * ((std::get< 1 >( phi ))[ 2 ])[ 2 ] } } };
        };
    }
    typedef Dune::FieldVector< double, 3 >  RRangeType;
    typedef Dune::FieldMatrix< double, 3, GridPartType::dimension >  RJacobianRangeType;
    typedef Dune::Fem::BoundaryIdProvider< typename GridPartType::GridType > BoundaryIdProviderType;
    typedef std::array<int,3> DirichletComponentType;

    bool hasDirichletBoundary () const
    {
      return false;
    }

    bool isDirichletIntersection ( const IntersectionType &intersection, DirichletComponentType &dirichletComponent ) const
    {
      return false;
    }

    template< class Point >
    void dirichlet ( int bndId, const Point &x, RRangeType &result ) const
    {}

    template< class Point >
    void dDirichlet ( int bndId, const Point &x, RJacobianRangeType &result ) const
    {}

    template< std::size_t i >
    const ConstantType< i > &constant () const
    {
      return *std::get< i >( constants_ );
    }

    template< std::size_t i >
    ConstantType< i > &constant ()
    {
      return *std::get< i >( constants_ );
    }

    const Condt &condt () const
    {
      return *std::get< 0 >( constants_ );
    }

    Condt &condt ()
    {
      return *std::get< 0 >( constants_ );
    }

    const EntityType &entity () const
    {
      return entity_;
    }

    const Geometry &geometry () const
    {
      return *geometry_;
    }
    EntityType entity_;
    IntersectionType intersection_;
    std::optional< Geometry > geometry_;
    ConstantTupleType constants_;
    std::tuple<  > coefficients_;
  };

} // namespace Integrands_dfb7d77952f419a7b0d1b99691ff3a05v1_3

PYBIND11_MODULE( integrands_dfb7d77952f419a7b0d1b99691ff3a05v1_3_981d293fc017fb4ad04dc4ea0cca632b, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::Fem::GeometryGridPart< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 3, 3 >, Dune::FemPy::GridPart< Dune::ALUGrid< 2, 3, Dune::simplex >::LeafGridView >, Dune::Fem::CodegenStorage> > > > GridPart;
  typedef Integrands_dfb7d77952f419a7b0d1b99691ff3a05v1_3::Integrands< GridPart > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_dfb7d77952f419a7b0d1b99691ff3a05v1_3::Integrands< GridPart >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_dfb7d77952f419a7b0d1b99691ff3a05v1_3_981d293fc017fb4ad04dc4ea0cca632b.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] () { return new Integrands(); } ) );
      cls.def_property( "dt", [] ( Integrands &self ) -> Integrands::Condt { return self.condt(); }, [] ( Integrands &self, const Integrands::Condt &v ) { self.condt() = v; } );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return false;});
  }
}
#endif // GuardIntegrands_dfb7d77952f419a7b0d1b99691ff3a05v1_3
