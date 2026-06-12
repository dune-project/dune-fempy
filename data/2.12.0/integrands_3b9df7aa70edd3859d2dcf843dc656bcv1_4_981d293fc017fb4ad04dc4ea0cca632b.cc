#ifndef GuardIntegrands_3b9df7aa70edd3859d2dcf843dc656bcv1_4
#define GuardIntegrands_3b9df7aa70edd3859d2dcf843dc656bcv1_4
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
#include <dune/python/alugrid/hierarchical.hh>
#include <dune/python/grid/gridview.hh>
#include <dune/python/grid/hierarchical.hh>
#include <tuple>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/extensions.h>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/geometry/edgelength.hh>
#include <dune/fempy/geometry/celldiameter.hh>
#include <dune/fempy/py/integrands.hh>

namespace Integrands_3b9df7aa70edd3859d2dcf843dc656bcv1_4
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

    void decouple ()
    {
      std::get< 0 >( constants_ ) = std::make_shared< Condt >( (constant<0>()) );
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
        assert(tmp0 == tmp0);
      const auto tmp1 = -1 * tmp0[ 0 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = std::get< 0 >( u );
        assert(tmp2 == tmp2);
      const auto tmp3 = tmp2[ 0 ] + tmp1;
        assert(tmp3 == tmp3);
      const auto tmp4 = -1 * tmp0[ 1 ];
        assert(tmp4 == tmp4);
      const auto tmp5 = tmp2[ 1 ] + tmp4;
        assert(tmp5 == tmp5);
      const auto tmp6 = -1 * tmp0[ 2 ];
        assert(tmp6 == tmp6);
      const auto tmp7 = tmp2[ 2 ] + tmp6;
        assert(tmp7 == tmp7);
      double tmp8 = constant< 0 >();
        assert(tmp8 == tmp8);
      const auto tmp9 = std::get< 1 >( u );
        assert(tmp9 == tmp9);
      const auto tmp10 = 0.5 * (tmp9[ 0 ])[ 0 ];
        assert(tmp10 == tmp10);
      const auto tmp11 = 0.5 + tmp10;
        assert(tmp11 == tmp11);
      const auto tmp12 = tmp11 * tmp8;
        assert(tmp12 == tmp12);
      const auto tmp13 = 0.5 * (tmp9[ 0 ])[ 1 ];
        assert(tmp13 == tmp13);
      const auto tmp14 = tmp13 * tmp8;
        assert(tmp14 == tmp14);
      const auto tmp15 = 0.5 * (tmp9[ 0 ])[ 2 ];
        assert(tmp15 == tmp15);
      const auto tmp16 = tmp15 * tmp8;
        assert(tmp16 == tmp16);
      const auto tmp17 = 0.5 * (tmp9[ 1 ])[ 0 ];
        assert(tmp17 == tmp17);
      const auto tmp18 = tmp17 * tmp8;
        assert(tmp18 == tmp18);
      const auto tmp19 = 0.5 * (tmp9[ 1 ])[ 1 ];
        assert(tmp19 == tmp19);
      const auto tmp20 = 0.5 + tmp19;
        assert(tmp20 == tmp20);
      const auto tmp21 = tmp20 * tmp8;
        assert(tmp21 == tmp21);
      const auto tmp22 = 0.5 * (tmp9[ 1 ])[ 2 ];
        assert(tmp22 == tmp22);
      const auto tmp23 = tmp22 * tmp8;
        assert(tmp23 == tmp23);
      const auto tmp24 = 0.5 * (tmp9[ 2 ])[ 0 ];
        assert(tmp24 == tmp24);
      const auto tmp25 = tmp24 * tmp8;
        assert(tmp25 == tmp25);
      const auto tmp26 = 0.5 * (tmp9[ 2 ])[ 1 ];
        assert(tmp26 == tmp26);
      const auto tmp27 = tmp26 * tmp8;
        assert(tmp27 == tmp27);
      const auto tmp28 = 0.5 * (tmp9[ 2 ])[ 2 ];
        assert(tmp28 == tmp28);
      const auto tmp29 = 0.5 + tmp28;
        assert(tmp29 == tmp29);
      const auto tmp30 = tmp29 * tmp8;
        assert(tmp30 == tmp30);
      return RangeValueType{ { tmp3, tmp5, tmp7 }, { { tmp12, tmp14, tmp16 }, { tmp18, tmp21, tmp23 }, { tmp25, tmp27, tmp30 } } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      double tmp0 = constant< 0 >();
        assert(tmp0 == tmp0);
      const auto tmp1 = 0.5 * tmp0;
        assert(tmp1 == tmp1);
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

} // namespace Integrands_3b9df7aa70edd3859d2dcf843dc656bcv1_4

PYBIND11_MODULE( integrands_3b9df7aa70edd3859d2dcf843dc656bcv1_4_981d293fc017fb4ad04dc4ea0cca632b, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::Fem::GeometryGridPart< Dune::Fem::AdaptiveDiscreteFunction< Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< Dune::Fem::FunctionSpace< double, double, 3, 3 >, Dune::FemPy::GridPart< Dune::ALUGrid< 2, 3, Dune::simplex >::LeafGridView >, Dune::Fem::CodegenStorage> > > > GridPart;
  typedef Integrands_3b9df7aa70edd3859d2dcf843dc656bcv1_4::Integrands< GridPart > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_3b9df7aa70edd3859d2dcf843dc656bcv1_4::Integrands< GridPart >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_3b9df7aa70edd3859d2dcf843dc656bcv1_4_981d293fc017fb4ad04dc4ea0cca632b.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] () { return new Integrands(); } ) );
      cls.def_property( "dt", [] ( Integrands &self ) -> Integrands::Condt { return self.condt(); }, [] ( Integrands &self, const Integrands::Condt &v ) { self.condt() = v; } );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return false;});
  }
}
#endif // GuardIntegrands_3b9df7aa70edd3859d2dcf843dc656bcv1_4
