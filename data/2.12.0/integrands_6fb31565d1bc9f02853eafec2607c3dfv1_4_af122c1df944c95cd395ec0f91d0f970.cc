#ifndef GuardIntegrands_6fb31565d1bc9f02853eafec2607c3dfv1_4
#define GuardIntegrands_6fb31565d1bc9f02853eafec2607c3dfv1_4
#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fempy/py/gridview.hh>
#include <dune/python/alugrid/hierarchical.hh>
#include <dune/python/grid/gridview.hh>
#include <dune/python/grid/hierarchical.hh>
#include <cmath>
#include <tuple>
#include <utility>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/extensions.h>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/geometry/edgelength.hh>
#include <dune/fempy/geometry/celldiameter.hh>
#include <dune/fempy/py/integrands.hh>

namespace Integrands_6fb31565d1bc9f02853eafec2607c3dfv1_4
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
    typedef std::tuple<> ConstantTupleType;
    typedef std::tuple<> CoefficientTupleType;
    static constexpr bool gridPartValid = true;
    template< std::size_t i >
    using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
    template< std::size_t i >
    using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;

    Integrands ( const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    {}

    void decouple ()
    {}

    bool init ( const EntityType &entity )
    {
      entity_[ (static_cast< std::size_t >( Side::in )) ] = entity;
      geometry_.emplace( this->entity().geometry() );
      maxCellEdgeLength_[ (static_cast< std::size_t >( Side::in )) ] = Dune::FemPy::maxEdgeLength( (geometry()) );
      return true;
    }

    void unbind ()
    {}

    bool init ( const IntersectionType &intersection )
    {
      intersection_ = intersection;
      const auto& intersectionGeom = intersection_.geometry();
      maxFacetEdgeLength_ = Dune::FemPy::maxEdgeLength( (intersectionGeom) );
      entity_[ (static_cast< std::size_t >( Side::in )) ] = (intersection.inside());
      geometry_.emplace( this->entity().geometry() );
      maxCellEdgeLength_[ (static_cast< std::size_t >( Side::in )) ] = Dune::FemPy::maxEdgeLength( (geometry()) );
      if( intersection.neighbor() )
      {
        entity_[ static_cast< std::size_t >( Side::out ) ] = intersection.outside();
        const Geometry& outsideGeom = entity_[ static_cast< std::size_t >( Side::out ) ].geometry();
      maxCellEdgeLength_[ (static_cast< std::size_t >( Side::out )) ] = Dune::FemPy::maxEdgeLength( (outsideGeom) );
      }
      return true;
    }
    typedef std::tuple< Dune::FieldMatrix< double, 1, 2 >, Dune::Fem::ExplicitFieldVector< Dune::FieldMatrix< double, 2, 2 >, 1 > > DomainValueType;
    typedef std::tuple< Dune::FieldVector< double, 1 > > RangeValueType;
    static constexpr bool _nonlinear = true;

    bool nonlinear () const
    {
      return _nonlinear;
    }

    template< class Point >
    RangeValueType interior ( const Point &x, const DomainValueType &u ) const
    {
      using std::pow;
      const auto tmp0 = std::get< 1 >( u );
        assert(tmp0 == tmp0);
      const auto tmp1 = ((tmp0[ 0 ])[ 0 ])[ 0 ] + ((tmp0[ 0 ])[ 1 ])[ 1 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = std::pow( tmp1, 2 );
        assert(tmp2 == tmp2);
      auto tmp3 = maxCellEdgeLength();
        assert(tmp3 == tmp3);
      const auto tmp4 = std::pow( tmp3, 2 );
        assert(tmp4 == tmp4);
      const auto tmp5 = tmp4 * tmp2;
        assert(tmp5 == tmp5);
      return RangeValueType{ { tmp5 } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      using std::pow;
      auto tmp0 = maxCellEdgeLength();
        assert(tmp0 == tmp0);
      const auto tmp1 = std::pow( tmp0, 2 );
        assert(tmp1 == tmp1);
      const auto tmp2 = std::get< 1 >( u );
        assert(tmp2 == tmp2);
      const auto tmp3 = ((tmp2[ 0 ])[ 0 ])[ 0 ] + ((tmp2[ 0 ])[ 1 ])[ 1 ];
        assert(tmp3 == tmp3);
      const auto tmp4 = 2 * tmp3;
        assert(tmp4 == tmp4);
      const auto tmp5 = tmp4 * tmp1;
        assert(tmp5 == tmp5);
      return [ tmp5 ] ( const DomainValueType &phi ) {
          return RangeValueType{ { tmp5 * (((std::get< 1 >( phi ))[ 0 ])[ 0 ])[ 0 ] + tmp5 * (((std::get< 1 >( phi ))[ 0 ])[ 1 ])[ 1 ] } };
        };
    }

    template< class Point >
    std::pair< RangeValueType, RangeValueType > skeleton ( const Point &xIn, const DomainValueType &uIn, const Point &xOut, const DomainValueType &uOut ) const
    {
      using std::pow;
      GlobalCoordinateType tmp0 = intersection_.unitOuterNormal( xIn.localPosition() );
        assert(tmp0 == tmp0);
      const auto tmp1 = std::get< 0 >( uOut );
        assert(tmp1 == tmp1);
      const auto tmp2 = -1 * (tmp1[ 0 ])[ 1 ];
        assert(tmp2 == tmp2);
      const auto tmp3 = std::get< 0 >( uIn );
        assert(tmp3 == tmp3);
      const auto tmp4 = (tmp3[ 0 ])[ 1 ] + tmp2;
        assert(tmp4 == tmp4);
      const auto tmp5 = tmp4 * tmp0[ 1 ];
        assert(tmp5 == tmp5);
      const auto tmp6 = -1 * (tmp1[ 0 ])[ 0 ];
        assert(tmp6 == tmp6);
      const auto tmp7 = (tmp3[ 0 ])[ 0 ] + tmp6;
        assert(tmp7 == tmp7);
      const auto tmp8 = tmp7 * tmp0[ 0 ];
        assert(tmp8 == tmp8);
      const auto tmp9 = tmp8 + tmp5;
        assert(tmp9 == tmp9);
      const auto tmp10 = std::pow( tmp9, 2 );
        assert(tmp10 == tmp10);
      auto tmp11 = maxFacetEdgeLength_;
        assert(tmp11 == tmp11);
      const auto tmp12 = tmp11 * tmp10;
        assert(tmp12 == tmp12);
      const auto tmp13 = 0.5 * tmp12;
        assert(tmp13 == tmp13);
      return std::make_pair( RangeValueType{ { tmp13 } }, RangeValueType{ { tmp13 } } );
    }

    template< class Point >
    auto linearizedSkeleton ( const Point &xIn, const DomainValueType &uIn, const Point &xOut, const DomainValueType &uOut ) const
    {
      GlobalCoordinateType tmp0 = intersection_.unitOuterNormal( xIn.localPosition() );
        assert(tmp0 == tmp0);
      const auto tmp1 = 2 * tmp0[ 0 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = std::get< 0 >( uOut );
        assert(tmp2 == tmp2);
      const auto tmp3 = -1 * (tmp2[ 0 ])[ 1 ];
        assert(tmp3 == tmp3);
      const auto tmp4 = std::get< 0 >( uIn );
        assert(tmp4 == tmp4);
      const auto tmp5 = (tmp4[ 0 ])[ 1 ] + tmp3;
        assert(tmp5 == tmp5);
      const auto tmp6 = tmp0[ 1 ] * tmp5;
        assert(tmp6 == tmp6);
      const auto tmp7 = -1 * (tmp2[ 0 ])[ 0 ];
        assert(tmp7 == tmp7);
      const auto tmp8 = (tmp4[ 0 ])[ 0 ] + tmp7;
        assert(tmp8 == tmp8);
      const auto tmp9 = tmp0[ 0 ] * tmp8;
        assert(tmp9 == tmp9);
      const auto tmp10 = tmp9 + tmp6;
        assert(tmp10 == tmp10);
      const auto tmp11 = tmp10 * tmp1;
        assert(tmp11 == tmp11);
      auto tmp12 = maxFacetEdgeLength_;
        assert(tmp12 == tmp12);
      const auto tmp13 = tmp12 * tmp11;
        assert(tmp13 == tmp13);
      const auto tmp14 = 0.5 * tmp13;
        assert(tmp14 == tmp14);
      const auto tmp15 = 2 * tmp0[ 1 ];
        assert(tmp15 == tmp15);
      const auto tmp16 = tmp10 * tmp15;
        assert(tmp16 == tmp16);
      const auto tmp17 = tmp12 * tmp16;
        assert(tmp17 == tmp17);
      const auto tmp18 = 0.5 * tmp17;
        assert(tmp18 == tmp18);
      const auto tmp19 = -1 * tmp0[ 0 ];
        assert(tmp19 == tmp19);
      const auto tmp20 = 2 * tmp19;
        assert(tmp20 == tmp20);
      const auto tmp21 = tmp10 * tmp20;
        assert(tmp21 == tmp21);
      const auto tmp22 = tmp12 * tmp21;
        assert(tmp22 == tmp22);
      const auto tmp23 = 0.5 * tmp22;
        assert(tmp23 == tmp23);
      const auto tmp24 = -1 * tmp0[ 1 ];
        assert(tmp24 == tmp24);
      const auto tmp25 = 2 * tmp24;
        assert(tmp25 == tmp25);
      const auto tmp26 = tmp10 * tmp25;
        assert(tmp26 == tmp26);
      const auto tmp27 = tmp12 * tmp26;
        assert(tmp27 == tmp27);
      const auto tmp28 = 0.5 * tmp27;
        assert(tmp28 == tmp28);
      return std::make_pair( [ tmp14, tmp18 ] ( const DomainValueType &phiIn ) {
          return std::make_pair( RangeValueType{ { tmp14 * ((std::get< 0 >( phiIn ))[ 0 ])[ 0 ] + tmp18 * ((std::get< 0 >( phiIn ))[ 0 ])[ 1 ] } }, RangeValueType{ { tmp14 * ((std::get< 0 >( phiIn ))[ 0 ])[ 0 ] + tmp18 * ((std::get< 0 >( phiIn ))[ 0 ])[ 1 ] } } );
        }, [ tmp23, tmp28 ] ( const DomainValueType &phiOut ) {
          return std::make_pair( RangeValueType{ { tmp23 * ((std::get< 0 >( phiOut ))[ 0 ])[ 0 ] + tmp28 * ((std::get< 0 >( phiOut ))[ 0 ])[ 1 ] } }, RangeValueType{ { tmp23 * ((std::get< 0 >( phiOut ))[ 0 ])[ 0 ] + tmp28 * ((std::get< 0 >( phiOut ))[ 0 ])[ 1 ] } } );
        } );
    }
    typedef Dune::FieldVector< double, 1 >  RRangeType;
    typedef Dune::FieldMatrix< double, 1, GridPartType::dimension >  RJacobianRangeType;
    typedef Dune::Fem::BoundaryIdProvider< typename GridPartType::GridType > BoundaryIdProviderType;
    typedef std::array<int,1> DirichletComponentType;

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

    const EntityType &entity () const
    {
      return entity_[ (static_cast< std::size_t >( Side::in )) ];
    }

    const Geometry &geometry () const
    {
      return *geometry_;
    }

    const ctype maxCellEdgeLength () const
    {
      return maxCellEdgeLength_[ (static_cast< std::size_t >( Side::in )) ];
    }
    std::array< EntityType, 2 > entity_;
    IntersectionType intersection_;
    std::optional< Geometry > geometry_;
    std::array< ctype, 2 > maxCellEdgeLength_;
    ctype maxFacetEdgeLength_;
    ConstantTupleType constants_;
    std::array< std::tuple<  >, 2 > coefficients_;
  };

} // namespace Integrands_6fb31565d1bc9f02853eafec2607c3dfv1_4

PYBIND11_MODULE( integrands_6fb31565d1bc9f02853eafec2607c3dfv1_4_af122c1df944c95cd395ec0f91d0f970, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > > > GridPart;
  typedef Integrands_6fb31565d1bc9f02853eafec2607c3dfv1_4::Integrands< GridPart > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_6fb31565d1bc9f02853eafec2607c3dfv1_4::Integrands< GridPart >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_6fb31565d1bc9f02853eafec2607c3dfv1_4_af122c1df944c95cd395ec0f91d0f970.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] () { return new Integrands(); } ) );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return false;});
  }
}
#endif // GuardIntegrands_6fb31565d1bc9f02853eafec2607c3dfv1_4
