#ifndef GuardIntegrands_5caf2b18d6da6400edecebe3f0a06e3cv1_4
#define GuardIntegrands_5caf2b18d6da6400edecebe3f0a06e3cv1_4
#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/yaspgrid.hh>
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

namespace Integrands_5caf2b18d6da6400edecebe3f0a06e3cv1_4
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
      return true;
    }

    void unbind ()
    {}

    bool init ( const IntersectionType &intersection )
    {
      intersection_ = intersection;
      entity_[ (static_cast< std::size_t >( Side::in )) ] = (intersection.inside());
      geometry_.emplace( this->entity().geometry() );
      if( intersection.neighbor() )
      {
        entity_[ static_cast< std::size_t >( Side::out ) ] = intersection.outside();
      }
      return true;
    }
    typedef std::tuple< Dune::FieldVector< double, 1 > > DomainValueType;
    typedef std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 1 > > RangeValueType;
    static constexpr bool _nonlinear = false;

    bool nonlinear () const
    {
      return _nonlinear;
    }

    template< class Point >
    RangeValueType interior ( const Point &x, const DomainValueType &u ) const
    {
      const auto tmp0 = std::get< 0 >( u );
        assert(tmp0 == tmp0);
      const auto tmp1 = -0.5 * tmp0[ 0 ];
        assert(tmp1 == tmp1);
      return RangeValueType{ { 0 }, { { tmp1 } } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      return [] ( const DomainValueType &phi ) {
          return RangeValueType{ { 0 }, { { -0.5 * (std::get< 0 >( phi ))[ 0 ] } } };
        };
    }

    template< class Point >
    RangeValueType boundary ( const Point &x, const DomainValueType &u ) const
    {
      using std::abs;
      GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
        assert(tmp0 == tmp0);
      GlobalCoordinateType tmp1 = intersection_.unitOuterNormal( x.localPosition() );
        assert(tmp1 == tmp1);
      const auto tmp2 = 0.5 * tmp1[ 0 ];
        assert(tmp2 == tmp2);
      const auto tmp3 = std::abs( tmp2 );
        assert(tmp3 == tmp3);
      const auto tmp4 = tmp2 + tmp3;
        assert(tmp4 == tmp4);
      const auto tmp5 = tmp4 / 2.0;
        assert(tmp5 == tmp5);
      const auto tmp6 = std::get< 0 >( u );
        assert(tmp6 == tmp6);
      const auto tmp7 = tmp6[ 0 ] * tmp5;
        assert(tmp7 == tmp7);
      const auto tmp8 = -1 * tmp5;
        assert(tmp8 == tmp8);
      const auto tmp9 = tmp2 + tmp8;
        assert(tmp9 == tmp9);
      const auto tmp10 = tmp9 + tmp7;
        assert(tmp10 == tmp10);
      const auto tmp11 = tmp10 * (tmp0[ 0 ] < 0.01 ? 1 : 0.0);
        assert(tmp11 == tmp11);
      return RangeValueType{ { tmp11 }, { { 0 } } };
    }

    template< class Point >
    auto linearizedBoundary ( const Point &x, const DomainValueType &u ) const
    {
      using std::abs;
      GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
        assert(tmp0 == tmp0);
      GlobalCoordinateType tmp1 = intersection_.unitOuterNormal( x.localPosition() );
        assert(tmp1 == tmp1);
      const auto tmp2 = 0.5 * tmp1[ 0 ];
        assert(tmp2 == tmp2);
      const auto tmp3 = std::abs( tmp2 );
        assert(tmp3 == tmp3);
      const auto tmp4 = tmp2 + tmp3;
        assert(tmp4 == tmp4);
      const auto tmp5 = tmp4 / 2.0;
        assert(tmp5 == tmp5);
      const auto tmp6 = tmp5 * (tmp0[ 0 ] < 0.01 ? 1 : 0.0);
        assert(tmp6 == tmp6);
      return [ tmp6 ] ( const DomainValueType &phi ) {
          return RangeValueType{ { tmp6 * (std::get< 0 >( phi ))[ 0 ] }, { { 0 } } };
        };
    }

    template< class Point >
    std::pair< RangeValueType, RangeValueType > skeleton ( const Point &xIn, const DomainValueType &uIn, const Point &xOut, const DomainValueType &uOut ) const
    {
      using std::abs;
      GlobalCoordinateType tmp0 = intersection_.unitOuterNormal( xIn.localPosition() );
        assert(tmp0 == tmp0);
      const auto tmp1 = 0.5 * tmp0[ 0 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = std::abs( tmp1 );
        assert(tmp2 == tmp2);
      const auto tmp3 = tmp1 + tmp2;
        assert(tmp3 == tmp3);
      const auto tmp4 = tmp3 / 2.0;
        assert(tmp4 == tmp4);
      const auto tmp5 = std::get< 0 >( uIn );
        assert(tmp5 == tmp5);
      const auto tmp6 = tmp5[ 0 ] * tmp4;
        assert(tmp6 == tmp6);
      const auto tmp7 = -1 * tmp0[ 0 ];
        assert(tmp7 == tmp7);
      const auto tmp8 = 0.5 * tmp7;
        assert(tmp8 == tmp8);
      const auto tmp9 = std::abs( tmp8 );
        assert(tmp9 == tmp9);
      const auto tmp10 = tmp8 + tmp9;
        assert(tmp10 == tmp10);
      const auto tmp11 = tmp10 / 2.0;
        assert(tmp11 == tmp11);
      const auto tmp12 = std::get< 0 >( uOut );
        assert(tmp12 == tmp12);
      const auto tmp13 = tmp12[ 0 ] * tmp11;
        assert(tmp13 == tmp13);
      const auto tmp14 = -1 * tmp13;
        assert(tmp14 == tmp14);
      const auto tmp15 = tmp14 + tmp6;
        assert(tmp15 == tmp15);
      const auto tmp16 = -1 * tmp15;
        assert(tmp16 == tmp16);
      return std::make_pair( RangeValueType{ { tmp15 }, { { 0 } } }, RangeValueType{ { tmp16 }, { { 0 } } } );
    }

    template< class Point >
    auto linearizedSkeleton ( const Point &xIn, const DomainValueType &uIn, const Point &xOut, const DomainValueType &uOut ) const
    {
      using std::abs;
      GlobalCoordinateType tmp0 = intersection_.unitOuterNormal( xIn.localPosition() );
        assert(tmp0 == tmp0);
      const auto tmp1 = 0.5 * tmp0[ 0 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = std::abs( tmp1 );
        assert(tmp2 == tmp2);
      const auto tmp3 = tmp1 + tmp2;
        assert(tmp3 == tmp3);
      const auto tmp4 = tmp3 / 2.0;
        assert(tmp4 == tmp4);
      const auto tmp5 = -1 * tmp4;
        assert(tmp5 == tmp5);
      const auto tmp6 = -1 * tmp0[ 0 ];
        assert(tmp6 == tmp6);
      const auto tmp7 = 0.5 * tmp6;
        assert(tmp7 == tmp7);
      const auto tmp8 = std::abs( tmp7 );
        assert(tmp8 == tmp8);
      const auto tmp9 = tmp7 + tmp8;
        assert(tmp9 == tmp9);
      const auto tmp10 = tmp9 / 2.0;
        assert(tmp10 == tmp10);
      const auto tmp11 = -1 * tmp10;
        assert(tmp11 == tmp11);
      const auto tmp12 = -1 * tmp11;
        assert(tmp12 == tmp12);
      return std::make_pair( [ tmp4, tmp5 ] ( const DomainValueType &phiIn ) {
          return std::make_pair( RangeValueType{ { tmp4 * (std::get< 0 >( phiIn ))[ 0 ] }, { { 0 } } }, RangeValueType{ { tmp5 * (std::get< 0 >( phiIn ))[ 0 ] }, { { 0 } } } );
        }, [ tmp11, tmp12 ] ( const DomainValueType &phiOut ) {
          return std::make_pair( RangeValueType{ { tmp11 * (std::get< 0 >( phiOut ))[ 0 ] }, { { 0 } } }, RangeValueType{ { tmp12 * (std::get< 0 >( phiOut ))[ 0 ] }, { { 0 } } } );
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
    std::array< EntityType, 2 > entity_;
    IntersectionType intersection_;
    std::optional< Geometry > geometry_;
    ConstantTupleType constants_;
    std::array< std::tuple<  >, 2 > coefficients_;
  };

} // namespace Integrands_5caf2b18d6da6400edecebe3f0a06e3cv1_4

PYBIND11_MODULE( integrands_5caf2b18d6da6400edecebe3f0a06e3cv1_4_e14d829faa8914e3ff98a63559ccacdd, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::YaspGrid< 1, Dune::EquidistantOffsetCoordinates< double, 1 > >::LeafGridView > GridPart;
  typedef Integrands_5caf2b18d6da6400edecebe3f0a06e3cv1_4::Integrands< GridPart > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_5caf2b18d6da6400edecebe3f0a06e3cv1_4::Integrands< GridPart >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_5caf2b18d6da6400edecebe3f0a06e3cv1_4_e14d829faa8914e3ff98a63559ccacdd.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] () { return new Integrands(); } ) );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return false;});
  }
}
#endif // GuardIntegrands_5caf2b18d6da6400edecebe3f0a06e3cv1_4
