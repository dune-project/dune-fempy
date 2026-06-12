#ifndef GuardIntegrands_32234f2faba0d01525bf0718b15669b7v1_4
#define GuardIntegrands_32234f2faba0d01525bf0718b15669b7v1_4
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
#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/function/simplegridfunction.hh>
#include <dune/fem/misc/gridfunctionview.hh>
#include <dune/fempy/py/integrands.hh>

namespace Integrands_32234f2faba0d01525bf0718b15669b7v1_4
{




  // Integrands
  // ----------

  template< class GridPart, class Coeffdual, class Coeffdual_h >
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
    typedef std::tuple< Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > >, Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > > > CoefficientFunctionSpaceTupleType;
    typedef std::tuple< Coeffdual, Coeffdual_h > CoefficientTupleType;
    template< std::size_t i >
    using CoefficientFunctionSpaceType = std::tuple_element_t< i, CoefficientFunctionSpaceTupleType >;
    template< std::size_t i >
    using CoefficientRangeType = typename CoefficientFunctionSpaceType< i >::RangeType;
    template< std::size_t i >
    using CoefficientJacobianRangeType = typename CoefficientFunctionSpaceType< i >::JacobianRangeType;
    static constexpr bool gridPartValid = Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffdual>>() && Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffdual_h>>();
    template< std::size_t i >
    using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
    template< std::size_t i >
    using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;

    Integrands ( const Coeffdual &coeffdual, const Coeffdual_h &coeffdual_h, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
      : coefficients_{{ std::tuple< Dune::Fem::ConstLocalFunction< Coeffdual >, Dune::Fem::ConstLocalFunction< Coeffdual_h > >( Dune::Fem::ConstLocalFunction< Coeffdual >( coeffdual ), Dune::Fem::ConstLocalFunction< Coeffdual_h >( coeffdual_h ) ), std::tuple< Dune::Fem::ConstLocalFunction< Coeffdual >, Dune::Fem::ConstLocalFunction< Coeffdual_h > >( Dune::Fem::ConstLocalFunction< Coeffdual >( coeffdual ), Dune::Fem::ConstLocalFunction< Coeffdual_h >( coeffdual_h ) ) }}
    {}

    void decouple ()
    {}

    bool init ( const EntityType &entity )
    {
      entity_[ (static_cast< std::size_t >( Side::in )) ] = entity;
      std::get< 0 >( coefficients_[ static_cast< std::size_t >( Side::in ) ] ).bind( this->entity() );
      std::get< 0 >( coefficients_[ static_cast< std::size_t >( Side::out ) ] ).bind( this->entity() );
      std::get< 1 >( coefficients_[ static_cast< std::size_t >( Side::in ) ] ).bind( this->entity() );
      std::get< 1 >( coefficients_[ static_cast< std::size_t >( Side::out ) ] ).bind( this->entity() );
      return true;
    }

    void unbind ()
    {
      std::get< 0 >( coefficients_[ static_cast< std::size_t >( Side::in ) ] ).unbind( );
      std::get< 0 >( coefficients_[ static_cast< std::size_t >( Side::out ) ] ).unbind( );
      std::get< 1 >( coefficients_[ static_cast< std::size_t >( Side::in ) ] ).unbind( );
      std::get< 1 >( coefficients_[ static_cast< std::size_t >( Side::out ) ] ).unbind( );
    }

    bool init ( const IntersectionType &intersection )
    {
      intersection_ = intersection;
      entity_[ (static_cast< std::size_t >( Side::in )) ] = (intersection.inside());
      std::get< 0 >( coefficients_[ static_cast< std::size_t >( Side::in ) ] ).bind( intersection_, Side::in  );
      std::get< 1 >( coefficients_[ static_cast< std::size_t >( Side::in ) ] ).bind( intersection_, Side::in  );
      if( intersection.neighbor() )
      {
        entity_[ static_cast< std::size_t >( Side::out ) ] = intersection.outside();
      std::get< 0 >( coefficients_[ static_cast< std::size_t >( Side::out ) ] ).bind( intersection_, Side::out );
      std::get< 1 >( coefficients_[ static_cast< std::size_t >( Side::out ) ] ).bind( intersection_, Side::out );
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
      using std::abs;
      typename CoefficientFunctionSpaceType< 1 >::RangeType tmp0 = evaluateCoefficient< 1 >( x );
        assert(tmp0 == tmp0);
      const auto tmp1 = -1 * tmp0[ 0 ];
        assert(tmp1 == tmp1);
      typename CoefficientFunctionSpaceType< 0 >::RangeType tmp2 = evaluateCoefficient< 0 >( x );
        assert(tmp2 == tmp2);
      const auto tmp3 = tmp2[ 0 ] + tmp1;
        assert(tmp3 == tmp3);
      const auto tmp4 = std::abs( tmp3 );
        assert(tmp4 == tmp4);
      const auto tmp5 = std::get< 1 >( u );
        assert(tmp5 == tmp5);
      const auto tmp6 = ((tmp5[ 0 ])[ 0 ])[ 0 ] + ((tmp5[ 0 ])[ 1 ])[ 1 ];
        assert(tmp6 == tmp6);
      const auto tmp7 = std::abs( tmp6 );
        assert(tmp7 == tmp7);
      const auto tmp8 = tmp7 * tmp4;
        assert(tmp8 == tmp8);
      return RangeValueType{ { tmp8 } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      using std::abs;
      const auto tmp0 = std::get< 1 >( u );
        assert(tmp0 == tmp0);
      const auto tmp1 = ((tmp0[ 0 ])[ 0 ])[ 0 ] + ((tmp0[ 0 ])[ 1 ])[ 1 ];
        assert(tmp1 == tmp1);
      typename CoefficientFunctionSpaceType< 1 >::RangeType tmp2 = evaluateCoefficient< 1 >( x );
        assert(tmp2 == tmp2);
      const auto tmp3 = -1 * tmp2[ 0 ];
        assert(tmp3 == tmp3);
      typename CoefficientFunctionSpaceType< 0 >::RangeType tmp4 = evaluateCoefficient< 0 >( x );
        assert(tmp4 == tmp4);
      const auto tmp5 = tmp4[ 0 ] + tmp3;
        assert(tmp5 == tmp5);
      const auto tmp6 = std::abs( tmp5 );
        assert(tmp6 == tmp6);
      const auto tmp7 = tmp6 * (tmp1 == 0.0 ? 0.0 : (tmp1 < 0.0 ? -1 : 1));
        assert(tmp7 == tmp7);
      return [ tmp7 ] ( const DomainValueType &phi ) {
          return RangeValueType{ { tmp7 * (((std::get< 1 >( phi ))[ 0 ])[ 0 ])[ 0 ] + tmp7 * (((std::get< 1 >( phi ))[ 0 ])[ 1 ])[ 1 ] } };
        };
    }

    template< class Point >
    std::pair< RangeValueType, RangeValueType > skeleton ( const Point &xIn, const DomainValueType &uIn, const Point &xOut, const DomainValueType &uOut ) const
    {
      using std::abs;
      typename CoefficientFunctionSpaceType< 1 >::RangeType tmp0 = evaluateCoefficient< 1, Side::out >( xOut );
        assert(tmp0 == tmp0);
      const auto tmp1 = -1 * tmp0[ 0 ];
        assert(tmp1 == tmp1);
      typename CoefficientFunctionSpaceType< 0 >::RangeType tmp2 = evaluateCoefficient< 0, Side::out >( xOut );
        assert(tmp2 == tmp2);
      const auto tmp3 = tmp2[ 0 ] + tmp1;
        assert(tmp3 == tmp3);
      typename CoefficientFunctionSpaceType< 1 >::RangeType tmp4 = evaluateCoefficient< 1, Side::in >( xIn );
        assert(tmp4 == tmp4);
      const auto tmp5 = -1 * tmp4[ 0 ];
        assert(tmp5 == tmp5);
      typename CoefficientFunctionSpaceType< 0 >::RangeType tmp6 = evaluateCoefficient< 0, Side::in >( xIn );
        assert(tmp6 == tmp6);
      const auto tmp7 = tmp6[ 0 ] + tmp5;
        assert(tmp7 == tmp7);
      const auto tmp8 = tmp7 + tmp3;
        assert(tmp8 == tmp8);
      const auto tmp9 = 0.5 * tmp8;
        assert(tmp9 == tmp9);
      const auto tmp10 = std::abs( tmp9 );
        assert(tmp10 == tmp10);
      GlobalCoordinateType tmp11 = intersection_.unitOuterNormal( xIn.localPosition() );
        assert(tmp11 == tmp11);
      const auto tmp12 = std::get< 0 >( uOut );
        assert(tmp12 == tmp12);
      const auto tmp13 = -1 * (tmp12[ 0 ])[ 1 ];
        assert(tmp13 == tmp13);
      const auto tmp14 = std::get< 0 >( uIn );
        assert(tmp14 == tmp14);
      const auto tmp15 = (tmp14[ 0 ])[ 1 ] + tmp13;
        assert(tmp15 == tmp15);
      const auto tmp16 = tmp15 * tmp11[ 1 ];
        assert(tmp16 == tmp16);
      const auto tmp17 = -1 * (tmp12[ 0 ])[ 0 ];
        assert(tmp17 == tmp17);
      const auto tmp18 = (tmp14[ 0 ])[ 0 ] + tmp17;
        assert(tmp18 == tmp18);
      const auto tmp19 = tmp18 * tmp11[ 0 ];
        assert(tmp19 == tmp19);
      const auto tmp20 = tmp19 + tmp16;
        assert(tmp20 == tmp20);
      const auto tmp21 = std::abs( tmp20 );
        assert(tmp21 == tmp21);
      const auto tmp22 = tmp21 * tmp10;
        assert(tmp22 == tmp22);
      const auto tmp23 = 0.5 * tmp22;
        assert(tmp23 == tmp23);
      return std::make_pair( RangeValueType{ { tmp23 } }, RangeValueType{ { tmp23 } } );
    }

    template< class Point >
    auto linearizedSkeleton ( const Point &xIn, const DomainValueType &uIn, const Point &xOut, const DomainValueType &uOut ) const
    {
      using std::abs;
      typename CoefficientFunctionSpaceType< 1 >::RangeType tmp0 = evaluateCoefficient< 1, Side::out >( xOut );
        assert(tmp0 == tmp0);
      const auto tmp1 = -1 * tmp0[ 0 ];
        assert(tmp1 == tmp1);
      typename CoefficientFunctionSpaceType< 0 >::RangeType tmp2 = evaluateCoefficient< 0, Side::out >( xOut );
        assert(tmp2 == tmp2);
      const auto tmp3 = tmp2[ 0 ] + tmp1;
        assert(tmp3 == tmp3);
      typename CoefficientFunctionSpaceType< 1 >::RangeType tmp4 = evaluateCoefficient< 1, Side::in >( xIn );
        assert(tmp4 == tmp4);
      const auto tmp5 = -1 * tmp4[ 0 ];
        assert(tmp5 == tmp5);
      typename CoefficientFunctionSpaceType< 0 >::RangeType tmp6 = evaluateCoefficient< 0, Side::in >( xIn );
        assert(tmp6 == tmp6);
      const auto tmp7 = tmp6[ 0 ] + tmp5;
        assert(tmp7 == tmp7);
      const auto tmp8 = tmp7 + tmp3;
        assert(tmp8 == tmp8);
      const auto tmp9 = 0.5 * tmp8;
        assert(tmp9 == tmp9);
      const auto tmp10 = std::abs( tmp9 );
        assert(tmp10 == tmp10);
      GlobalCoordinateType tmp11 = intersection_.unitOuterNormal( xIn.localPosition() );
        assert(tmp11 == tmp11);
      const auto tmp12 = std::get< 0 >( uOut );
        assert(tmp12 == tmp12);
      const auto tmp13 = -1 * (tmp12[ 0 ])[ 1 ];
        assert(tmp13 == tmp13);
      const auto tmp14 = std::get< 0 >( uIn );
        assert(tmp14 == tmp14);
      const auto tmp15 = (tmp14[ 0 ])[ 1 ] + tmp13;
        assert(tmp15 == tmp15);
      const auto tmp16 = tmp15 * tmp11[ 1 ];
        assert(tmp16 == tmp16);
      const auto tmp17 = -1 * (tmp12[ 0 ])[ 0 ];
        assert(tmp17 == tmp17);
      const auto tmp18 = (tmp14[ 0 ])[ 0 ] + tmp17;
        assert(tmp18 == tmp18);
      const auto tmp19 = tmp18 * tmp11[ 0 ];
        assert(tmp19 == tmp19);
      const auto tmp20 = tmp19 + tmp16;
        assert(tmp20 == tmp20);
      const auto tmp21 = tmp11[ 0 ] * (tmp20 == 0.0 ? 0.0 : (tmp20 < 0.0 ? -1 : 1));
        assert(tmp21 == tmp21);
      const auto tmp22 = tmp21 * tmp10;
        assert(tmp22 == tmp22);
      const auto tmp23 = 0.5 * tmp22;
        assert(tmp23 == tmp23);
      const auto tmp24 = tmp11[ 1 ] * (tmp20 == 0.0 ? 0.0 : (tmp20 < 0.0 ? -1 : 1));
        assert(tmp24 == tmp24);
      const auto tmp25 = tmp24 * tmp10;
        assert(tmp25 == tmp25);
      const auto tmp26 = 0.5 * tmp25;
        assert(tmp26 == tmp26);
      const auto tmp27 = -1 * tmp11[ 0 ];
        assert(tmp27 == tmp27);
      const auto tmp28 = tmp27 * (tmp20 == 0.0 ? 0.0 : (tmp20 < 0.0 ? -1 : 1));
        assert(tmp28 == tmp28);
      const auto tmp29 = tmp28 * tmp10;
        assert(tmp29 == tmp29);
      const auto tmp30 = 0.5 * tmp29;
        assert(tmp30 == tmp30);
      const auto tmp31 = -1 * tmp11[ 1 ];
        assert(tmp31 == tmp31);
      const auto tmp32 = tmp31 * (tmp20 == 0.0 ? 0.0 : (tmp20 < 0.0 ? -1 : 1));
        assert(tmp32 == tmp32);
      const auto tmp33 = tmp32 * tmp10;
        assert(tmp33 == tmp33);
      const auto tmp34 = 0.5 * tmp33;
        assert(tmp34 == tmp34);
      return std::make_pair( [ tmp23, tmp26 ] ( const DomainValueType &phiIn ) {
          return std::make_pair( RangeValueType{ { tmp23 * ((std::get< 0 >( phiIn ))[ 0 ])[ 0 ] + tmp26 * ((std::get< 0 >( phiIn ))[ 0 ])[ 1 ] } }, RangeValueType{ { tmp23 * ((std::get< 0 >( phiIn ))[ 0 ])[ 0 ] + tmp26 * ((std::get< 0 >( phiIn ))[ 0 ])[ 1 ] } } );
        }, [ tmp30, tmp34 ] ( const DomainValueType &phiOut ) {
          return std::make_pair( RangeValueType{ { tmp30 * ((std::get< 0 >( phiOut ))[ 0 ])[ 0 ] + tmp34 * ((std::get< 0 >( phiOut ))[ 0 ])[ 1 ] } }, RangeValueType{ { tmp30 * ((std::get< 0 >( phiOut ))[ 0 ])[ 0 ] + tmp34 * ((std::get< 0 >( phiOut ))[ 0 ])[ 1 ] } } );
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

    template< std::size_t i, Side side, class Point >
    typename CoefficientFunctionSpaceType< i >::RangeType evaluateCoefficient ( const Point &x ) const
    {
      typename CoefficientFunctionSpaceType< i >::RangeType result;
      std::get< i >( coefficients_[ static_cast< std::size_t >( side ) ] ).evaluate( x, result );
      return result;
    }

    template< std::size_t i, class Point >
    typename CoefficientFunctionSpaceType< i >::RangeType evaluateCoefficient ( const Point &x ) const
    {
      return evaluateCoefficient< i, Side::in >( x );
    }

    template< std::size_t i, Side side, class Point >
    typename CoefficientFunctionSpaceType< i >::JacobianRangeType jacobianCoefficient ( const Point &x ) const
    {
      typename CoefficientFunctionSpaceType< i >::JacobianRangeType result;
      std::get< i >( coefficients_[ static_cast< std::size_t >( side ) ] ).jacobian( x, result );
      return result;
    }

    template< std::size_t i, class Point >
    typename CoefficientFunctionSpaceType< i >::JacobianRangeType jacobianCoefficient ( const Point &x ) const
    {
      return jacobianCoefficient< i, Side::in >( x );
    }

    template< std::size_t i, Side side, class Point >
    typename CoefficientFunctionSpaceType< i >::HessianRangeType hessianCoefficient ( const Point &x ) const
    {
      typename CoefficientFunctionSpaceType< i >::HessianRangeType result;
      std::get< i >( coefficients_[ static_cast< std::size_t >( side ) ] ).hessian( x, result );
      return result;
    }

    template< std::size_t i, class Point >
    typename CoefficientFunctionSpaceType< i >::HessianRangeType hessianCoefficient ( const Point &x ) const
    {
      return hessianCoefficient< i, Side::in >( x );
    }
    std::array< EntityType, 2 > entity_;
    IntersectionType intersection_;
    ConstantTupleType constants_;
    std::array< std::tuple< Dune::Fem::ConstLocalFunction< Coeffdual >, Dune::Fem::ConstLocalFunction< Coeffdual_h > >, 2 > coefficients_;
  };

} // namespace Integrands_32234f2faba0d01525bf0718b15669b7v1_4

PYBIND11_MODULE( integrands_32234f2faba0d01525bf0718b15669b7v1_4_af122c1df944c95cd395ec0f91d0f970, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > > > GridPart;
  typedef Integrands_32234f2faba0d01525bf0718b15669b7v1_4::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > >>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > >>, Dune::FieldVector< double, 1 > > > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_32234f2faba0d01525bf0718b15669b7v1_4::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > >>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > >>, Dune::FieldVector< double, 1 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_32234f2faba0d01525bf0718b15669b7v1_4_af122c1df944c95cd395ec0f91d0f970.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] ( const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > >>, Dune::FieldVector< double, 1 > > &coeffdual, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > >>, Dune::FieldVector< double, 1 > > &coeffdual_h ) { return new Integrands( coeffdual, coeffdual_h ); } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >() );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return false;});
  }
}
#endif // GuardIntegrands_32234f2faba0d01525bf0718b15669b7v1_4
