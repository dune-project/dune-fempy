#ifndef GuardIntegrands_a7dfd384341a2a13a6d4112ae7603090v1_4
#define GuardIntegrands_a7dfd384341a2a13a6d4112ae7603090v1_4
#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/python/alugrid/hierarchical.hh>
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

namespace Integrands_a7dfd384341a2a13a6d4112ae7603090v1_4
{




  // Integrands
  // ----------

  template< class GridPart, class Coeffdiscrete_solution >
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
    typedef std::tuple< Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > > > CoefficientFunctionSpaceTupleType;
    typedef std::tuple< Coeffdiscrete_solution > CoefficientTupleType;
    template< std::size_t i >
    using CoefficientFunctionSpaceType = std::tuple_element_t< i, CoefficientFunctionSpaceTupleType >;
    template< std::size_t i >
    using CoefficientRangeType = typename CoefficientFunctionSpaceType< i >::RangeType;
    template< std::size_t i >
    using CoefficientJacobianRangeType = typename CoefficientFunctionSpaceType< i >::JacobianRangeType;
    static constexpr bool gridPartValid = Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffdiscrete_solution>>();
    template< std::size_t i >
    using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
    template< std::size_t i >
    using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;

    Integrands ( const Coeffdiscrete_solution &coeffdiscrete_solution, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
      : coefficients_{{ std::tuple< Dune::Fem::ConstLocalFunction< Coeffdiscrete_solution > >( Dune::Fem::ConstLocalFunction< Coeffdiscrete_solution >( coeffdiscrete_solution ) ), std::tuple< Dune::Fem::ConstLocalFunction< Coeffdiscrete_solution > >( Dune::Fem::ConstLocalFunction< Coeffdiscrete_solution >( coeffdiscrete_solution ) ) }}
    {}

    void decouple ()
    {}

    bool init ( const EntityType &entity )
    {
      entity_[ (static_cast< std::size_t >( Side::in )) ] = entity;
      std::get< 0 >( coefficients_[ static_cast< std::size_t >( Side::in ) ] ).bind( this->entity() );
      std::get< 0 >( coefficients_[ static_cast< std::size_t >( Side::out ) ] ).bind( this->entity() );
      return true;
    }

    void unbind ()
    {
      std::get< 0 >( coefficients_[ static_cast< std::size_t >( Side::in ) ] ).unbind( );
      std::get< 0 >( coefficients_[ static_cast< std::size_t >( Side::out ) ] ).unbind( );
    }

    bool init ( const IntersectionType &intersection )
    {
      intersection_ = intersection;
      entity_[ (static_cast< std::size_t >( Side::in )) ] = (intersection.inside());
      std::get< 0 >( coefficients_[ static_cast< std::size_t >( Side::in ) ] ).bind( intersection_, Side::in  );
      if( intersection.neighbor() )
      {
        entity_[ static_cast< std::size_t >( Side::out ) ] = intersection.outside();
      std::get< 0 >( coefficients_[ static_cast< std::size_t >( Side::out ) ] ).bind( intersection_, Side::out );
      }
      return true;
    }
    typedef std::tuple< Dune::FieldVector< double, 1 > > DomainValueType;
    typedef std::tuple< Dune::FieldVector< double, 1 > > RangeValueType;
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
      return RangeValueType{ { tmp0[ 0 ] } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      return [] ( const DomainValueType &phi ) {
          return RangeValueType{ { (std::get< 0 >( phi ))[ 0 ] } };
        };
    }

    template< class Point >
    std::pair< RangeValueType, RangeValueType > skeleton ( const Point &xIn, const DomainValueType &uIn, const Point &xOut, const DomainValueType &uOut ) const
    {
      using std::pow;
      GlobalCoordinateType tmp0 = intersection_.unitOuterNormal( xIn.localPosition() );
        assert(tmp0 == tmp0);
      const auto tmp1 = -1 * tmp0[ 1 ];
        assert(tmp1 == tmp1);
      typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp2 = jacobianCoefficient< 0, Side::out >( xOut );
        assert(tmp2 == tmp2);
      const auto tmp3 = (tmp2[ 0 ])[ 1 ] * tmp1;
        assert(tmp3 == tmp3);
      const auto tmp4 = -1 * tmp0[ 0 ];
        assert(tmp4 == tmp4);
      const auto tmp5 = (tmp2[ 0 ])[ 0 ] * tmp4;
        assert(tmp5 == tmp5);
      const auto tmp6 = tmp5 + tmp3;
        assert(tmp6 == tmp6);
      typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp7 = jacobianCoefficient< 0, Side::in >( xIn );
        assert(tmp7 == tmp7);
      const auto tmp8 = (tmp7[ 0 ])[ 1 ] * tmp0[ 1 ];
        assert(tmp8 == tmp8);
      const auto tmp9 = tmp0[ 0 ] * (tmp7[ 0 ])[ 0 ];
        assert(tmp9 == tmp9);
      const auto tmp10 = tmp9 + tmp8;
        assert(tmp10 == tmp10);
      const auto tmp11 = tmp10 + tmp6;
        assert(tmp11 == tmp11);
      const auto tmp12 = std::pow( tmp11, 2 );
        assert(tmp12 == tmp12);
      const auto tmp13 = 0.5 * tmp12;
        assert(tmp13 == tmp13);
      const auto tmp14 = -1 * tmp13;
        assert(tmp14 == tmp14);
      return std::make_pair( RangeValueType{ { tmp14 } }, RangeValueType{ { tmp14 } } );
    }

    template< class Point >
    auto linearizedSkeleton ( const Point &xIn, const DomainValueType &uIn, const Point &xOut, const DomainValueType &uOut ) const
    {
      return std::make_pair( [] ( const DomainValueType &phiIn ) {
          return std::make_pair( RangeValueType{ 0 }, RangeValueType{ 0 } );
        }, [] ( const DomainValueType &phiOut ) {
          return std::make_pair( RangeValueType{ 0 }, RangeValueType{ 0 } );
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
    std::array< std::tuple< Dune::Fem::ConstLocalFunction< Coeffdiscrete_solution > >, 2 > coefficients_;
  };

} // namespace Integrands_a7dfd384341a2a13a6d4112ae7603090v1_4

PYBIND11_MODULE( integrands_a7dfd384341a2a13a6d4112ae7603090v1_4_b1dc94f21d524d06c9380ea86cedb77b, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView > GridPart;
  typedef Integrands_a7dfd384341a2a13a6d4112ae7603090v1_4::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 1 > > > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_a7dfd384341a2a13a6d4112ae7603090v1_4::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 1 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_a7dfd384341a2a13a6d4112ae7603090v1_4_b1dc94f21d524d06c9380ea86cedb77b.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] ( const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffdiscrete_solution ) { return new Integrands( coeffdiscrete_solution ); } ), pybind11::keep_alive< 1, 2 >() );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return false;});
  }
}
#endif // GuardIntegrands_a7dfd384341a2a13a6d4112ae7603090v1_4
