#ifndef GuardModelImpl_nva984a5a51ea74a78d4c550e9740f76c1v1_3_e6a10f5bcaae79ae9c94877a15456f3a
#define GuardModelImpl_nva984a5a51ea74a78d4c550e9740f76c1v1_3_e6a10f5bcaae79ae9c94877a15456f3a
#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/grid/spgrid.hh>
#include <dune/grid/spgrid/dgfparser.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/fem/misc/boundaryidprovider.hh>

#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/extensions.h>

#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/function/virtualizedgridfunction.hh>

#include <dune/fem/schemes/conservationlawmodel.hh>

namespace ModelImpl_nva984a5a51ea74a78d4c550e9740f76c1v1_3_e6a10f5bcaae79ae9c94877a15456f3a
{




  // Model
  // -----

  template< class GridPart, class Velocity0 >
  struct Model
  {
    typedef GridPart GridPartType;
    typedef typename GridPart::template Codim< 0 >::EntityType EntityType;
    typedef typename GridPart::IntersectionType IntersectionType;
    typedef Dune::Fem::FunctionSpace< typename GridPartType::ctype, double, GridPartType::dimensionworld, 1 > DFunctionSpaceType;
    typedef typename DFunctionSpaceType::DomainFieldType DDomainFieldType;
    typedef typename DFunctionSpaceType::RangeFieldType DRangeFieldType;
    typedef typename DFunctionSpaceType::DomainType DDomainType;
    typedef typename DFunctionSpaceType::RangeType DRangeType;
    typedef typename DFunctionSpaceType::JacobianRangeType DJacobianRangeType;
    typedef typename DFunctionSpaceType::HessianRangeType DHessianRangeType;
    static const int dimDomain = GridPartType::dimensionworld;
    static const int dimD = 1;
    typedef Dune::Fem::FunctionSpace< typename GridPartType::ctype, double, GridPartType::dimensionworld, 1 > RFunctionSpaceType;
    typedef typename RFunctionSpaceType::DomainFieldType RDomainFieldType;
    typedef typename RFunctionSpaceType::RangeFieldType RRangeFieldType;
    typedef typename RFunctionSpaceType::DomainType RDomainType;
    typedef typename RFunctionSpaceType::RangeType RRangeType;
    typedef typename RFunctionSpaceType::JacobianRangeType RJacobianRangeType;
    typedef typename RFunctionSpaceType::HessianRangeType RHessianRangeType;
    static const int dimR = 1;
    static const int dimLocal = GridPartType::dimension;
    template< std::size_t i >
    using CoefficientType = std::tuple_element_t< i, std::tuple< Dune::Fem::ConstLocalFunction<Velocity0>  > >;
    typedef std::tuple< typename CoefficientType<0>::FunctionSpaceType > CoefficientFunctionSpaceTupleType;

    Model ( const Velocity0 &coefficient0, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
      : coefficients_(CoefficientType<0>(coefficient0) )
    {}

    bool init ( const EntityType &entity ) const
    {
      {
        entity_ = &entity;
        std::get< 0 >( coefficients_).bind( this->entity() );
      }
      return true;
    }

    void unbind () const
    {
      {
        std::get< 0 >( coefficients_).unbind( );
        entity_ = nullptr;
      }
    }

    const EntityType &entity () const
    {
      return *entity_;
    }

    std::string name () const
    {
      return "Model";
    }
    typedef Dune::Fem::BoundaryIdProvider< typename GridPartType::GridType > BoundaryIdProviderType;
    static const bool symmetric = false;

    template< class Point >
    void source ( const Point &x, const DRangeType &u, const DJacobianRangeType &du, RRangeType &result ) const
    {
      result[ 0 ] = 0.0;
    }

    template< class Point >
    void linSource ( const DRangeType &ubar, const DJacobianRangeType &dubar, const Point &x, const DRangeType &u, const DJacobianRangeType &du, RRangeType &result ) const
    {
      result[ 0 ] = 0.0;
    }

    template< class Point >
    void flux ( const Point &x, const DRangeType &u, const DJacobianRangeType &du, RJacobianRangeType &result ) const
    {
      const auto tmp0 = [ this ] ( auto x ) {
          typename std::tuple_element_t< 0, CoefficientFunctionSpaceTupleType >::RangeType result;
          std::get< 0 >( coefficients_ ).evaluate( x, result );
          return result;
        }( (x) );
        assert(tmp0 == tmp0);
      const auto tmp1 = u[ 0 ] * tmp0[ 0 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = u[ 0 ] * tmp0[ 1 ];
        assert(tmp2 == tmp2);
      (result[ 0 ])[ 0 ] = tmp1;
      (result[ 0 ])[ 1 ] = tmp2;
    }

    template< class Point >
    void linFlux ( const DRangeType &ubar, const DJacobianRangeType &dubar, const Point &x, const DRangeType &u, const DJacobianRangeType &du, RJacobianRangeType &result ) const
    {
      const auto tmp0 = [ this ] ( auto x ) {
          typename std::tuple_element_t< 0, CoefficientFunctionSpaceTupleType >::RangeType result;
          std::get< 0 >( coefficients_ ).evaluate( x, result );
          return result;
        }( (x) );
        assert(tmp0 == tmp0);
      const auto tmp1 = u[ 0 ] * tmp0[ 0 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = u[ 0 ] * tmp0[ 1 ];
        assert(tmp2 == tmp2);
      (result[ 0 ])[ 0 ] = tmp1;
      (result[ 0 ])[ 1 ] = tmp2;
    }

    template< class Point >
    void fluxDivergence ( const Point &x, const DRangeType &u, const DJacobianRangeType &du, const DHessianRangeType &d2u, RRangeType &result ) const
    {
      const auto tmp0 = [ this ] ( auto x ) {
          typename std::tuple_element_t< 0, CoefficientFunctionSpaceTupleType >::JacobianRangeType result;
          std::get< 0 >( coefficients_ ).jacobian( x, result );
          return result;
        }( (x) );
        assert(tmp0 == tmp0);
      const auto tmp1 = u[ 0 ] * (tmp0[ 1 ])[ 1 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = [ this ] ( auto x ) {
          typename std::tuple_element_t< 0, CoefficientFunctionSpaceTupleType >::RangeType result;
          std::get< 0 >( coefficients_ ).evaluate( x, result );
          return result;
        }( (x) );
        assert(tmp2 == tmp2);
      const auto tmp3 = (du[ 0 ])[ 1 ] * tmp2[ 1 ];
        assert(tmp3 == tmp3);
      const auto tmp4 = tmp3 + tmp1;
        assert(tmp4 == tmp4);
      const auto tmp5 = u[ 0 ] * (tmp0[ 0 ])[ 0 ];
        assert(tmp5 == tmp5);
      const auto tmp6 = tmp2[ 0 ] * (du[ 0 ])[ 0 ];
        assert(tmp6 == tmp6);
      const auto tmp7 = tmp6 + tmp5;
        assert(tmp7 == tmp7);
      const auto tmp8 = tmp7 + tmp4;
        assert(tmp8 == tmp8);
      const auto tmp9 = -1 * tmp8;
        assert(tmp9 == tmp9);
      result[ 0 ] = tmp9;
    }

    template< class Point >
    [[deprecated]] void diffusiveFlux ( const Point &x, const DRangeType &u, const DJacobianRangeType &du, RJacobianRangeType &result ) const
    {
      flux(x, u, du, result );
    }

    template< class Point >
    [[deprecated]] void linDiffusiveFlux ( const DRangeType &ubar, const DJacobianRangeType &dubar, const Point &x, const DRangeType &u, const DJacobianRangeType &du, RJacobianRangeType &result ) const
    {
      linFlux( ubar, dubar, x, u, du, result );
    }

    template< std::size_t i >
    const CoefficientType< i > &coefficient () const
    {
      return std::get< i >( coefficients_ );
    }

    template< std::size_t i >
    CoefficientType< i > &coefficient ()
    {
      return std::get< i >( coefficients_ );
    }

  private:
    mutable const EntityType *entity_ = nullptr;
    mutable std::tuple< Dune::Fem::ConstLocalFunction<Velocity0>  > coefficients_;

  public:
    static constexpr bool hasGamma = false;

    template< class State >
    inline Dune::FieldVector<double, 2> pressureTemperature ( const State &u ) const
    {
      Dune::FieldVector<double, 2> result;
      result = Dune::FieldVector<double, 2>( (0) );
      return result;
    }

    template< class Entity, class Point >
    inline RRangeType background ( const double &t, const Entity &entity, const Point &x ) const
    {
      RRangeType result;
      result = RRangeType( (0) );
      return result;
    }

    template< class Entity, class Point >
    inline double maxWaveSpeed ( const double &t, const Entity &entity, const Point &x, const DDomainType &normal, const DRangeType &u ) const
    {
      double result;
      using std::abs;
      const auto tmp0 = [ this ] ( auto x ) {
          typename std::tuple_element_t< 0, CoefficientFunctionSpaceTupleType >::RangeType result;
          std::get< 0 >( coefficients_ ).evaluate( x, result );
          return result;
        }( (x) );
        assert(tmp0 == tmp0);
      const auto tmp1 = normal[ 1 ] * tmp0[ 1 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = normal[ 0 ] * tmp0[ 0 ];
        assert(tmp2 == tmp2);
      const auto tmp3 = tmp2 + tmp1;
        assert(tmp3 == tmp3);
      const auto tmp4 = std::abs( tmp3 );
        assert(tmp4 == tmp4);
      result = tmp4;
      return result;
    }

    template< class Entity, class Point >
    inline DDomainType velocity ( const double &t, const Entity &entity, const Point &x, const DRangeType &u ) const
    {
      DDomainType result;
      const auto tmp0 = [ this ] ( auto x ) {
          typename std::tuple_element_t< 0, CoefficientFunctionSpaceTupleType >::RangeType result;
          std::get< 0 >( coefficients_ ).evaluate( x, result );
          return result;
        }( (x) );
        assert(tmp0 == tmp0);
      result[ 0 ] = tmp0[ 0 ];
      result[ 1 ] = tmp0[ 1 ];
      return result;
    }

    template< class Entity, class Point >
    inline double diffusionTimeStep ( const Entity& entity, const Point &x, const DRangeType &u ) const
    {
      double result;
      result = double( (0) );
      return result;
    }
    static constexpr bool hasPhysical = false;

    template< class Entity, class Point >
    inline double physical ( const Entity &entity, const Point &x, const DRangeType &u ) const
    {
      double result;
      const auto tmp0 = (u[ 0 ] <= 1 ? 1 : 0.0) * (u[ 0 ] >= 0.0 ? 1 : 0.0);
        assert(tmp0 == tmp0);
      result = tmp0;
      return result;
    }

    template< class Intersection, class Point >
    inline double jump ( const Intersection& it, const Point &x, const DRangeType &u, const DRangeType &w ) const
    {
      double result;
      result = double( (0) );
      return result;
    }

    template< class Entity, class Point >
    inline void adjustAverageValue ( const Entity& entity, const Point &x, DRangeType &u ) const
    {
      {}
    }

    template< class Entity, class Point >
    inline bool boundaryFlux ( const int bndId, const double &t, const Entity& entity, const Point &x, const DDomainType &normal, const DRangeType &u, RRangeType &result ) const
    {
      switch( bndId )
      {
      default:
        {
          return false;
        }
      }
    }

    template< class Entity, class Point >
    inline bool diffusionBoundaryFlux ( const int bndId, const double &t, const Entity& entity, const Point &x, const DDomainType &normal, const DRangeType &u, const DJacobianRangeType &du, RRangeType &result ) const
    {
      switch( bndId )
      {
      default:
        {
          return false;
        }
      }
    }

    template< class Entity, class Point >
    inline bool hasBoundaryValue ( const int bndId, const double &t, const Entity& entity, const Point &x, const DRangeType &u, RRangeType &result ) const
    {
      switch( bndId )
      {
      case 1:
        {
          result[ 0 ] = 1;
          return true;
        }
        break;
      case 2:
        {
          result[ 0 ] = 1;
          return true;
        }
        break;
      case 3:
        {
          result[ 0 ] = 1;
          return true;
        }
        break;
      case 4:
        {
          result[ 0 ] = 1;
          return true;
        }
        break;
      default:
        {
          return false;
        }
      }
    }

    template< class Entity, class Point >
    inline bool boundaryValue ( const int bndId, const double &t, const Entity& entity, const Point &x, const DDomainType &normal, const DRangeType &u, RRangeType &result ) const
    {
      switch( bndId )
      {
      case 1:
        {
          result[ 0 ] = 0.0;
          return true;
        }
        break;
      case 2:
        {
          result[ 0 ] = 0.0;
          return true;
        }
        break;
      case 3:
        {
          result[ 0 ] = 0.0;
          return true;
        }
        break;
      case 4:
        {
          result[ 0 ] = 0.0;
          return true;
        }
        break;
      default:
        {
          return false;
        }
      }
    }

    template< class LimitedRange >
    inline void limitedRange ( LimitedRange& limRange ) const
    {
      {}
    }

    template< class LimitedRange >
    void obtainBounds ( LimitedRange& globalMin, LimitedRange& globalMax ) const
    {
      globalMin = {0};
      globalMax = {1};
    }
  };

} // namespace ModelImpl_nva984a5a51ea74a78d4c550e9740f76c1v1_3_e6a10f5bcaae79ae9c94877a15456f3a

namespace ModelImpl_nva984a5a51ea74a78d4c550e9740f76c1v1_3_e6a10f5bcaae79ae9c94877a15456f3a
{

} // namespace ModelImpl_nva984a5a51ea74a78d4c550e9740f76c1v1_3_e6a10f5bcaae79ae9c94877a15456f3a

PYBIND11_MODULE( modelFemDgAdv_nva984a5a51ea74a78d4c550e9740f76c1v1_3_e6a10f5bcaae79ae9c94877a15456f3a, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView > GridPart;
  typedef ModelImpl_nva984a5a51ea74a78d4c550e9740f76c1v1_3_e6a10f5bcaae79ae9c94877a15456f3a::Model< GridPart, Dune::FemPy::VirtualizedGridFunction<typename Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, Dune::FieldVector< double, 2 > > > Model;
  typedef Model ModelWrapper;
  auto cls = Dune::Python::insertClass<ModelWrapper>(module,"Model",Dune::Python::GenerateTypeName("ModelImpl_nva984a5a51ea74a78d4c550e9740f76c1v1_3_e6a10f5bcaae79ae9c94877a15456f3a::Model< typename Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction<typename Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, Dune::FieldVector< double, 2 > > >"),Dune::Python::IncludeFiles({"python/dune/generated/modelFemDgAdv_nva984a5a51ea74a78d4c550e9740f76c1v1_3_e6a10f5bcaae79ae9c94877a15456f3a.cc"})).first;
  cls.def_property_readonly( "dimRange", [] ( ModelWrapper & ) { return 1; } );
  cls.def_property_readonly( "hasDirichletBoundary", [] ( ModelWrapper& ) -> bool { return false;});



  cls.def( pybind11::init( [] ( const Dune::FemPy::VirtualizedGridFunction< GridPart, Dune::FieldVector< double, 2 > > &coefficient0 ) {
    return new  ModelWrapper( coefficient0 );
    }), pybind11::keep_alive< 1, 2 >() );
}
#endif // GuardModelImpl_nva984a5a51ea74a78d4c550e9740f76c1v1_3_e6a10f5bcaae79ae9c94877a15456f3a
