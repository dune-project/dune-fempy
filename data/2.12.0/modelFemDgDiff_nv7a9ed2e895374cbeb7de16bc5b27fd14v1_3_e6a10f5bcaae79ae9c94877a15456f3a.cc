#ifndef GuardModelImpl_nv7a9ed2e895374cbeb7de16bc5b27fd14v1_3_e6a10f5bcaae79ae9c94877a15456f3a
#define GuardModelImpl_nv7a9ed2e895374cbeb7de16bc5b27fd14v1_3_e6a10f5bcaae79ae9c94877a15456f3a
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

namespace ModelImpl_nv7a9ed2e895374cbeb7de16bc5b27fd14v1_3_e6a10f5bcaae79ae9c94877a15456f3a
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
    using ConstantType = typename std::tuple_element_t< i, std::tuple< std::shared_ptr< double > > >::element_type;
    static const std::size_t numConstants = 1;
    template< std::size_t i >
    using CoefficientType = std::tuple_element_t< i, std::tuple< Dune::Fem::ConstLocalFunction<Velocity0>  > >;
    typedef std::tuple< typename CoefficientType<0>::FunctionSpaceType > CoefficientFunctionSpaceTupleType;

    Model ( const Velocity0 &coefficient0, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
      : coefficients_(CoefficientType<0>(coefficient0) )
    {
      std::get< 0 >( constants_ ) = std::make_shared< double >();
    }

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
      double tmp0 = constant< 0 >();
        assert(tmp0 == tmp0);
      const auto tmp1 = -1 * (du[ 0 ])[ 0 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = (du[ 0 ])[ 0 ] + tmp1;
        assert(tmp2 == tmp2);
      const auto tmp3 = tmp2 * tmp0;
        assert(tmp3 == tmp3);
      const auto tmp4 = -1 * (du[ 0 ])[ 1 ];
        assert(tmp4 == tmp4);
      const auto tmp5 = (du[ 0 ])[ 1 ] + tmp4;
        assert(tmp5 == tmp5);
      const auto tmp6 = tmp5 * tmp0;
        assert(tmp6 == tmp6);
      (result[ 0 ])[ 0 ] = tmp3;
      (result[ 0 ])[ 1 ] = tmp6;
    }

    template< class Point >
    void linFlux ( const DRangeType &ubar, const DJacobianRangeType &dubar, const Point &x, const DRangeType &u, const DJacobianRangeType &du, RJacobianRangeType &result ) const
    {
      double tmp0 = constant< 0 >();
        assert(tmp0 == tmp0);
      const auto tmp1 = -1 * (du[ 0 ])[ 0 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = (du[ 0 ])[ 0 ] + tmp1;
        assert(tmp2 == tmp2);
      const auto tmp3 = tmp2 * tmp0;
        assert(tmp3 == tmp3);
      const auto tmp4 = -1 * (du[ 0 ])[ 1 ];
        assert(tmp4 == tmp4);
      const auto tmp5 = (du[ 0 ])[ 1 ] + tmp4;
        assert(tmp5 == tmp5);
      const auto tmp6 = tmp5 * tmp0;
        assert(tmp6 == tmp6);
      (result[ 0 ])[ 0 ] = tmp3;
      (result[ 0 ])[ 1 ] = tmp6;
    }

    template< class Point >
    void fluxDivergence ( const Point &x, const DRangeType &u, const DJacobianRangeType &du, const DHessianRangeType &d2u, RRangeType &result ) const
    {
      double tmp0 = constant< 0 >();
        assert(tmp0 == tmp0);
      const auto tmp1 = -1 * ((d2u[ 0 ])[ 1 ])[ 1 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = ((d2u[ 0 ])[ 1 ])[ 1 ] + tmp1;
        assert(tmp2 == tmp2);
      const auto tmp3 = tmp2 * tmp0;
        assert(tmp3 == tmp3);
      const auto tmp4 = -1 * ((d2u[ 0 ])[ 0 ])[ 0 ];
        assert(tmp4 == tmp4);
      const auto tmp5 = ((d2u[ 0 ])[ 0 ])[ 0 ] + tmp4;
        assert(tmp5 == tmp5);
      const auto tmp6 = tmp5 * tmp0;
        assert(tmp6 == tmp6);
      const auto tmp7 = tmp6 + tmp3;
        assert(tmp7 == tmp7);
      const auto tmp8 = -1 * tmp7;
        assert(tmp8 == tmp8);
      result[ 0 ] = tmp8;
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
    const ConstantType< i > &constant () const
    {
      return *std::get< i >( constants_ );
    }

    template< std::size_t i >
    ConstantType< i > &constant ()
    {
      return *std::get< i >( constants_ );
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

    const double &time () const
    {
      return *std::get< 0 >( constants_ );
    }

    double &time ()
    {
      return *std::get< 0 >( constants_ );
    }

  private:
    mutable const EntityType *entity_ = nullptr;
    mutable std::tuple< std::shared_ptr< double > > constants_;
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

} // namespace ModelImpl_nv7a9ed2e895374cbeb7de16bc5b27fd14v1_3_e6a10f5bcaae79ae9c94877a15456f3a

namespace ModelImpl_nv7a9ed2e895374cbeb7de16bc5b27fd14v1_3_e6a10f5bcaae79ae9c94877a15456f3a
{

  typedef ModelImpl_nv7a9ed2e895374cbeb7de16bc5b27fd14v1_3_e6a10f5bcaae79ae9c94877a15456f3a::Model< typename Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction<typename Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, Dune::FieldVector< double, 2 > > > ModelType;

  std::size_t renumberConstants ( pybind11::handle &obj )
  {
    std::string id = pybind11::str( obj );
    if( pybind11::hasattr(obj,"name") ) id = pybind11::str(obj.attr("name"));
    if (id == "time") return 0;
    throw pybind11::value_error("coefficient '" + id + "' has not been registered");
  }

  template< std::size_t i >
  void setConstant ( ModelType &model, pybind11::handle value )
  {
    model.template constant< i >() = value.template cast< typename ModelType::ConstantType< i > >();
  }

  template< std::size_t... i >
  auto DUNE_PRIVATE defSetConstant ( std::index_sequence< i... > )
  {
    typedef std::function< void( ModelType &model, pybind11::handle ) > Dispatch;
    std::array< Dispatch, sizeof...( i ) > dispatch = {{ Dispatch( setConstant< i > )... }};

    return [ dispatch ] ( ModelImpl_nv7a9ed2e895374cbeb7de16bc5b27fd14v1_3_e6a10f5bcaae79ae9c94877a15456f3a::Model< typename Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction<typename Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, Dune::FieldVector< double, 2 > > > &model, pybind11::handle coeff, pybind11::handle value ) {
        std::size_t k = renumberConstants( coeff );
        if( k >= dispatch.size() )
          throw std::range_error( "No such coefficient: "+std::to_string(k)+" >= "+std::to_string(dispatch.size()) );
        dispatch[ k ]( model, value );
        return k;
      };
  }

} // namespace ModelImpl_nv7a9ed2e895374cbeb7de16bc5b27fd14v1_3_e6a10f5bcaae79ae9c94877a15456f3a

PYBIND11_MODULE( modelFemDgDiff_nv7a9ed2e895374cbeb7de16bc5b27fd14v1_3_e6a10f5bcaae79ae9c94877a15456f3a, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView > GridPart;
  typedef ModelImpl_nv7a9ed2e895374cbeb7de16bc5b27fd14v1_3_e6a10f5bcaae79ae9c94877a15456f3a::Model< GridPart, Dune::FemPy::VirtualizedGridFunction<typename Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, Dune::FieldVector< double, 2 > > > Model;
  typedef Model ModelWrapper;
  auto cls = Dune::Python::insertClass<ModelWrapper>(module,"Model",Dune::Python::GenerateTypeName("ModelImpl_nv7a9ed2e895374cbeb7de16bc5b27fd14v1_3_e6a10f5bcaae79ae9c94877a15456f3a::Model< typename Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction<typename Dune::FemPy::GridPart< Dune::SPGrid< double, 2, Dune::SPIsotropicRefinement >::LeafGridView >, Dune::FieldVector< double, 2 > > >"),Dune::Python::IncludeFiles({"python/dune/generated/modelFemDgDiff_nv7a9ed2e895374cbeb7de16bc5b27fd14v1_3_e6a10f5bcaae79ae9c94877a15456f3a.cc"})).first;
  cls.def_property_readonly( "dimRange", [] ( ModelWrapper & ) { return 1; } );
  cls.def_property_readonly( "hasDirichletBoundary", [] ( ModelWrapper& ) -> bool { return false;});

  cls.def_property( "time", [] ( ModelWrapper &self ) { return self.template constant<0>(); }, [] ( ModelWrapper &self, typename ModelWrapper::ConstantType<0>& value) { self.template constant<0>() = value; });

  cls.def( "setConstant",ModelImpl_nv7a9ed2e895374cbeb7de16bc5b27fd14v1_3_e6a10f5bcaae79ae9c94877a15456f3a::defSetConstant( std::make_index_sequence< Model::numConstants >() ) );

  cls.def( pybind11::init( [] ( const Dune::FemPy::VirtualizedGridFunction< GridPart, Dune::FieldVector< double, 2 > > &coefficient0 ) {
    return new  ModelWrapper( coefficient0 );
    }), pybind11::keep_alive< 1, 2 >() );
}
#endif // GuardModelImpl_nv7a9ed2e895374cbeb7de16bc5b27fd14v1_3_e6a10f5bcaae79ae9c94877a15456f3a
