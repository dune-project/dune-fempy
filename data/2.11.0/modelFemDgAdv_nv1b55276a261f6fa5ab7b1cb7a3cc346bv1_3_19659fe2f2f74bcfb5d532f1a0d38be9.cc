#ifndef GuardModelImpl_nv1b55276a261f6fa5ab7b1cb7a3cc346bv1_3_19659fe2f2f74bcfb5d532f1a0d38be9
#define GuardModelImpl_nv1b55276a261f6fa5ab7b1cb7a3cc346bv1_3_19659fe2f2f74bcfb5d532f1a0d38be9
#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/fem/misc/boundaryidprovider.hh>

#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/extensions.h>

#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/function/virtualizedgridfunction.hh>

#include <dune/fem/schemes/conservationlawmodel.hh>

namespace ModelImpl_nv1b55276a261f6fa5ab7b1cb7a3cc346bv1_3_19659fe2f2f74bcfb5d532f1a0d38be9
{




  // Model
  // -----

  template< class GridPart, class StreamFunction >
  struct Model
  {
    typedef GridPart GridPartType;
    typedef typename GridPart::template Codim< 0 >::EntityType EntityType;
    typedef typename GridPart::IntersectionType IntersectionType;
    typedef Dune::Fem::FunctionSpace< typename GridPartType::ctype, double, GridPartType::dimensionworld, 3 > DFunctionSpaceType;
    typedef typename DFunctionSpaceType::DomainFieldType DDomainFieldType;
    typedef typename DFunctionSpaceType::RangeFieldType DRangeFieldType;
    typedef typename DFunctionSpaceType::DomainType DDomainType;
    typedef typename DFunctionSpaceType::RangeType DRangeType;
    typedef typename DFunctionSpaceType::JacobianRangeType DJacobianRangeType;
    typedef typename DFunctionSpaceType::HessianRangeType DHessianRangeType;
    static const int dimDomain = GridPartType::dimensionworld;
    static const int dimD = 3;
    typedef Dune::Fem::FunctionSpace< typename GridPartType::ctype, double, GridPartType::dimensionworld, 3 > RFunctionSpaceType;
    typedef typename RFunctionSpaceType::DomainFieldType RDomainFieldType;
    typedef typename RFunctionSpaceType::RangeFieldType RRangeFieldType;
    typedef typename RFunctionSpaceType::DomainType RDomainType;
    typedef typename RFunctionSpaceType::RangeType RRangeType;
    typedef typename RFunctionSpaceType::JacobianRangeType RJacobianRangeType;
    typedef typename RFunctionSpaceType::HessianRangeType RHessianRangeType;
    static const int dimR = 3;
    static const int dimLocal = GridPartType::dimension;
    template< std::size_t i >
    using ConstantType = typename std::tuple_element_t< i, std::tuple< std::shared_ptr< double > > >::element_type;
    static const std::size_t numConstants = 1;
    template< std::size_t i >
    using CoefficientType = std::tuple_element_t< i, std::tuple< Dune::Fem::ConstLocalFunction<StreamFunction>  > >;
    typedef std::tuple< typename CoefficientType<0>::FunctionSpaceType > CoefficientFunctionSpaceTupleType;

    Model ( const StreamFunction &coefficient0, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
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
      auto tmp0 = entity().geometry().global( Dune::Fem::coordinate( x ) );
      const auto tmp1 = -0.6283185307179586 + tmp0[ 1 ];
      const auto tmp2 = tmp1 * tmp1;
      const auto tmp3 = -0.6283185307179586 + tmp0[ 0 ];
      const auto tmp4 = tmp3 * tmp3;
      const auto tmp5 = tmp4 + tmp2;
      double tmp6 = constant< 0 >();
      const auto tmp7 = u[ 0 ] * u[ 1 ];
      const auto tmp8 = 10 * tmp7;
      const auto tmp9 = -1 * tmp8;
      const auto tmp10 = tmp9 + (tmp6 < 5 ? (tmp5 < 0.2 ? 1 : 0.0) : 0.0);
      const auto tmp11 = -5.654866776461628 + tmp0[ 1 ];
      const auto tmp12 = tmp11 * tmp11;
      const auto tmp13 = -5.654866776461628 + tmp0[ 0 ];
      const auto tmp14 = tmp13 * tmp13;
      const auto tmp15 = tmp14 + tmp12;
      const auto tmp16 = tmp9 + (tmp6 < 5 ? (tmp15 < 0.2 ? 1 : 0.0) : 0.0);
      const auto tmp17 = -2 * u[ 0 ];
      const auto tmp18 = u[ 1 ] * tmp17;
      const auto tmp19 = 10 * tmp18;
      const auto tmp20 = -1 * tmp19;
      result[ 0 ] = tmp10;
      result[ 1 ] = tmp16;
      result[ 2 ] = tmp20;
    }

    template< class Point >
    void linSource ( const DRangeType &ubar, const DJacobianRangeType &dubar, const Point &x, const DRangeType &u, const DJacobianRangeType &du, RRangeType &result ) const
    {
      const auto tmp0 = u[ 0 ] * ubar[ 1 ];
      const auto tmp1 = ubar[ 0 ] * u[ 1 ];
      const auto tmp2 = tmp1 + tmp0;
      const auto tmp3 = 10 * tmp2;
      const auto tmp4 = -1 * tmp3;
      const auto tmp5 = -2 * ubar[ 0 ];
      const auto tmp6 = u[ 1 ] * tmp5;
      const auto tmp7 = -2 * u[ 0 ];
      const auto tmp8 = ubar[ 1 ] * tmp7;
      const auto tmp9 = tmp8 + tmp6;
      const auto tmp10 = 10 * tmp9;
      const auto tmp11 = -1 * tmp10;
      result[ 0 ] = tmp4;
      result[ 1 ] = tmp4;
      result[ 2 ] = tmp11;
    }

    template< class Point >
    void flux ( const Point &x, const DRangeType &u, const DJacobianRangeType &du, RJacobianRangeType &result ) const
    {
      const auto tmp0 = [ this ] ( auto x ) {
          typename std::tuple_element_t< 0, CoefficientFunctionSpaceTupleType >::JacobianRangeType result;
          std::get< 0 >( coefficients_ ).jacobian( x, result );
          return result;
        }( (x) );
      const auto tmp1 = -1 * (tmp0[ 0 ])[ 1 ];
      const auto tmp2 = u[ 0 ] * tmp1;
      const auto tmp3 = u[ 0 ] * (tmp0[ 0 ])[ 0 ];
      const auto tmp4 = u[ 1 ] * tmp1;
      const auto tmp5 = (tmp0[ 0 ])[ 0 ] * u[ 1 ];
      const auto tmp6 = u[ 2 ] * tmp1;
      const auto tmp7 = (tmp0[ 0 ])[ 0 ] * u[ 2 ];
      (result[ 0 ])[ 0 ] = tmp2;
      (result[ 0 ])[ 1 ] = tmp3;
      (result[ 1 ])[ 0 ] = tmp4;
      (result[ 1 ])[ 1 ] = tmp5;
      (result[ 2 ])[ 0 ] = tmp6;
      (result[ 2 ])[ 1 ] = tmp7;
    }

    template< class Point >
    void linFlux ( const DRangeType &ubar, const DJacobianRangeType &dubar, const Point &x, const DRangeType &u, const DJacobianRangeType &du, RJacobianRangeType &result ) const
    {
      const auto tmp0 = [ this ] ( auto x ) {
          typename std::tuple_element_t< 0, CoefficientFunctionSpaceTupleType >::JacobianRangeType result;
          std::get< 0 >( coefficients_ ).jacobian( x, result );
          return result;
        }( (x) );
      const auto tmp1 = -1 * (tmp0[ 0 ])[ 1 ];
      const auto tmp2 = u[ 0 ] * tmp1;
      const auto tmp3 = u[ 0 ] * (tmp0[ 0 ])[ 0 ];
      const auto tmp4 = u[ 1 ] * tmp1;
      const auto tmp5 = (tmp0[ 0 ])[ 0 ] * u[ 1 ];
      const auto tmp6 = u[ 2 ] * tmp1;
      const auto tmp7 = (tmp0[ 0 ])[ 0 ] * u[ 2 ];
      (result[ 0 ])[ 0 ] = tmp2;
      (result[ 0 ])[ 1 ] = tmp3;
      (result[ 1 ])[ 0 ] = tmp4;
      (result[ 1 ])[ 1 ] = tmp5;
      (result[ 2 ])[ 0 ] = tmp6;
      (result[ 2 ])[ 1 ] = tmp7;
    }

    template< class Point >
    void fluxDivergence ( const Point &x, const DRangeType &u, const DJacobianRangeType &du, const DHessianRangeType &d2u, RRangeType &result ) const
    {
      const auto tmp0 = [ this ] ( auto x ) {
          typename std::tuple_element_t< 0, CoefficientFunctionSpaceTupleType >::HessianRangeType result;
          std::get< 0 >( coefficients_ ).hessian( x, result );
          return result;
        }( (x) );
      const auto tmp1 = -1 * ((tmp0[ 0 ])[ 1 ])[ 0 ];
      const auto tmp2 = u[ 0 ] * tmp1;
      const auto tmp3 = [ this ] ( auto x ) {
          typename std::tuple_element_t< 0, CoefficientFunctionSpaceTupleType >::JacobianRangeType result;
          std::get< 0 >( coefficients_ ).jacobian( x, result );
          return result;
        }( (x) );
      const auto tmp4 = -1 * (tmp3[ 0 ])[ 1 ];
      const auto tmp5 = (du[ 0 ])[ 0 ] * tmp4;
      const auto tmp6 = tmp5 + tmp2;
      const auto tmp7 = (tmp3[ 0 ])[ 0 ] * (du[ 0 ])[ 1 ];
      const auto tmp8 = u[ 0 ] * ((tmp0[ 0 ])[ 0 ])[ 1 ];
      const auto tmp9 = tmp8 + tmp7;
      const auto tmp10 = tmp9 + tmp6;
      const auto tmp11 = -1 * tmp10;
      auto tmp12 = entity().geometry().global( Dune::Fem::coordinate( x ) );
      const auto tmp13 = -0.6283185307179586 + tmp12[ 1 ];
      const auto tmp14 = tmp13 * tmp13;
      const auto tmp15 = -0.6283185307179586 + tmp12[ 0 ];
      const auto tmp16 = tmp15 * tmp15;
      const auto tmp17 = tmp16 + tmp14;
      double tmp18 = constant< 0 >();
      const auto tmp19 = u[ 0 ] * u[ 1 ];
      const auto tmp20 = 10 * tmp19;
      const auto tmp21 = -1 * tmp20;
      const auto tmp22 = tmp21 + (tmp18 < 5 ? (tmp17 < 0.2 ? 1 : 0.0) : 0.0);
      const auto tmp23 = tmp22 + tmp11;
      const auto tmp24 = u[ 1 ] * tmp1;
      const auto tmp25 = (du[ 1 ])[ 0 ] * tmp4;
      const auto tmp26 = tmp25 + tmp24;
      const auto tmp27 = (tmp3[ 0 ])[ 0 ] * (du[ 1 ])[ 1 ];
      const auto tmp28 = ((tmp0[ 0 ])[ 0 ])[ 1 ] * u[ 1 ];
      const auto tmp29 = tmp28 + tmp27;
      const auto tmp30 = tmp29 + tmp26;
      const auto tmp31 = -1 * tmp30;
      const auto tmp32 = -5.654866776461628 + tmp12[ 1 ];
      const auto tmp33 = tmp32 * tmp32;
      const auto tmp34 = -5.654866776461628 + tmp12[ 0 ];
      const auto tmp35 = tmp34 * tmp34;
      const auto tmp36 = tmp35 + tmp33;
      const auto tmp37 = tmp21 + (tmp18 < 5 ? (tmp36 < 0.2 ? 1 : 0.0) : 0.0);
      const auto tmp38 = tmp37 + tmp31;
      const auto tmp39 = -2 * u[ 0 ];
      const auto tmp40 = u[ 1 ] * tmp39;
      const auto tmp41 = 10 * tmp40;
      const auto tmp42 = -1 * tmp41;
      const auto tmp43 = u[ 2 ] * tmp1;
      const auto tmp44 = (du[ 2 ])[ 0 ] * tmp4;
      const auto tmp45 = tmp44 + tmp43;
      const auto tmp46 = (tmp3[ 0 ])[ 0 ] * (du[ 2 ])[ 1 ];
      const auto tmp47 = ((tmp0[ 0 ])[ 0 ])[ 1 ] * u[ 2 ];
      const auto tmp48 = tmp47 + tmp46;
      const auto tmp49 = tmp48 + tmp45;
      const auto tmp50 = -1 * tmp49;
      const auto tmp51 = tmp50 + tmp42;
      result[ 0 ] = tmp23;
      result[ 1 ] = tmp38;
      result[ 2 ] = tmp51;
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
    mutable std::tuple< Dune::Fem::ConstLocalFunction<StreamFunction>  > coefficients_;

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
          typename std::tuple_element_t< 0, CoefficientFunctionSpaceTupleType >::JacobianRangeType result;
          std::get< 0 >( coefficients_ ).jacobian( x, result );
          return result;
        }( (x) );
      const auto tmp1 = -1 * (tmp0[ 0 ])[ 1 ];
      const auto tmp2 = normal[ 0 ] * tmp1;
      const auto tmp3 = (tmp0[ 0 ])[ 0 ] * normal[ 1 ];
      const auto tmp4 = tmp3 + tmp2;
      const auto tmp5 = std::abs( tmp4 );
      result = tmp5;
      return result;
    }

    template< class Entity, class Point >
    inline DDomainType velocity ( const double &t, const Entity &entity, const Point &x, const DRangeType &u ) const
    {
      DDomainType result;
      const auto tmp0 = [ this ] ( auto x ) {
          typename std::tuple_element_t< 0, CoefficientFunctionSpaceTupleType >::JacobianRangeType result;
          std::get< 0 >( coefficients_ ).jacobian( x, result );
          return result;
        }( (x) );
      const auto tmp1 = -1 * (tmp0[ 0 ])[ 1 ];
      result[ 0 ] = tmp1;
      result[ 1 ] = (tmp0[ 0 ])[ 0 ];
      return result;
    }

    template< class Entity, class Point >
    inline double diffusionTimeStep ( const Entity& entity, const Point &x, const DRangeType &u ) const
    {
      double result;
      result = double( (0) );
      return result;
    }
    static constexpr bool hasPhysical = true;

    template< class Entity, class Point >
    inline double physical ( const Entity &entity, const Point &x, const DRangeType &u ) const
    {
      double result;
      const auto tmp0 = (u[ 0 ] >= 0.0 ? 1 : 0.0) * (u[ 1 ] >= 0.0 ? 1 : 0.0);
      const auto tmp1 = tmp0 * (u[ 2 ] >= 0.0 ? 1 : 0.0);
      result = tmp1;
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
    inline bool diffusionBoundaryFlux ( const int bndId, const double &t, const Entity& entity, const Point &x, const DDomainType &normal, const DRangeType &u, const DJacobianRangeType &jac, RRangeType &result ) const
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
          result[ 1 ] = 0.0;
          result[ 2 ] = 0.0;
          return true;
        }
        break;
      case 2:
        {
          result[ 0 ] = 0.0;
          result[ 1 ] = 0.0;
          result[ 2 ] = 0.0;
          return true;
        }
        break;
      case 3:
        {
          result[ 0 ] = 0.0;
          result[ 1 ] = 0.0;
          result[ 2 ] = 0.0;
          return true;
        }
        break;
      case 4:
        {
          result[ 0 ] = 0.0;
          result[ 1 ] = 0.0;
          result[ 2 ] = 0.0;
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
      globalMin = {std::numeric_limits<double>::min(), std::numeric_limits<double>::min(), std::numeric_limits<double>::min()};
      globalMax = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
    }
  };

} // namespace ModelImpl_nv1b55276a261f6fa5ab7b1cb7a3cc346bv1_3_19659fe2f2f74bcfb5d532f1a0d38be9

namespace ModelImpl_nv1b55276a261f6fa5ab7b1cb7a3cc346bv1_3_19659fe2f2f74bcfb5d532f1a0d38be9
{

  typedef ModelImpl_nv1b55276a261f6fa5ab7b1cb7a3cc346bv1_3_19659fe2f2f74bcfb5d532f1a0d38be9::Model< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction<typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::FieldVector< double, 1 > > > ModelType;

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

    return [ dispatch ] ( ModelImpl_nv1b55276a261f6fa5ab7b1cb7a3cc346bv1_3_19659fe2f2f74bcfb5d532f1a0d38be9::Model< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction<typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::FieldVector< double, 1 > > > &model, pybind11::handle coeff, pybind11::handle value ) {
        std::size_t k = renumberConstants( coeff );
        if( k >= dispatch.size() )
          throw std::range_error( "No such coefficient: "+std::to_string(k)+" >= "+std::to_string(dispatch.size()) );
        dispatch[ k ]( model, value );
        return k;
      };
  }

} // namespace ModelImpl_nv1b55276a261f6fa5ab7b1cb7a3cc346bv1_3_19659fe2f2f74bcfb5d532f1a0d38be9

PYBIND11_MODULE( modelFemDgAdv_nv1b55276a261f6fa5ab7b1cb7a3cc346bv1_3_19659fe2f2f74bcfb5d532f1a0d38be9, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView > GridPart;
  typedef ModelImpl_nv1b55276a261f6fa5ab7b1cb7a3cc346bv1_3_19659fe2f2f74bcfb5d532f1a0d38be9::Model< GridPart, Dune::FemPy::VirtualizedGridFunction<typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::FieldVector< double, 1 > > > Model;
  typedef Model ModelWrapper;
  auto cls = Dune::Python::insertClass<ModelWrapper>(module,"Model",Dune::Python::GenerateTypeName("ModelImpl_nv1b55276a261f6fa5ab7b1cb7a3cc346bv1_3_19659fe2f2f74bcfb5d532f1a0d38be9::Model< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction<typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::FieldVector< double, 1 > > >"),Dune::Python::IncludeFiles({"python/dune/generated/modelFemDgAdv_nv1b55276a261f6fa5ab7b1cb7a3cc346bv1_3_19659fe2f2f74bcfb5d532f1a0d38be9.cc"})).first;
  cls.def_property_readonly( "dimRange", [] ( ModelWrapper & ) { return 3; } );
  cls.def_property_readonly( "hasDirichletBoundary", [] ( ModelWrapper& ) -> bool { return false;});

  cls.def_property( "time", [] ( ModelWrapper &self ) { return self.template constant<0>(); }, [] ( ModelWrapper &self, typename ModelWrapper::ConstantType<0>& value) { self.template constant<0>() = value; });

  cls.def( "setConstant",ModelImpl_nv1b55276a261f6fa5ab7b1cb7a3cc346bv1_3_19659fe2f2f74bcfb5d532f1a0d38be9::defSetConstant( std::make_index_sequence< Model::numConstants >() ) );

  cls.def( pybind11::init( [] ( const Dune::FemPy::VirtualizedGridFunction< GridPart, Dune::FieldVector< double, 1 > > &coefficient0 ) {
    return new  ModelWrapper( coefficient0 );
    }), pybind11::keep_alive< 1, 2 >() );
}
#endif // GuardModelImpl_nv1b55276a261f6fa5ab7b1cb7a3cc346bv1_3_19659fe2f2f74bcfb5d532f1a0d38be9
