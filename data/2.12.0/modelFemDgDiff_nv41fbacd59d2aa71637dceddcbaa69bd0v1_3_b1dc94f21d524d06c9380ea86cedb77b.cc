#ifndef GuardModelImpl_nv41fbacd59d2aa71637dceddcbaa69bd0v1_3_b1dc94f21d524d06c9380ea86cedb77b
#define GuardModelImpl_nv41fbacd59d2aa71637dceddcbaa69bd0v1_3_b1dc94f21d524d06c9380ea86cedb77b
#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/python/alugrid/hierarchical.hh>
#include <dune/python/grid/hierarchical.hh>

#include <dune/fem/misc/boundaryidprovider.hh>

#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/extensions.h>

#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fem/schemes/conservationlawmodel.hh>

namespace ModelImpl_nv41fbacd59d2aa71637dceddcbaa69bd0v1_3_b1dc94f21d524d06c9380ea86cedb77b
{




  // Model
  // -----

  template< class GridPart >
  struct Model
  {
    typedef GridPart GridPartType;
    typedef typename GridPart::template Codim< 0 >::EntityType EntityType;
    typedef typename GridPart::IntersectionType IntersectionType;
    typedef Dune::Fem::FunctionSpace< typename GridPartType::ctype, double, GridPartType::dimensionworld, 4 > DFunctionSpaceType;
    typedef typename DFunctionSpaceType::DomainFieldType DDomainFieldType;
    typedef typename DFunctionSpaceType::RangeFieldType DRangeFieldType;
    typedef typename DFunctionSpaceType::DomainType DDomainType;
    typedef typename DFunctionSpaceType::RangeType DRangeType;
    typedef typename DFunctionSpaceType::JacobianRangeType DJacobianRangeType;
    typedef typename DFunctionSpaceType::HessianRangeType DHessianRangeType;
    static const int dimDomain = GridPartType::dimensionworld;
    static const int dimD = 4;
    typedef Dune::Fem::FunctionSpace< typename GridPartType::ctype, double, GridPartType::dimensionworld, 4 > RFunctionSpaceType;
    typedef typename RFunctionSpaceType::DomainFieldType RDomainFieldType;
    typedef typename RFunctionSpaceType::RangeFieldType RRangeFieldType;
    typedef typename RFunctionSpaceType::DomainType RDomainType;
    typedef typename RFunctionSpaceType::RangeType RRangeType;
    typedef typename RFunctionSpaceType::JacobianRangeType RJacobianRangeType;
    typedef typename RFunctionSpaceType::HessianRangeType RHessianRangeType;
    static const int dimR = 4;
    static const int dimLocal = GridPartType::dimension;
    template< std::size_t i >
    using ConstantType = typename std::tuple_element_t< i, std::tuple< std::shared_ptr< double > > >::element_type;
    static const std::size_t numConstants = 1;

    Model ( const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    {
      std::get< 0 >( constants_ ) = std::make_shared< double >();
    }

    bool init ( const EntityType &entity ) const
    {
      {
        entity_ = &entity;
      }
      return true;
    }

    void unbind () const
    {
      {
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
      result[ 1 ] = 0.0;
      result[ 2 ] = 0.0;
      result[ 3 ] = 0.0;
    }

    template< class Point >
    void linSource ( const DRangeType &ubar, const DJacobianRangeType &dubar, const Point &x, const DRangeType &u, const DJacobianRangeType &du, RRangeType &result ) const
    {
      result[ 0 ] = 0.0;
      result[ 1 ] = 0.0;
      result[ 2 ] = 0.0;
      result[ 3 ] = 0.0;
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
      const auto tmp7 = -1 * (du[ 1 ])[ 0 ];
        assert(tmp7 == tmp7);
      const auto tmp8 = (du[ 1 ])[ 0 ] + tmp7;
        assert(tmp8 == tmp8);
      const auto tmp9 = tmp8 * tmp0;
        assert(tmp9 == tmp9);
      const auto tmp10 = -1 * (du[ 1 ])[ 1 ];
        assert(tmp10 == tmp10);
      const auto tmp11 = (du[ 1 ])[ 1 ] + tmp10;
        assert(tmp11 == tmp11);
      const auto tmp12 = tmp11 * tmp0;
        assert(tmp12 == tmp12);
      const auto tmp13 = -1 * (du[ 2 ])[ 0 ];
        assert(tmp13 == tmp13);
      const auto tmp14 = (du[ 2 ])[ 0 ] + tmp13;
        assert(tmp14 == tmp14);
      const auto tmp15 = tmp14 * tmp0;
        assert(tmp15 == tmp15);
      const auto tmp16 = -1 * (du[ 2 ])[ 1 ];
        assert(tmp16 == tmp16);
      const auto tmp17 = (du[ 2 ])[ 1 ] + tmp16;
        assert(tmp17 == tmp17);
      const auto tmp18 = tmp17 * tmp0;
        assert(tmp18 == tmp18);
      const auto tmp19 = -1 * (du[ 3 ])[ 0 ];
        assert(tmp19 == tmp19);
      const auto tmp20 = (du[ 3 ])[ 0 ] + tmp19;
        assert(tmp20 == tmp20);
      const auto tmp21 = tmp20 * tmp0;
        assert(tmp21 == tmp21);
      const auto tmp22 = -1 * (du[ 3 ])[ 1 ];
        assert(tmp22 == tmp22);
      const auto tmp23 = (du[ 3 ])[ 1 ] + tmp22;
        assert(tmp23 == tmp23);
      const auto tmp24 = tmp23 * tmp0;
        assert(tmp24 == tmp24);
      (result[ 0 ])[ 0 ] = tmp3;
      (result[ 0 ])[ 1 ] = tmp6;
      (result[ 1 ])[ 0 ] = tmp9;
      (result[ 1 ])[ 1 ] = tmp12;
      (result[ 2 ])[ 0 ] = tmp15;
      (result[ 2 ])[ 1 ] = tmp18;
      (result[ 3 ])[ 0 ] = tmp21;
      (result[ 3 ])[ 1 ] = tmp24;
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
      const auto tmp7 = -1 * (du[ 1 ])[ 0 ];
        assert(tmp7 == tmp7);
      const auto tmp8 = (du[ 1 ])[ 0 ] + tmp7;
        assert(tmp8 == tmp8);
      const auto tmp9 = tmp8 * tmp0;
        assert(tmp9 == tmp9);
      const auto tmp10 = -1 * (du[ 1 ])[ 1 ];
        assert(tmp10 == tmp10);
      const auto tmp11 = (du[ 1 ])[ 1 ] + tmp10;
        assert(tmp11 == tmp11);
      const auto tmp12 = tmp11 * tmp0;
        assert(tmp12 == tmp12);
      const auto tmp13 = -1 * (du[ 2 ])[ 0 ];
        assert(tmp13 == tmp13);
      const auto tmp14 = (du[ 2 ])[ 0 ] + tmp13;
        assert(tmp14 == tmp14);
      const auto tmp15 = tmp14 * tmp0;
        assert(tmp15 == tmp15);
      const auto tmp16 = -1 * (du[ 2 ])[ 1 ];
        assert(tmp16 == tmp16);
      const auto tmp17 = (du[ 2 ])[ 1 ] + tmp16;
        assert(tmp17 == tmp17);
      const auto tmp18 = tmp17 * tmp0;
        assert(tmp18 == tmp18);
      const auto tmp19 = -1 * (du[ 3 ])[ 0 ];
        assert(tmp19 == tmp19);
      const auto tmp20 = (du[ 3 ])[ 0 ] + tmp19;
        assert(tmp20 == tmp20);
      const auto tmp21 = tmp20 * tmp0;
        assert(tmp21 == tmp21);
      const auto tmp22 = -1 * (du[ 3 ])[ 1 ];
        assert(tmp22 == tmp22);
      const auto tmp23 = (du[ 3 ])[ 1 ] + tmp22;
        assert(tmp23 == tmp23);
      const auto tmp24 = tmp23 * tmp0;
        assert(tmp24 == tmp24);
      (result[ 0 ])[ 0 ] = tmp3;
      (result[ 0 ])[ 1 ] = tmp6;
      (result[ 1 ])[ 0 ] = tmp9;
      (result[ 1 ])[ 1 ] = tmp12;
      (result[ 2 ])[ 0 ] = tmp15;
      (result[ 2 ])[ 1 ] = tmp18;
      (result[ 3 ])[ 0 ] = tmp21;
      (result[ 3 ])[ 1 ] = tmp24;
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
      const auto tmp9 = -1 * ((d2u[ 1 ])[ 1 ])[ 1 ];
        assert(tmp9 == tmp9);
      const auto tmp10 = ((d2u[ 1 ])[ 1 ])[ 1 ] + tmp9;
        assert(tmp10 == tmp10);
      const auto tmp11 = tmp10 * tmp0;
        assert(tmp11 == tmp11);
      const auto tmp12 = -1 * ((d2u[ 1 ])[ 0 ])[ 0 ];
        assert(tmp12 == tmp12);
      const auto tmp13 = ((d2u[ 1 ])[ 0 ])[ 0 ] + tmp12;
        assert(tmp13 == tmp13);
      const auto tmp14 = tmp13 * tmp0;
        assert(tmp14 == tmp14);
      const auto tmp15 = tmp14 + tmp11;
        assert(tmp15 == tmp15);
      const auto tmp16 = -1 * tmp15;
        assert(tmp16 == tmp16);
      const auto tmp17 = -1 * ((d2u[ 2 ])[ 1 ])[ 1 ];
        assert(tmp17 == tmp17);
      const auto tmp18 = ((d2u[ 2 ])[ 1 ])[ 1 ] + tmp17;
        assert(tmp18 == tmp18);
      const auto tmp19 = tmp18 * tmp0;
        assert(tmp19 == tmp19);
      const auto tmp20 = -1 * ((d2u[ 2 ])[ 0 ])[ 0 ];
        assert(tmp20 == tmp20);
      const auto tmp21 = ((d2u[ 2 ])[ 0 ])[ 0 ] + tmp20;
        assert(tmp21 == tmp21);
      const auto tmp22 = tmp21 * tmp0;
        assert(tmp22 == tmp22);
      const auto tmp23 = tmp22 + tmp19;
        assert(tmp23 == tmp23);
      const auto tmp24 = -1 * tmp23;
        assert(tmp24 == tmp24);
      const auto tmp25 = -1 * ((d2u[ 3 ])[ 1 ])[ 1 ];
        assert(tmp25 == tmp25);
      const auto tmp26 = ((d2u[ 3 ])[ 1 ])[ 1 ] + tmp25;
        assert(tmp26 == tmp26);
      const auto tmp27 = tmp26 * tmp0;
        assert(tmp27 == tmp27);
      const auto tmp28 = -1 * ((d2u[ 3 ])[ 0 ])[ 0 ];
        assert(tmp28 == tmp28);
      const auto tmp29 = ((d2u[ 3 ])[ 0 ])[ 0 ] + tmp28;
        assert(tmp29 == tmp29);
      const auto tmp30 = tmp29 * tmp0;
        assert(tmp30 == tmp30);
      const auto tmp31 = tmp30 + tmp27;
        assert(tmp31 == tmp31);
      const auto tmp32 = -1 * tmp31;
        assert(tmp32 == tmp32);
      result[ 0 ] = tmp8;
      result[ 1 ] = tmp16;
      result[ 2 ] = tmp24;
      result[ 3 ] = tmp32;
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

  public:
    static constexpr bool hasGamma = true;
    static constexpr double gamma = 1.4;

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
      using std::sqrt;
      const auto tmp0 = u[ 2 ] / u[ 0 ];
        assert(tmp0 == tmp0);
      const auto tmp1 = tmp0 * tmp0;
        assert(tmp1 == tmp1);
      const auto tmp2 = u[ 1 ] / u[ 0 ];
        assert(tmp2 == tmp2);
      const auto tmp3 = tmp2 * tmp2;
        assert(tmp3 == tmp3);
      const auto tmp4 = tmp3 + tmp1;
        assert(tmp4 == tmp4);
      const auto tmp5 = u[ 0 ] * tmp4;
        assert(tmp5 == tmp5);
      const auto tmp6 = tmp5 / 2;
        assert(tmp6 == tmp6);
      const auto tmp7 = -1 * tmp6;
        assert(tmp7 == tmp7);
      const auto tmp8 = u[ 3 ] + tmp7;
        assert(tmp8 == tmp8);
      const auto tmp9 = 0.3999999999999999 * tmp8;
        assert(tmp9 == tmp9);
      const auto tmp10 = 1.4 * tmp9;
        assert(tmp10 == tmp10);
      const auto tmp11 = tmp10 / u[ 0 ];
        assert(tmp11 == tmp11);
      const auto tmp12 = std::sqrt( tmp11 );
        assert(tmp12 == tmp12);
      const auto tmp13 = normal[ 1 ] * tmp0;
        assert(tmp13 == tmp13);
      const auto tmp14 = normal[ 0 ] * tmp2;
        assert(tmp14 == tmp14);
      const auto tmp15 = tmp14 + tmp13;
        assert(tmp15 == tmp15);
      const auto tmp16 = std::abs( tmp15 );
        assert(tmp16 == tmp16);
      const auto tmp17 = tmp16 + tmp12;
        assert(tmp17 == tmp17);
      result = tmp17;
      return result;
    }

    template< class Entity, class Point >
    inline DDomainType velocity ( const double &t, const Entity &entity, const Point &x, const DRangeType &u ) const
    {
      DDomainType result;
      const auto tmp0 = u[ 1 ] / u[ 0 ];
        assert(tmp0 == tmp0);
      const auto tmp1 = u[ 2 ] / u[ 0 ];
        assert(tmp1 == tmp1);
      result[ 0 ] = tmp0;
      result[ 1 ] = tmp1;
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
      const auto tmp0 = u[ 2 ] / u[ 0 ];
        assert(tmp0 == tmp0);
      const auto tmp1 = tmp0 * tmp0;
        assert(tmp1 == tmp1);
      const auto tmp2 = u[ 1 ] / u[ 0 ];
        assert(tmp2 == tmp2);
      const auto tmp3 = tmp2 * tmp2;
        assert(tmp3 == tmp3);
      const auto tmp4 = tmp3 + tmp1;
        assert(tmp4 == tmp4);
      const auto tmp5 = u[ 0 ] * tmp4;
        assert(tmp5 == tmp5);
      const auto tmp6 = tmp5 / 2;
        assert(tmp6 == tmp6);
      const auto tmp7 = -1 * tmp6;
        assert(tmp7 == tmp7);
      const auto tmp8 = u[ 3 ] + tmp7;
        assert(tmp8 == tmp8);
      const auto tmp9 = 0.3999999999999999 * tmp8;
        assert(tmp9 == tmp9);
      result = (u[ 0 ] > 1e-08 ? (tmp9 > 1e-08 ? 1 : 0.0) : 0.0);
      return result;
    }

    template< class Intersection, class Point >
    inline double jump ( const Intersection& it, const Point &x, const DRangeType &u, const DRangeType &w ) const
    {
      double result;
      const auto tmp0 = w[ 2 ] / w[ 0 ];
        assert(tmp0 == tmp0);
      const auto tmp1 = tmp0 * tmp0;
        assert(tmp1 == tmp1);
      const auto tmp2 = w[ 1 ] / w[ 0 ];
        assert(tmp2 == tmp2);
      const auto tmp3 = tmp2 * tmp2;
        assert(tmp3 == tmp3);
      const auto tmp4 = tmp3 + tmp1;
        assert(tmp4 == tmp4);
      const auto tmp5 = w[ 0 ] * tmp4;
        assert(tmp5 == tmp5);
      const auto tmp6 = tmp5 / 2;
        assert(tmp6 == tmp6);
      const auto tmp7 = -1 * tmp6;
        assert(tmp7 == tmp7);
      const auto tmp8 = w[ 3 ] + tmp7;
        assert(tmp8 == tmp8);
      const auto tmp9 = 0.3999999999999999 * tmp8;
        assert(tmp9 == tmp9);
      const auto tmp10 = u[ 2 ] / u[ 0 ];
        assert(tmp10 == tmp10);
      const auto tmp11 = tmp10 * tmp10;
        assert(tmp11 == tmp11);
      const auto tmp12 = u[ 1 ] / u[ 0 ];
        assert(tmp12 == tmp12);
      const auto tmp13 = tmp12 * tmp12;
        assert(tmp13 == tmp13);
      const auto tmp14 = tmp13 + tmp11;
        assert(tmp14 == tmp14);
      const auto tmp15 = u[ 0 ] * tmp14;
        assert(tmp15 == tmp15);
      const auto tmp16 = tmp15 / 2;
        assert(tmp16 == tmp16);
      const auto tmp17 = -1 * tmp16;
        assert(tmp17 == tmp17);
      const auto tmp18 = u[ 3 ] + tmp17;
        assert(tmp18 == tmp18);
      const auto tmp19 = 0.3999999999999999 * tmp18;
        assert(tmp19 == tmp19);
      const auto tmp20 = tmp19 + tmp9;
        assert(tmp20 == tmp20);
      const auto tmp21 = 0.5 * tmp20;
        assert(tmp21 == tmp21);
      const auto tmp22 = -1 * tmp9;
        assert(tmp22 == tmp22);
      const auto tmp23 = tmp19 + tmp22;
        assert(tmp23 == tmp23);
      const auto tmp24 = tmp23 / tmp21;
        assert(tmp24 == tmp24);
      result = tmp24;
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
          result[ 0 ] = u[ 0 ];
          result[ 1 ] = u[ 1 ];
          result[ 2 ] = u[ 2 ];
          result[ 3 ] = u[ 3 ];
          return true;
        }
        break;
      case 2:
        {
          result[ 0 ] = u[ 0 ];
          result[ 1 ] = u[ 1 ];
          result[ 2 ] = u[ 2 ];
          result[ 3 ] = u[ 3 ];
          return true;
        }
        break;
      case 3:
        {
          result[ 0 ] = u[ 0 ];
          result[ 1 ] = u[ 1 ];
          result[ 2 ] = u[ 2 ];
          result[ 3 ] = u[ 3 ];
          return true;
        }
        break;
      case 4:
        {
          result[ 0 ] = u[ 0 ];
          result[ 1 ] = u[ 1 ];
          result[ 2 ] = u[ 2 ];
          result[ 3 ] = u[ 3 ];
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
      globalMin = {std::numeric_limits<double>::min(), std::numeric_limits<double>::min(), std::numeric_limits<double>::min(), std::numeric_limits<double>::min()};
      globalMax = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
    }
  };

} // namespace ModelImpl_nv41fbacd59d2aa71637dceddcbaa69bd0v1_3_b1dc94f21d524d06c9380ea86cedb77b

namespace ModelImpl_nv41fbacd59d2aa71637dceddcbaa69bd0v1_3_b1dc94f21d524d06c9380ea86cedb77b
{

  typedef ModelImpl_nv41fbacd59d2aa71637dceddcbaa69bd0v1_3_b1dc94f21d524d06c9380ea86cedb77b::Model< typename Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView > > ModelType;

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

    return [ dispatch ] ( ModelImpl_nv41fbacd59d2aa71637dceddcbaa69bd0v1_3_b1dc94f21d524d06c9380ea86cedb77b::Model< typename Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView > > &model, pybind11::handle coeff, pybind11::handle value ) {
        std::size_t k = renumberConstants( coeff );
        if( k >= dispatch.size() )
          throw std::range_error( "No such coefficient: "+std::to_string(k)+" >= "+std::to_string(dispatch.size()) );
        dispatch[ k ]( model, value );
        return k;
      };
  }

} // namespace ModelImpl_nv41fbacd59d2aa71637dceddcbaa69bd0v1_3_b1dc94f21d524d06c9380ea86cedb77b

PYBIND11_MODULE( modelFemDgDiff_nv41fbacd59d2aa71637dceddcbaa69bd0v1_3_b1dc94f21d524d06c9380ea86cedb77b, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView > GridPart;
  typedef ModelImpl_nv41fbacd59d2aa71637dceddcbaa69bd0v1_3_b1dc94f21d524d06c9380ea86cedb77b::Model< GridPart > Model;
  typedef Model ModelWrapper;
  auto cls = Dune::Python::insertClass<ModelWrapper>(module,"Model",Dune::Python::GenerateTypeName("ModelImpl_nv41fbacd59d2aa71637dceddcbaa69bd0v1_3_b1dc94f21d524d06c9380ea86cedb77b::Model< typename Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView > >"),Dune::Python::IncludeFiles({"python/dune/generated/modelFemDgDiff_nv41fbacd59d2aa71637dceddcbaa69bd0v1_3_b1dc94f21d524d06c9380ea86cedb77b.cc"})).first;
  cls.def_property_readonly( "dimRange", [] ( ModelWrapper & ) { return 4; } );
  cls.def_property_readonly( "hasDirichletBoundary", [] ( ModelWrapper& ) -> bool { return false;});

  cls.def_property( "time", [] ( ModelWrapper &self ) { return self.template constant<0>(); }, [] ( ModelWrapper &self, typename ModelWrapper::ConstantType<0>& value) { self.template constant<0>() = value; });

  cls.def( "setConstant",ModelImpl_nv41fbacd59d2aa71637dceddcbaa69bd0v1_3_b1dc94f21d524d06c9380ea86cedb77b::defSetConstant( std::make_index_sequence< Model::numConstants >() ) );

  cls.def( pybind11::init( [] (  ) {
    return new  ModelWrapper();
    } ) );
}
#endif // GuardModelImpl_nv41fbacd59d2aa71637dceddcbaa69bd0v1_3_b1dc94f21d524d06c9380ea86cedb77b
