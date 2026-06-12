#ifndef GuardModelImpl_nv62309d681344dcb62d9fe22c916fb8d2v1_3_19659fe2f2f74bcfb5d532f1a0d38be9
#define GuardModelImpl_nv62309d681344dcb62d9fe22c916fb8d2v1_3_19659fe2f2f74bcfb5d532f1a0d38be9
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

namespace ModelImpl_nv62309d681344dcb62d9fe22c916fb8d2v1_3_19659fe2f2f74bcfb5d532f1a0d38be9
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
    using CoefficientType = std::tuple_element_t< i, std::tuple< Dune::Fem::ConstLocalFunction<StreamFunction>  > >;
    typedef std::tuple< typename CoefficientType<0>::FunctionSpaceType > CoefficientFunctionSpaceTupleType;

    Model ( const StreamFunction &coefficient0, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
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
      result[ 1 ] = 0.0;
      result[ 2 ] = 0.0;
    }

    template< class Point >
    void linSource ( const DRangeType &ubar, const DJacobianRangeType &dubar, const Point &x, const DRangeType &u, const DJacobianRangeType &du, RRangeType &result ) const
    {
      result[ 0 ] = 0.0;
      result[ 1 ] = 0.0;
      result[ 2 ] = 0.0;
    }

    template< class Point >
    void flux ( const Point &x, const DRangeType &u, const DJacobianRangeType &du, RJacobianRangeType &result ) const
    {
      const auto tmp0 = 0.02 * (du[ 0 ])[ 0 ];
        assert(tmp0 == tmp0);
      const auto tmp1 = 0.02 * (du[ 0 ])[ 1 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = 0.02 * (du[ 1 ])[ 0 ];
        assert(tmp2 == tmp2);
      const auto tmp3 = 0.02 * (du[ 1 ])[ 1 ];
        assert(tmp3 == tmp3);
      const auto tmp4 = 0.02 * (du[ 2 ])[ 0 ];
        assert(tmp4 == tmp4);
      const auto tmp5 = 0.02 * (du[ 2 ])[ 1 ];
        assert(tmp5 == tmp5);
      (result[ 0 ])[ 0 ] = tmp0;
      (result[ 0 ])[ 1 ] = tmp1;
      (result[ 1 ])[ 0 ] = tmp2;
      (result[ 1 ])[ 1 ] = tmp3;
      (result[ 2 ])[ 0 ] = tmp4;
      (result[ 2 ])[ 1 ] = tmp5;
    }

    template< class Point >
    void linFlux ( const DRangeType &ubar, const DJacobianRangeType &dubar, const Point &x, const DRangeType &u, const DJacobianRangeType &du, RJacobianRangeType &result ) const
    {
      const auto tmp0 = 0.02 * (du[ 0 ])[ 0 ];
        assert(tmp0 == tmp0);
      const auto tmp1 = 0.02 * (du[ 0 ])[ 1 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = 0.02 * (du[ 1 ])[ 0 ];
        assert(tmp2 == tmp2);
      const auto tmp3 = 0.02 * (du[ 1 ])[ 1 ];
        assert(tmp3 == tmp3);
      const auto tmp4 = 0.02 * (du[ 2 ])[ 0 ];
        assert(tmp4 == tmp4);
      const auto tmp5 = 0.02 * (du[ 2 ])[ 1 ];
        assert(tmp5 == tmp5);
      (result[ 0 ])[ 0 ] = tmp0;
      (result[ 0 ])[ 1 ] = tmp1;
      (result[ 1 ])[ 0 ] = tmp2;
      (result[ 1 ])[ 1 ] = tmp3;
      (result[ 2 ])[ 0 ] = tmp4;
      (result[ 2 ])[ 1 ] = tmp5;
    }

    template< class Point >
    void fluxDivergence ( const Point &x, const DRangeType &u, const DJacobianRangeType &du, const DHessianRangeType &d2u, RRangeType &result ) const
    {
      const auto tmp0 = 0.02 * ((d2u[ 0 ])[ 1 ])[ 1 ];
        assert(tmp0 == tmp0);
      const auto tmp1 = 0.02 * ((d2u[ 0 ])[ 0 ])[ 0 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = tmp1 + tmp0;
        assert(tmp2 == tmp2);
      const auto tmp3 = -1 * tmp2;
        assert(tmp3 == tmp3);
      const auto tmp4 = 0.02 * ((d2u[ 1 ])[ 1 ])[ 1 ];
        assert(tmp4 == tmp4);
      const auto tmp5 = 0.02 * ((d2u[ 1 ])[ 0 ])[ 0 ];
        assert(tmp5 == tmp5);
      const auto tmp6 = tmp5 + tmp4;
        assert(tmp6 == tmp6);
      const auto tmp7 = -1 * tmp6;
        assert(tmp7 == tmp7);
      const auto tmp8 = 0.02 * ((d2u[ 2 ])[ 1 ])[ 1 ];
        assert(tmp8 == tmp8);
      const auto tmp9 = 0.02 * ((d2u[ 2 ])[ 0 ])[ 0 ];
        assert(tmp9 == tmp9);
      const auto tmp10 = tmp9 + tmp8;
        assert(tmp10 == tmp10);
      const auto tmp11 = -1 * tmp10;
        assert(tmp11 == tmp11);
      result[ 0 ] = tmp3;
      result[ 1 ] = tmp7;
      result[ 2 ] = tmp11;
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
        assert(tmp0 == tmp0);
      const auto tmp1 = -1 * (tmp0[ 0 ])[ 1 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = normal[ 0 ] * tmp1;
        assert(tmp2 == tmp2);
      const auto tmp3 = (tmp0[ 0 ])[ 0 ] * normal[ 1 ];
        assert(tmp3 == tmp3);
      const auto tmp4 = tmp3 + tmp2;
        assert(tmp4 == tmp4);
      const auto tmp5 = std::abs( tmp4 );
        assert(tmp5 == tmp5);
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
        assert(tmp0 == tmp0);
      const auto tmp1 = -1 * (tmp0[ 0 ])[ 1 ];
        assert(tmp1 == tmp1);
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
        assert(tmp0 == tmp0);
      const auto tmp1 = tmp0 * (u[ 2 ] >= 0.0 ? 1 : 0.0);
        assert(tmp1 == tmp1);
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

} // namespace ModelImpl_nv62309d681344dcb62d9fe22c916fb8d2v1_3_19659fe2f2f74bcfb5d532f1a0d38be9

namespace ModelImpl_nv62309d681344dcb62d9fe22c916fb8d2v1_3_19659fe2f2f74bcfb5d532f1a0d38be9
{

} // namespace ModelImpl_nv62309d681344dcb62d9fe22c916fb8d2v1_3_19659fe2f2f74bcfb5d532f1a0d38be9

PYBIND11_MODULE( modelFemDgDiff_nv62309d681344dcb62d9fe22c916fb8d2v1_3_19659fe2f2f74bcfb5d532f1a0d38be9, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView > GridPart;
  typedef ModelImpl_nv62309d681344dcb62d9fe22c916fb8d2v1_3_19659fe2f2f74bcfb5d532f1a0d38be9::Model< GridPart, Dune::FemPy::VirtualizedGridFunction<typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::FieldVector< double, 1 > > > Model;
  typedef Model ModelWrapper;
  auto cls = Dune::Python::insertClass<ModelWrapper>(module,"Model",Dune::Python::GenerateTypeName("ModelImpl_nv62309d681344dcb62d9fe22c916fb8d2v1_3_19659fe2f2f74bcfb5d532f1a0d38be9::Model< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::FemPy::VirtualizedGridFunction<typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView >, Dune::FieldVector< double, 1 > > >"),Dune::Python::IncludeFiles({"python/dune/generated/modelFemDgDiff_nv62309d681344dcb62d9fe22c916fb8d2v1_3_19659fe2f2f74bcfb5d532f1a0d38be9.cc"})).first;
  cls.def_property_readonly( "dimRange", [] ( ModelWrapper & ) { return 3; } );
  cls.def_property_readonly( "hasDirichletBoundary", [] ( ModelWrapper& ) -> bool { return false;});



  cls.def( pybind11::init( [] ( const Dune::FemPy::VirtualizedGridFunction< GridPart, Dune::FieldVector< double, 1 > > &coefficient0 ) {
    return new  ModelWrapper( coefficient0 );
    }), pybind11::keep_alive< 1, 2 >() );
}
#endif // GuardModelImpl_nv62309d681344dcb62d9fe22c916fb8d2v1_3_19659fe2f2f74bcfb5d532f1a0d38be9
