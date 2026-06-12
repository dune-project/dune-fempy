#ifndef GuardModelImpl_719110241111e05224e2f1c3078e038fv1_3_f084d41ce9ff164524e5c3ea6f8082cc
#define GuardModelImpl_719110241111e05224e2f1c3078e038fv1_3_f084d41ce9ff164524e5c3ea6f8082cc
#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/python/grid/hierarchical.hh>
#include <dune/vem/misc/grid.hh>

#include <dune/fem/misc/boundaryidprovider.hh>

#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/extensions.h>

#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fem/schemes/conservationlawmodel.hh>

namespace ModelImpl_719110241111e05224e2f1c3078e038fv1_3_f084d41ce9ff164524e5c3ea6f8082cc
{




  // Model
  // -----

  template< class GridPart >
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

    Model ( const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    {}

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
      using std::cos;
      using std::sin;
      auto tmp0 = entity().geometry().global( Dune::Fem::coordinate( x ) );
        assert(tmp0 == tmp0);
      const auto tmp1 = tmp0[ 1 ] * tmp0[ 1 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = tmp0[ 0 ] * tmp0[ 0 ];
        assert(tmp2 == tmp2);
      const auto tmp3 = tmp2 + tmp1;
        assert(tmp3 == tmp3);
      const auto tmp4 = std::sin( tmp3 );
        assert(tmp4 == tmp4);
      const auto tmp5 = 1 + tmp4;
        assert(tmp5 == tmp5);
      const auto tmp6 = u[ 0 ] * tmp5;
        assert(tmp6 == tmp6);
      const auto tmp7 = 3.141592653589793 * tmp0[ 0 ];
        assert(tmp7 == tmp7);
      const auto tmp8 = tmp0[ 1 ] * tmp7;
        assert(tmp8 == tmp8);
      const auto tmp9 = std::cos( tmp8 );
        assert(tmp9 == tmp9);
      const auto tmp10 = tmp0[ 0 ] * tmp0[ 1 ];
        assert(tmp10 == tmp10);
      const auto tmp11 = tmp10 * tmp9;
        assert(tmp11 == tmp11);
      const auto tmp12 = tmp5 * tmp11;
        assert(tmp12 == tmp12);
      const auto tmp13 = -1 * tmp4;
        assert(tmp13 == tmp13);
      const auto tmp14 = tmp0[ 1 ] + tmp0[ 1 ];
        assert(tmp14 == tmp14);
      const auto tmp15 = tmp14 * tmp13;
        assert(tmp15 == tmp15);
      const auto tmp16 = 0.9 * tmp15;
        assert(tmp16 == tmp16);
      const auto tmp17 = -1 * tmp16;
        assert(tmp17 == tmp17);
      const auto tmp18 = tmp0[ 0 ] * tmp9;
        assert(tmp18 == tmp18);
      const auto tmp19 = std::sin( tmp8 );
        assert(tmp19 == tmp19);
      const auto tmp20 = -1 * tmp19;
        assert(tmp20 == tmp20);
      const auto tmp21 = tmp7 * tmp20;
        assert(tmp21 == tmp21);
      const auto tmp22 = tmp10 * tmp21;
        assert(tmp22 == tmp22);
      const auto tmp23 = tmp22 + tmp18;
        assert(tmp23 == tmp23);
      const auto tmp24 = tmp23 * tmp17;
        assert(tmp24 == tmp24);
      const auto tmp25 = std::cos( tmp3 );
        assert(tmp25 == tmp25);
      const auto tmp26 = 0.9 * tmp25;
        assert(tmp26 == tmp26);
      const auto tmp27 = -1 * tmp26;
        assert(tmp27 == tmp27);
      const auto tmp28 = 1 + tmp27;
        assert(tmp28 == tmp28);
      const auto tmp29 = tmp0[ 0 ] * tmp21;
        assert(tmp29 == tmp29);
      const auto tmp30 = tmp7 * tmp9;
        assert(tmp30 == tmp30);
      const auto tmp31 = -1 * tmp30;
        assert(tmp31 == tmp31);
      const auto tmp32 = tmp7 * tmp31;
        assert(tmp32 == tmp32);
      const auto tmp33 = tmp10 * tmp32;
        assert(tmp33 == tmp33);
      const auto tmp34 = tmp33 + tmp29;
        assert(tmp34 == tmp34);
      const auto tmp35 = tmp34 + tmp29;
        assert(tmp35 == tmp35);
      const auto tmp36 = tmp35 * tmp28;
        assert(tmp36 == tmp36);
      const auto tmp37 = tmp36 + tmp24;
        assert(tmp37 == tmp37);
      const auto tmp38 = tmp0[ 0 ] + tmp0[ 0 ];
        assert(tmp38 == tmp38);
      const auto tmp39 = tmp38 * tmp13;
        assert(tmp39 == tmp39);
      const auto tmp40 = 0.9 * tmp39;
        assert(tmp40 == tmp40);
      const auto tmp41 = -1 * tmp40;
        assert(tmp41 == tmp41);
      const auto tmp42 = tmp0[ 1 ] * tmp9;
        assert(tmp42 == tmp42);
      const auto tmp43 = 3.141592653589793 * tmp0[ 1 ];
        assert(tmp43 == tmp43);
      const auto tmp44 = tmp43 * tmp20;
        assert(tmp44 == tmp44);
      const auto tmp45 = tmp10 * tmp44;
        assert(tmp45 == tmp45);
      const auto tmp46 = tmp45 + tmp42;
        assert(tmp46 == tmp46);
      const auto tmp47 = tmp46 * tmp41;
        assert(tmp47 == tmp47);
      const auto tmp48 = tmp0[ 1 ] * tmp44;
        assert(tmp48 == tmp48);
      const auto tmp49 = tmp43 * tmp9;
        assert(tmp49 == tmp49);
      const auto tmp50 = -1 * tmp49;
        assert(tmp50 == tmp50);
      const auto tmp51 = tmp43 * tmp50;
        assert(tmp51 == tmp51);
      const auto tmp52 = tmp10 * tmp51;
        assert(tmp52 == tmp52);
      const auto tmp53 = tmp52 + tmp48;
        assert(tmp53 == tmp53);
      const auto tmp54 = tmp53 + tmp48;
        assert(tmp54 == tmp54);
      const auto tmp55 = tmp54 * tmp28;
        assert(tmp55 == tmp55);
      const auto tmp56 = tmp55 + tmp47;
        assert(tmp56 == tmp56);
      const auto tmp57 = tmp56 + tmp37;
        assert(tmp57 == tmp57);
      const auto tmp58 = -1 * tmp57;
        assert(tmp58 == tmp58);
      const auto tmp59 = tmp58 + tmp12;
        assert(tmp59 == tmp59);
      const auto tmp60 = -1 * tmp59;
        assert(tmp60 == tmp60);
      const auto tmp61 = tmp60 + tmp6;
        assert(tmp61 == tmp61);
      result[ 0 ] = tmp61;
    }

    template< class Point >
    void linSource ( const DRangeType &ubar, const DJacobianRangeType &dubar, const Point &x, const DRangeType &u, const DJacobianRangeType &du, RRangeType &result ) const
    {
      using std::sin;
      auto tmp0 = entity().geometry().global( Dune::Fem::coordinate( x ) );
        assert(tmp0 == tmp0);
      const auto tmp1 = tmp0[ 1 ] * tmp0[ 1 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = tmp0[ 0 ] * tmp0[ 0 ];
        assert(tmp2 == tmp2);
      const auto tmp3 = tmp2 + tmp1;
        assert(tmp3 == tmp3);
      const auto tmp4 = std::sin( tmp3 );
        assert(tmp4 == tmp4);
      const auto tmp5 = 1 + tmp4;
        assert(tmp5 == tmp5);
      const auto tmp6 = u[ 0 ] * tmp5;
        assert(tmp6 == tmp6);
      result[ 0 ] = tmp6;
    }

    template< class Point >
    void flux ( const Point &x, const DRangeType &u, const DJacobianRangeType &du, RJacobianRangeType &result ) const
    {
      using std::cos;
      auto tmp0 = entity().geometry().global( Dune::Fem::coordinate( x ) );
        assert(tmp0 == tmp0);
      const auto tmp1 = tmp0[ 1 ] * tmp0[ 1 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = tmp0[ 0 ] * tmp0[ 0 ];
        assert(tmp2 == tmp2);
      const auto tmp3 = tmp2 + tmp1;
        assert(tmp3 == tmp3);
      const auto tmp4 = std::cos( tmp3 );
        assert(tmp4 == tmp4);
      const auto tmp5 = 0.9 * tmp4;
        assert(tmp5 == tmp5);
      const auto tmp6 = -1 * tmp5;
        assert(tmp6 == tmp6);
      const auto tmp7 = 1 + tmp6;
        assert(tmp7 == tmp7);
      const auto tmp8 = (du[ 0 ])[ 0 ] * tmp7;
        assert(tmp8 == tmp8);
      const auto tmp9 = (du[ 0 ])[ 1 ] * tmp7;
        assert(tmp9 == tmp9);
      (result[ 0 ])[ 0 ] = tmp8;
      (result[ 0 ])[ 1 ] = tmp9;
    }

    template< class Point >
    void linFlux ( const DRangeType &ubar, const DJacobianRangeType &dubar, const Point &x, const DRangeType &u, const DJacobianRangeType &du, RJacobianRangeType &result ) const
    {
      using std::cos;
      auto tmp0 = entity().geometry().global( Dune::Fem::coordinate( x ) );
        assert(tmp0 == tmp0);
      const auto tmp1 = tmp0[ 1 ] * tmp0[ 1 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = tmp0[ 0 ] * tmp0[ 0 ];
        assert(tmp2 == tmp2);
      const auto tmp3 = tmp2 + tmp1;
        assert(tmp3 == tmp3);
      const auto tmp4 = std::cos( tmp3 );
        assert(tmp4 == tmp4);
      const auto tmp5 = 0.9 * tmp4;
        assert(tmp5 == tmp5);
      const auto tmp6 = -1 * tmp5;
        assert(tmp6 == tmp6);
      const auto tmp7 = 1 + tmp6;
        assert(tmp7 == tmp7);
      const auto tmp8 = (du[ 0 ])[ 0 ] * tmp7;
        assert(tmp8 == tmp8);
      const auto tmp9 = (du[ 0 ])[ 1 ] * tmp7;
        assert(tmp9 == tmp9);
      (result[ 0 ])[ 0 ] = tmp8;
      (result[ 0 ])[ 1 ] = tmp9;
    }

    template< class Point >
    void fluxDivergence ( const Point &x, const DRangeType &u, const DJacobianRangeType &du, const DHessianRangeType &d2u, RRangeType &result ) const
    {
      using std::cos;
      using std::sin;
      auto tmp0 = entity().geometry().global( Dune::Fem::coordinate( x ) );
        assert(tmp0 == tmp0);
      const auto tmp1 = tmp0[ 1 ] * tmp0[ 1 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = tmp0[ 0 ] * tmp0[ 0 ];
        assert(tmp2 == tmp2);
      const auto tmp3 = tmp2 + tmp1;
        assert(tmp3 == tmp3);
      const auto tmp4 = std::sin( tmp3 );
        assert(tmp4 == tmp4);
      const auto tmp5 = -1 * tmp4;
        assert(tmp5 == tmp5);
      const auto tmp6 = tmp0[ 1 ] + tmp0[ 1 ];
        assert(tmp6 == tmp6);
      const auto tmp7 = tmp6 * tmp5;
        assert(tmp7 == tmp7);
      const auto tmp8 = 0.9 * tmp7;
        assert(tmp8 == tmp8);
      const auto tmp9 = -1 * tmp8;
        assert(tmp9 == tmp9);
      const auto tmp10 = (du[ 0 ])[ 1 ] * tmp9;
        assert(tmp10 == tmp10);
      const auto tmp11 = std::cos( tmp3 );
        assert(tmp11 == tmp11);
      const auto tmp12 = 0.9 * tmp11;
        assert(tmp12 == tmp12);
      const auto tmp13 = -1 * tmp12;
        assert(tmp13 == tmp13);
      const auto tmp14 = 1 + tmp13;
        assert(tmp14 == tmp14);
      const auto tmp15 = ((d2u[ 0 ])[ 1 ])[ 1 ] * tmp14;
        assert(tmp15 == tmp15);
      const auto tmp16 = tmp15 + tmp10;
        assert(tmp16 == tmp16);
      const auto tmp17 = tmp0[ 0 ] + tmp0[ 0 ];
        assert(tmp17 == tmp17);
      const auto tmp18 = tmp17 * tmp5;
        assert(tmp18 == tmp18);
      const auto tmp19 = 0.9 * tmp18;
        assert(tmp19 == tmp19);
      const auto tmp20 = -1 * tmp19;
        assert(tmp20 == tmp20);
      const auto tmp21 = (du[ 0 ])[ 0 ] * tmp20;
        assert(tmp21 == tmp21);
      const auto tmp22 = ((d2u[ 0 ])[ 0 ])[ 0 ] * tmp14;
        assert(tmp22 == tmp22);
      const auto tmp23 = tmp22 + tmp21;
        assert(tmp23 == tmp23);
      const auto tmp24 = tmp23 + tmp16;
        assert(tmp24 == tmp24);
      const auto tmp25 = -1 * tmp24;
        assert(tmp25 == tmp25);
      const auto tmp26 = 1 + tmp4;
        assert(tmp26 == tmp26);
      const auto tmp27 = u[ 0 ] * tmp26;
        assert(tmp27 == tmp27);
      const auto tmp28 = 3.141592653589793 * tmp0[ 0 ];
        assert(tmp28 == tmp28);
      const auto tmp29 = tmp0[ 1 ] * tmp28;
        assert(tmp29 == tmp29);
      const auto tmp30 = std::cos( tmp29 );
        assert(tmp30 == tmp30);
      const auto tmp31 = tmp0[ 0 ] * tmp0[ 1 ];
        assert(tmp31 == tmp31);
      const auto tmp32 = tmp31 * tmp30;
        assert(tmp32 == tmp32);
      const auto tmp33 = tmp26 * tmp32;
        assert(tmp33 == tmp33);
      const auto tmp34 = tmp0[ 0 ] * tmp30;
        assert(tmp34 == tmp34);
      const auto tmp35 = std::sin( tmp29 );
        assert(tmp35 == tmp35);
      const auto tmp36 = -1 * tmp35;
        assert(tmp36 == tmp36);
      const auto tmp37 = tmp28 * tmp36;
        assert(tmp37 == tmp37);
      const auto tmp38 = tmp31 * tmp37;
        assert(tmp38 == tmp38);
      const auto tmp39 = tmp38 + tmp34;
        assert(tmp39 == tmp39);
      const auto tmp40 = tmp39 * tmp9;
        assert(tmp40 == tmp40);
      const auto tmp41 = tmp0[ 0 ] * tmp37;
        assert(tmp41 == tmp41);
      const auto tmp42 = tmp28 * tmp30;
        assert(tmp42 == tmp42);
      const auto tmp43 = -1 * tmp42;
        assert(tmp43 == tmp43);
      const auto tmp44 = tmp28 * tmp43;
        assert(tmp44 == tmp44);
      const auto tmp45 = tmp31 * tmp44;
        assert(tmp45 == tmp45);
      const auto tmp46 = tmp45 + tmp41;
        assert(tmp46 == tmp46);
      const auto tmp47 = tmp46 + tmp41;
        assert(tmp47 == tmp47);
      const auto tmp48 = tmp47 * tmp14;
        assert(tmp48 == tmp48);
      const auto tmp49 = tmp48 + tmp40;
        assert(tmp49 == tmp49);
      const auto tmp50 = tmp0[ 1 ] * tmp30;
        assert(tmp50 == tmp50);
      const auto tmp51 = 3.141592653589793 * tmp0[ 1 ];
        assert(tmp51 == tmp51);
      const auto tmp52 = tmp51 * tmp36;
        assert(tmp52 == tmp52);
      const auto tmp53 = tmp31 * tmp52;
        assert(tmp53 == tmp53);
      const auto tmp54 = tmp53 + tmp50;
        assert(tmp54 == tmp54);
      const auto tmp55 = tmp54 * tmp20;
        assert(tmp55 == tmp55);
      const auto tmp56 = tmp0[ 1 ] * tmp52;
        assert(tmp56 == tmp56);
      const auto tmp57 = tmp51 * tmp30;
        assert(tmp57 == tmp57);
      const auto tmp58 = -1 * tmp57;
        assert(tmp58 == tmp58);
      const auto tmp59 = tmp51 * tmp58;
        assert(tmp59 == tmp59);
      const auto tmp60 = tmp31 * tmp59;
        assert(tmp60 == tmp60);
      const auto tmp61 = tmp60 + tmp56;
        assert(tmp61 == tmp61);
      const auto tmp62 = tmp61 + tmp56;
        assert(tmp62 == tmp62);
      const auto tmp63 = tmp62 * tmp14;
        assert(tmp63 == tmp63);
      const auto tmp64 = tmp63 + tmp55;
        assert(tmp64 == tmp64);
      const auto tmp65 = tmp64 + tmp49;
        assert(tmp65 == tmp65);
      const auto tmp66 = -1 * tmp65;
        assert(tmp66 == tmp66);
      const auto tmp67 = tmp66 + tmp33;
        assert(tmp67 == tmp67);
      const auto tmp68 = -1 * tmp67;
        assert(tmp68 == tmp68);
      const auto tmp69 = tmp68 + tmp27;
        assert(tmp69 == tmp69);
      const auto tmp70 = tmp69 + tmp25;
        assert(tmp70 == tmp70);
      result[ 0 ] = tmp70;
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

    template< class Point >
    void alpha ( const Point &x, const DRangeType &u, RRangeType &result ) const
    {
      result[ 0 ] = 0.0;
    }

    template< class Point >
    void linAlpha ( const DRangeType &ubar, const Point &x, const DRangeType &u, RRangeType &result ) const
    {
      result[ 0 ] = 0.0;
    }

    bool hasNeumanBoundary () const
    {
      return false;
    }
    typedef std::array<int,1> DirichletComponentType;

    bool hasDirichletBoundary () const
    {
      return true;
    }

    bool isDirichletIntersection ( const IntersectionType &intersection, DirichletComponentType &dirichletComponent ) const
    {
      const int bndId = BoundaryIdProviderType::boundaryId( intersection );
      {
        std::fill( dirichletComponent.begin(), dirichletComponent.end(), bndId );
      }
      switch( bndId )
      {
      case 1:
        {
          return true;
        }
        break;
      case 2:
        {
          return true;
        }
        break;
      case 3:
        {
          return true;
        }
        break;
      case 4:
        {
          return true;
        }
        break;
      default:
        {
          return false;
        }
      }
    }

    template< class Point >
    void dirichlet ( int bndId, const Point &x, RRangeType &result ) const
    {
      switch( bndId )
      {
      case 1:
        {
          using std::cos;
          auto tmp0 = entity().geometry().global( Dune::Fem::coordinate( x ) );
            assert(tmp0 == tmp0);
          const auto tmp1 = 3.141592653589793 * tmp0[ 0 ];
            assert(tmp1 == tmp1);
          const auto tmp2 = tmp0[ 1 ] * tmp1;
            assert(tmp2 == tmp2);
          const auto tmp3 = std::cos( tmp2 );
            assert(tmp3 == tmp3);
          const auto tmp4 = tmp0[ 0 ] * tmp0[ 1 ];
            assert(tmp4 == tmp4);
          const auto tmp5 = tmp4 * tmp3;
            assert(tmp5 == tmp5);
          result[ 0 ] = tmp5;
        }
        break;
      case 2:
        {
          using std::cos;
          auto tmp0 = entity().geometry().global( Dune::Fem::coordinate( x ) );
            assert(tmp0 == tmp0);
          const auto tmp1 = 3.141592653589793 * tmp0[ 0 ];
            assert(tmp1 == tmp1);
          const auto tmp2 = tmp0[ 1 ] * tmp1;
            assert(tmp2 == tmp2);
          const auto tmp3 = std::cos( tmp2 );
            assert(tmp3 == tmp3);
          const auto tmp4 = tmp0[ 0 ] * tmp0[ 1 ];
            assert(tmp4 == tmp4);
          const auto tmp5 = tmp4 * tmp3;
            assert(tmp5 == tmp5);
          result[ 0 ] = tmp5;
        }
        break;
      case 3:
        {
          using std::cos;
          auto tmp0 = entity().geometry().global( Dune::Fem::coordinate( x ) );
            assert(tmp0 == tmp0);
          const auto tmp1 = 3.141592653589793 * tmp0[ 0 ];
            assert(tmp1 == tmp1);
          const auto tmp2 = tmp0[ 1 ] * tmp1;
            assert(tmp2 == tmp2);
          const auto tmp3 = std::cos( tmp2 );
            assert(tmp3 == tmp3);
          const auto tmp4 = tmp0[ 0 ] * tmp0[ 1 ];
            assert(tmp4 == tmp4);
          const auto tmp5 = tmp4 * tmp3;
            assert(tmp5 == tmp5);
          result[ 0 ] = tmp5;
        }
        break;
      case 4:
        {
          using std::cos;
          auto tmp0 = entity().geometry().global( Dune::Fem::coordinate( x ) );
            assert(tmp0 == tmp0);
          const auto tmp1 = 3.141592653589793 * tmp0[ 0 ];
            assert(tmp1 == tmp1);
          const auto tmp2 = tmp0[ 1 ] * tmp1;
            assert(tmp2 == tmp2);
          const auto tmp3 = std::cos( tmp2 );
            assert(tmp3 == tmp3);
          const auto tmp4 = tmp0[ 0 ] * tmp0[ 1 ];
            assert(tmp4 == tmp4);
          const auto tmp5 = tmp4 * tmp3;
            assert(tmp5 == tmp5);
          result[ 0 ] = tmp5;
        }
        break;
      default:
        {
          result = RRangeType( 0 );
        }
      }
    }

  private:
    mutable const EntityType *entity_ = nullptr;

  public:
    template< class Point, class DRangeType >
    RRangeType penalty ( const Point &x, const DRangeType &u ) const
    {
      RRangeType result;
      using std::cos;
      auto tmp0 = entity().geometry().global( Dune::Fem::coordinate( x ) );
        assert(tmp0 == tmp0);
      const auto tmp1 = tmp0[ 1 ] * tmp0[ 1 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = tmp0[ 0 ] * tmp0[ 0 ];
        assert(tmp2 == tmp2);
      const auto tmp3 = tmp2 + tmp1;
        assert(tmp3 == tmp3);
      const auto tmp4 = std::cos( tmp3 );
        assert(tmp4 == tmp4);
      const auto tmp5 = 0.9 * tmp4;
        assert(tmp5 == tmp5);
      const auto tmp6 = -1 * tmp5;
        assert(tmp6 == tmp6);
      const auto tmp7 = 1 + tmp6;
        assert(tmp7 == tmp7);
      result[ 0 ] = tmp7;
      return result;
    }

    template< class Point, class DRangeType >
    RRangeType linPenalty ( const Point &x, const DRangeType &u ) const
    {
      RRangeType result;
      result[ 0 ] = 0.0;
      return result;
    }
  };

} // namespace ModelImpl_719110241111e05224e2f1c3078e038fv1_3_f084d41ce9ff164524e5c3ea6f8082cc

namespace ModelImpl_719110241111e05224e2f1c3078e038fv1_3_f084d41ce9ff164524e5c3ea6f8082cc
{

} // namespace ModelImpl_719110241111e05224e2f1c3078e038fv1_3_f084d41ce9ff164524e5c3ea6f8082cc

PYBIND11_MODULE( modelDG_719110241111e05224e2f1c3078e038fv1_3_f084d41ce9ff164524e5c3ea6f8082cc, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView > GridPart;
  typedef ModelImpl_719110241111e05224e2f1c3078e038fv1_3_f084d41ce9ff164524e5c3ea6f8082cc::Model< GridPart > Model;
  typedef DGConservationLawModelWrapper< Model > ModelWrapper;
  typedef typename ModelWrapper::Base ModelBase;
  // export abstract base class
  if( !pybind11::already_registered< ModelBase >() )
    pybind11::class_< ModelBase >( module, "ModelBase" );

  // actual wrapper class for model derived from abstract base
  auto cls = Dune::Python::insertClass<ModelWrapper,ModelBase>(module,"Model",Dune::Python::GenerateTypeName("DGConservationLawModelWrapper< Model >"),Dune::Python::IncludeFiles({})).first;
  cls.def_property_readonly( "dimRange", [] ( ModelWrapper & ) { return 1; } );
  cls.def_property_readonly( "hasDirichletBoundary", [] ( ModelWrapper& ) -> bool { return true;});



  cls.def( pybind11::init( [] (  ) {
    return new  ModelWrapper();
    } ) );
}
#endif // GuardModelImpl_719110241111e05224e2f1c3078e038fv1_3_f084d41ce9ff164524e5c3ea6f8082cc
