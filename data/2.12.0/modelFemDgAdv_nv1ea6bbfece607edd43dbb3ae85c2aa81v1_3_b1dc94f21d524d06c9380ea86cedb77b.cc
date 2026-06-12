#ifndef GuardModelImpl_nv1ea6bbfece607edd43dbb3ae85c2aa81v1_3_b1dc94f21d524d06c9380ea86cedb77b
#define GuardModelImpl_nv1ea6bbfece607edd43dbb3ae85c2aa81v1_3_b1dc94f21d524d06c9380ea86cedb77b
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

namespace ModelImpl_nv1ea6bbfece607edd43dbb3ae85c2aa81v1_3_b1dc94f21d524d06c9380ea86cedb77b
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
      const auto tmp0 = u[ 1 ] / u[ 0 ];
        assert(tmp0 == tmp0);
      const auto tmp1 = u[ 0 ] * tmp0;
        assert(tmp1 == tmp1);
      const auto tmp2 = u[ 2 ] / u[ 0 ];
        assert(tmp2 == tmp2);
      const auto tmp3 = u[ 0 ] * tmp2;
        assert(tmp3 == tmp3);
      const auto tmp4 = tmp1 * tmp0;
        assert(tmp4 == tmp4);
      const auto tmp5 = tmp2 * tmp2;
        assert(tmp5 == tmp5);
      const auto tmp6 = tmp0 * tmp0;
        assert(tmp6 == tmp6);
      const auto tmp7 = tmp6 + tmp5;
        assert(tmp7 == tmp7);
      const auto tmp8 = u[ 0 ] * tmp7;
        assert(tmp8 == tmp8);
      const auto tmp9 = tmp8 / 2;
        assert(tmp9 == tmp9);
      const auto tmp10 = -1 * tmp9;
        assert(tmp10 == tmp10);
      const auto tmp11 = u[ 3 ] + tmp10;
        assert(tmp11 == tmp11);
      const auto tmp12 = 0.3999999999999999 * tmp11;
        assert(tmp12 == tmp12);
      const auto tmp13 = tmp12 + tmp4;
        assert(tmp13 == tmp13);
      const auto tmp14 = tmp1 * tmp2;
        assert(tmp14 == tmp14);
      const auto tmp15 = tmp3 * tmp2;
        assert(tmp15 == tmp15);
      const auto tmp16 = tmp12 + tmp15;
        assert(tmp16 == tmp16);
      const auto tmp17 = u[ 3 ] + tmp12;
        assert(tmp17 == tmp17);
      const auto tmp18 = tmp17 * tmp0;
        assert(tmp18 == tmp18);
      const auto tmp19 = tmp17 * tmp2;
        assert(tmp19 == tmp19);
      (result[ 0 ])[ 0 ] = tmp1;
      (result[ 0 ])[ 1 ] = tmp3;
      (result[ 1 ])[ 0 ] = tmp13;
      (result[ 1 ])[ 1 ] = tmp14;
      (result[ 2 ])[ 0 ] = tmp14;
      (result[ 2 ])[ 1 ] = tmp16;
      (result[ 3 ])[ 0 ] = tmp18;
      (result[ 3 ])[ 1 ] = tmp19;
    }

    template< class Point >
    void linFlux ( const DRangeType &ubar, const DJacobianRangeType &dubar, const Point &x, const DRangeType &u, const DJacobianRangeType &du, RJacobianRangeType &result ) const
    {
      const auto tmp0 = ubar[ 1 ] / ubar[ 0 ];
        assert(tmp0 == tmp0);
      const auto tmp1 = u[ 0 ] * tmp0;
        assert(tmp1 == tmp1);
      const auto tmp2 = -1 * tmp1;
        assert(tmp2 == tmp2);
      const auto tmp3 = u[ 1 ] + tmp2;
        assert(tmp3 == tmp3);
      const auto tmp4 = tmp3 / ubar[ 0 ];
        assert(tmp4 == tmp4);
      const auto tmp5 = ubar[ 0 ] * tmp4;
        assert(tmp5 == tmp5);
      const auto tmp6 = tmp1 + tmp5;
        assert(tmp6 == tmp6);
      const auto tmp7 = ubar[ 2 ] / ubar[ 0 ];
        assert(tmp7 == tmp7);
      const auto tmp8 = u[ 0 ] * tmp7;
        assert(tmp8 == tmp8);
      const auto tmp9 = -1 * tmp8;
        assert(tmp9 == tmp9);
      const auto tmp10 = u[ 2 ] + tmp9;
        assert(tmp10 == tmp10);
      const auto tmp11 = tmp10 / ubar[ 0 ];
        assert(tmp11 == tmp11);
      const auto tmp12 = ubar[ 0 ] * tmp11;
        assert(tmp12 == tmp12);
      const auto tmp13 = tmp8 + tmp12;
        assert(tmp13 == tmp13);
      const auto tmp14 = tmp7 * tmp7;
        assert(tmp14 == tmp14);
      const auto tmp15 = tmp0 * tmp0;
        assert(tmp15 == tmp15);
      const auto tmp16 = tmp15 + tmp14;
        assert(tmp16 == tmp16);
      const auto tmp17 = u[ 0 ] * tmp16;
        assert(tmp17 == tmp17);
      const auto tmp18 = tmp7 * tmp11;
        assert(tmp18 == tmp18);
      const auto tmp19 = tmp18 + tmp18;
        assert(tmp19 == tmp19);
      const auto tmp20 = tmp0 * tmp4;
        assert(tmp20 == tmp20);
      const auto tmp21 = tmp20 + tmp20;
        assert(tmp21 == tmp21);
      const auto tmp22 = tmp21 + tmp19;
        assert(tmp22 == tmp22);
      const auto tmp23 = ubar[ 0 ] * tmp22;
        assert(tmp23 == tmp23);
      const auto tmp24 = tmp23 + tmp17;
        assert(tmp24 == tmp24);
      const auto tmp25 = tmp24 / 2;
        assert(tmp25 == tmp25);
      const auto tmp26 = -1 * tmp25;
        assert(tmp26 == tmp26);
      const auto tmp27 = u[ 3 ] + tmp26;
        assert(tmp27 == tmp27);
      const auto tmp28 = 0.3999999999999999 * tmp27;
        assert(tmp28 == tmp28);
      const auto tmp29 = ubar[ 0 ] * tmp0;
        assert(tmp29 == tmp29);
      const auto tmp30 = tmp29 * tmp4;
        assert(tmp30 == tmp30);
      const auto tmp31 = tmp6 * tmp0;
        assert(tmp31 == tmp31);
      const auto tmp32 = tmp31 + tmp30;
        assert(tmp32 == tmp32);
      const auto tmp33 = tmp32 + tmp28;
        assert(tmp33 == tmp33);
      const auto tmp34 = tmp29 * tmp11;
        assert(tmp34 == tmp34);
      const auto tmp35 = tmp6 * tmp7;
        assert(tmp35 == tmp35);
      const auto tmp36 = tmp35 + tmp34;
        assert(tmp36 == tmp36);
      const auto tmp37 = ubar[ 0 ] * tmp7;
        assert(tmp37 == tmp37);
      const auto tmp38 = tmp37 * tmp11;
        assert(tmp38 == tmp38);
      const auto tmp39 = tmp13 * tmp7;
        assert(tmp39 == tmp39);
      const auto tmp40 = tmp39 + tmp38;
        assert(tmp40 == tmp40);
      const auto tmp41 = tmp40 + tmp28;
        assert(tmp41 == tmp41);
      const auto tmp42 = ubar[ 0 ] * tmp16;
        assert(tmp42 == tmp42);
      const auto tmp43 = tmp42 / 2;
        assert(tmp43 == tmp43);
      const auto tmp44 = -1 * tmp43;
        assert(tmp44 == tmp44);
      const auto tmp45 = ubar[ 3 ] + tmp44;
        assert(tmp45 == tmp45);
      const auto tmp46 = 0.3999999999999999 * tmp45;
        assert(tmp46 == tmp46);
      const auto tmp47 = ubar[ 3 ] + tmp46;
        assert(tmp47 == tmp47);
      const auto tmp48 = tmp47 * tmp4;
        assert(tmp48 == tmp48);
      const auto tmp49 = u[ 3 ] + tmp28;
        assert(tmp49 == tmp49);
      const auto tmp50 = tmp49 * tmp0;
        assert(tmp50 == tmp50);
      const auto tmp51 = tmp50 + tmp48;
        assert(tmp51 == tmp51);
      const auto tmp52 = tmp47 * tmp11;
        assert(tmp52 == tmp52);
      const auto tmp53 = tmp49 * tmp7;
        assert(tmp53 == tmp53);
      const auto tmp54 = tmp53 + tmp52;
        assert(tmp54 == tmp54);
      (result[ 0 ])[ 0 ] = tmp6;
      (result[ 0 ])[ 1 ] = tmp13;
      (result[ 1 ])[ 0 ] = tmp33;
      (result[ 1 ])[ 1 ] = tmp36;
      (result[ 2 ])[ 0 ] = tmp36;
      (result[ 2 ])[ 1 ] = tmp41;
      (result[ 3 ])[ 0 ] = tmp51;
      (result[ 3 ])[ 1 ] = tmp54;
    }

    template< class Point >
    void fluxDivergence ( const Point &x, const DRangeType &u, const DJacobianRangeType &du, const DHessianRangeType &d2u, RRangeType &result ) const
    {
      const auto tmp0 = u[ 2 ] / u[ 0 ];
        assert(tmp0 == tmp0);
      const auto tmp1 = (du[ 0 ])[ 1 ] * tmp0;
        assert(tmp1 == tmp1);
      const auto tmp2 = -1 * tmp1;
        assert(tmp2 == tmp2);
      const auto tmp3 = (du[ 2 ])[ 1 ] + tmp2;
        assert(tmp3 == tmp3);
      const auto tmp4 = tmp3 / u[ 0 ];
        assert(tmp4 == tmp4);
      const auto tmp5 = u[ 0 ] * tmp4;
        assert(tmp5 == tmp5);
      const auto tmp6 = tmp1 + tmp5;
        assert(tmp6 == tmp6);
      const auto tmp7 = u[ 1 ] / u[ 0 ];
        assert(tmp7 == tmp7);
      const auto tmp8 = (du[ 0 ])[ 0 ] * tmp7;
        assert(tmp8 == tmp8);
      const auto tmp9 = -1 * tmp8;
        assert(tmp9 == tmp9);
      const auto tmp10 = (du[ 1 ])[ 0 ] + tmp9;
        assert(tmp10 == tmp10);
      const auto tmp11 = tmp10 / u[ 0 ];
        assert(tmp11 == tmp11);
      const auto tmp12 = u[ 0 ] * tmp11;
        assert(tmp12 == tmp12);
      const auto tmp13 = tmp8 + tmp12;
        assert(tmp13 == tmp13);
      const auto tmp14 = tmp13 + tmp6;
        assert(tmp14 == tmp14);
      const auto tmp15 = -1 * tmp14;
        assert(tmp15 == tmp15);
      const auto tmp16 = u[ 0 ] * tmp7;
        assert(tmp16 == tmp16);
      const auto tmp17 = tmp16 * tmp4;
        assert(tmp17 == tmp17);
      const auto tmp18 = (du[ 0 ])[ 1 ] * tmp7;
        assert(tmp18 == tmp18);
      const auto tmp19 = -1 * tmp18;
        assert(tmp19 == tmp19);
      const auto tmp20 = (du[ 1 ])[ 1 ] + tmp19;
        assert(tmp20 == tmp20);
      const auto tmp21 = tmp20 / u[ 0 ];
        assert(tmp21 == tmp21);
      const auto tmp22 = u[ 0 ] * tmp21;
        assert(tmp22 == tmp22);
      const auto tmp23 = tmp18 + tmp22;
        assert(tmp23 == tmp23);
      const auto tmp24 = tmp23 * tmp0;
        assert(tmp24 == tmp24);
      const auto tmp25 = tmp24 + tmp17;
        assert(tmp25 == tmp25);
      const auto tmp26 = tmp0 * tmp0;
        assert(tmp26 == tmp26);
      const auto tmp27 = tmp7 * tmp7;
        assert(tmp27 == tmp27);
      const auto tmp28 = tmp27 + tmp26;
        assert(tmp28 == tmp28);
      const auto tmp29 = (du[ 0 ])[ 0 ] * tmp28;
        assert(tmp29 == tmp29);
      const auto tmp30 = (du[ 0 ])[ 0 ] * tmp0;
        assert(tmp30 == tmp30);
      const auto tmp31 = -1 * tmp30;
        assert(tmp31 == tmp31);
      const auto tmp32 = (du[ 2 ])[ 0 ] + tmp31;
        assert(tmp32 == tmp32);
      const auto tmp33 = tmp32 / u[ 0 ];
        assert(tmp33 == tmp33);
      const auto tmp34 = tmp0 * tmp33;
        assert(tmp34 == tmp34);
      const auto tmp35 = tmp34 + tmp34;
        assert(tmp35 == tmp35);
      const auto tmp36 = tmp7 * tmp11;
        assert(tmp36 == tmp36);
      const auto tmp37 = tmp36 + tmp36;
        assert(tmp37 == tmp37);
      const auto tmp38 = tmp37 + tmp35;
        assert(tmp38 == tmp38);
      const auto tmp39 = u[ 0 ] * tmp38;
        assert(tmp39 == tmp39);
      const auto tmp40 = tmp39 + tmp29;
        assert(tmp40 == tmp40);
      const auto tmp41 = tmp40 / 2;
        assert(tmp41 == tmp41);
      const auto tmp42 = -1 * tmp41;
        assert(tmp42 == tmp42);
      const auto tmp43 = (du[ 3 ])[ 0 ] + tmp42;
        assert(tmp43 == tmp43);
      const auto tmp44 = 0.3999999999999999 * tmp43;
        assert(tmp44 == tmp44);
      const auto tmp45 = tmp16 * tmp11;
        assert(tmp45 == tmp45);
      const auto tmp46 = tmp13 * tmp7;
        assert(tmp46 == tmp46);
      const auto tmp47 = tmp46 + tmp45;
        assert(tmp47 == tmp47);
      const auto tmp48 = tmp47 + tmp44;
        assert(tmp48 == tmp48);
      const auto tmp49 = tmp48 + tmp25;
        assert(tmp49 == tmp49);
      const auto tmp50 = -1 * tmp49;
        assert(tmp50 == tmp50);
      const auto tmp51 = tmp16 * tmp33;
        assert(tmp51 == tmp51);
      const auto tmp52 = tmp13 * tmp0;
        assert(tmp52 == tmp52);
      const auto tmp53 = tmp52 + tmp51;
        assert(tmp53 == tmp53);
      const auto tmp54 = (du[ 0 ])[ 1 ] * tmp28;
        assert(tmp54 == tmp54);
      const auto tmp55 = tmp0 * tmp4;
        assert(tmp55 == tmp55);
      const auto tmp56 = tmp55 + tmp55;
        assert(tmp56 == tmp56);
      const auto tmp57 = tmp7 * tmp21;
        assert(tmp57 == tmp57);
      const auto tmp58 = tmp57 + tmp57;
        assert(tmp58 == tmp58);
      const auto tmp59 = tmp58 + tmp56;
        assert(tmp59 == tmp59);
      const auto tmp60 = u[ 0 ] * tmp59;
        assert(tmp60 == tmp60);
      const auto tmp61 = tmp60 + tmp54;
        assert(tmp61 == tmp61);
      const auto tmp62 = tmp61 / 2;
        assert(tmp62 == tmp62);
      const auto tmp63 = -1 * tmp62;
        assert(tmp63 == tmp63);
      const auto tmp64 = (du[ 3 ])[ 1 ] + tmp63;
        assert(tmp64 == tmp64);
      const auto tmp65 = 0.3999999999999999 * tmp64;
        assert(tmp65 == tmp65);
      const auto tmp66 = u[ 0 ] * tmp0;
        assert(tmp66 == tmp66);
      const auto tmp67 = tmp66 * tmp4;
        assert(tmp67 == tmp67);
      const auto tmp68 = tmp6 * tmp0;
        assert(tmp68 == tmp68);
      const auto tmp69 = tmp68 + tmp67;
        assert(tmp69 == tmp69);
      const auto tmp70 = tmp69 + tmp65;
        assert(tmp70 == tmp70);
      const auto tmp71 = tmp70 + tmp53;
        assert(tmp71 == tmp71);
      const auto tmp72 = -1 * tmp71;
        assert(tmp72 == tmp72);
      const auto tmp73 = u[ 0 ] * tmp28;
        assert(tmp73 == tmp73);
      const auto tmp74 = tmp73 / 2;
        assert(tmp74 == tmp74);
      const auto tmp75 = -1 * tmp74;
        assert(tmp75 == tmp75);
      const auto tmp76 = u[ 3 ] + tmp75;
        assert(tmp76 == tmp76);
      const auto tmp77 = 0.3999999999999999 * tmp76;
        assert(tmp77 == tmp77);
      const auto tmp78 = u[ 3 ] + tmp77;
        assert(tmp78 == tmp78);
      const auto tmp79 = tmp78 * tmp4;
        assert(tmp79 == tmp79);
      const auto tmp80 = (du[ 3 ])[ 1 ] + tmp65;
        assert(tmp80 == tmp80);
      const auto tmp81 = tmp80 * tmp0;
        assert(tmp81 == tmp81);
      const auto tmp82 = tmp81 + tmp79;
        assert(tmp82 == tmp82);
      const auto tmp83 = tmp78 * tmp11;
        assert(tmp83 == tmp83);
      const auto tmp84 = (du[ 3 ])[ 0 ] + tmp44;
        assert(tmp84 == tmp84);
      const auto tmp85 = tmp84 * tmp7;
        assert(tmp85 == tmp85);
      const auto tmp86 = tmp85 + tmp83;
        assert(tmp86 == tmp86);
      const auto tmp87 = tmp86 + tmp82;
        assert(tmp87 == tmp87);
      const auto tmp88 = -1 * tmp87;
        assert(tmp88 == tmp88);
      result[ 0 ] = tmp15;
      result[ 1 ] = tmp50;
      result[ 2 ] = tmp72;
      result[ 3 ] = tmp88;
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

  private:
    mutable const EntityType *entity_ = nullptr;

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

} // namespace ModelImpl_nv1ea6bbfece607edd43dbb3ae85c2aa81v1_3_b1dc94f21d524d06c9380ea86cedb77b

namespace ModelImpl_nv1ea6bbfece607edd43dbb3ae85c2aa81v1_3_b1dc94f21d524d06c9380ea86cedb77b
{

} // namespace ModelImpl_nv1ea6bbfece607edd43dbb3ae85c2aa81v1_3_b1dc94f21d524d06c9380ea86cedb77b

PYBIND11_MODULE( modelFemDgAdv_nv1ea6bbfece607edd43dbb3ae85c2aa81v1_3_b1dc94f21d524d06c9380ea86cedb77b, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView > GridPart;
  typedef ModelImpl_nv1ea6bbfece607edd43dbb3ae85c2aa81v1_3_b1dc94f21d524d06c9380ea86cedb77b::Model< GridPart > Model;
  typedef Model ModelWrapper;
  auto cls = Dune::Python::insertClass<ModelWrapper>(module,"Model",Dune::Python::GenerateTypeName("ModelImpl_nv1ea6bbfece607edd43dbb3ae85c2aa81v1_3_b1dc94f21d524d06c9380ea86cedb77b::Model< typename Dune::FemPy::GridPart< Dune::ALUGrid< 2, 2, Dune::simplex >::LeafGridView > >"),Dune::Python::IncludeFiles({"python/dune/generated/modelFemDgAdv_nv1ea6bbfece607edd43dbb3ae85c2aa81v1_3_b1dc94f21d524d06c9380ea86cedb77b.cc"})).first;
  cls.def_property_readonly( "dimRange", [] ( ModelWrapper & ) { return 4; } );
  cls.def_property_readonly( "hasDirichletBoundary", [] ( ModelWrapper& ) -> bool { return false;});



  cls.def( pybind11::init( [] (  ) {
    return new  ModelWrapper();
    } ) );
}
#endif // GuardModelImpl_nv1ea6bbfece607edd43dbb3ae85c2aa81v1_3_b1dc94f21d524d06c9380ea86cedb77b
