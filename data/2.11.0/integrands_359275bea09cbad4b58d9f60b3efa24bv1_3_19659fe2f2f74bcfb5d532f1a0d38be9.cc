#ifndef GuardIntegrands_359275bea09cbad4b58d9f60b3efa24bv1_3
#define GuardIntegrands_359275bea09cbad4b58d9f60b3efa24bv1_3
#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/python/grid/hierarchical.hh>
#include <cmath>
#include <tuple>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/extensions.h>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/geometry/edgelength.hh>
#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/function/simplegridfunction.hh>
#include <dune/fem/misc/gridfunctionview.hh>
#include <dune/fempy/py/integrands.hh>

namespace Integrands_359275bea09cbad4b58d9f60b3efa24bv1_3
{




  // Integrands
  // ----------

  template< class GridPart, class Coeffprevious >
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
    typedef double Condt;
    typedef double Cont;
    typedef std::tuple< std::shared_ptr< Condt >, std::shared_ptr< Cont > > ConstantTupleType;
    template< std::size_t i >
    using ConstantsRangeType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;
    typedef std::tuple< Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > > > CoefficientFunctionSpaceTupleType;
    typedef std::tuple< Coeffprevious > CoefficientTupleType;
    template< std::size_t i >
    using CoefficientFunctionSpaceType = std::tuple_element_t< i, CoefficientFunctionSpaceTupleType >;
    template< std::size_t i >
    using CoefficientRangeType = typename CoefficientFunctionSpaceType< i >::RangeType;
    template< std::size_t i >
    using CoefficientJacobianRangeType = typename CoefficientFunctionSpaceType< i >::JacobianRangeType;
    static constexpr bool gridPartValid = Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffprevious>>();
    template< std::size_t i >
    using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
    template< std::size_t i >
    using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;

    Integrands ( const Coeffprevious &coeffprevious, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
      : coefficients_( Dune::Fem::ConstLocalFunction< Coeffprevious >( coeffprevious ) )
    {
      std::get< 0 >( constants_ ) = std::make_shared< Condt >( (Condt(0)) );
      std::get< 1 >( constants_ ) = std::make_shared< Cont >( (Cont(0)) );
    }

    bool init ( const EntityType &entity )
    {
      entity_ = entity;
      geometry_.emplace( this->entity().geometry() );
      std::get< 0 >( coefficients_ ).bind( this->entity() );
      return true;
    }

    void unbind ()
    {
      std::get< 0 >( coefficients_ ).unbind();
    }

    bool init ( const IntersectionType &intersection )
    {
      intersection_ = intersection;
      return (intersection.boundary() && init( intersection.inside() ));
    }
    typedef std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 > > DomainValueType;
    typedef std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 > > RangeValueType;
    static constexpr bool _nonlinear = true;

    bool nonlinear () const
    {
      return _nonlinear;
    }

    template< class Point >
    RangeValueType interior ( const Point &x, const DomainValueType &u ) const
    {
      using std::exp;
      using std::pow;
      using std::sqrt;
      double tmp0 = constant< 0 >();
      typename CoefficientFunctionSpaceType< 0 >::RangeType tmp1 = evaluateCoefficient< 0 >( x );
      const auto tmp2 = -1 * tmp1[ 0 ];
      const auto tmp3 = std::get< 0 >( u );
      const auto tmp4 = tmp3[ 0 ] + tmp2;
      const auto tmp5 = tmp4 / tmp0;
      double tmp6 = constant< 1 >();
      const auto tmp7 = tmp0 + tmp6;
      const auto tmp8 = -2 * tmp7;
      const auto tmp9 = std::exp( tmp8 );
      const auto tmp10 = -2 * tmp9;
      GlobalCoordinateType tmp11 = geometry().global( Dune::Fem::coordinate( x ) );
      const auto tmp12 = std::pow( tmp11[ 0 ], 3 );
      const auto tmp13 = std::pow( tmp11[ 1 ], 3 );
      const auto tmp14 = -1 * tmp13;
      const auto tmp15 = tmp14 + tmp12;
      const auto tmp16 = 0.3333333333333333 * tmp15;
      const auto tmp17 = -1 * tmp16;
      const auto tmp18 = std::pow( tmp11[ 1 ], 2 );
      const auto tmp19 = std::pow( tmp11[ 0 ], 2 );
      const auto tmp20 = tmp19 + tmp18;
      const auto tmp21 = 0.5 * tmp20;
      const auto tmp22 = tmp21 + tmp17;
      const auto tmp23 = 1 + tmp22;
      const auto tmp24 = -1 + tmp23;
      const auto tmp25 = tmp24 * tmp10;
      const auto tmp26 = 3 * tmp19;
      const auto tmp27 = 0.3333333333333333 * tmp26;
      const auto tmp28 = -1 * tmp27;
      const auto tmp29 = 2 * tmp11[ 0 ];
      const auto tmp30 = 0.5 * tmp29;
      const auto tmp31 = tmp30 + tmp28;
      const auto tmp32 = tmp31 * tmp9;
      const auto tmp33 = tmp32 * tmp32;
      const auto tmp34 = 3 * tmp18;
      const auto tmp35 = -1 * tmp34;
      const auto tmp36 = 0.3333333333333333 * tmp35;
      const auto tmp37 = -1 * tmp36;
      const auto tmp38 = 2 * tmp11[ 1 ];
      const auto tmp39 = 0.5 * tmp38;
      const auto tmp40 = tmp39 + tmp37;
      const auto tmp41 = tmp40 * tmp9;
      const auto tmp42 = tmp41 * tmp41;
      const auto tmp43 = tmp42 + tmp33;
      const auto tmp44 = std::sqrt( tmp43 );
      const auto tmp45 = 4 * tmp44;
      const auto tmp46 = 1 + tmp45;
      const auto tmp47 = std::sqrt( tmp46 );
      const auto tmp48 = 1 + tmp47;
      const auto tmp49 = 2 * tmp47;
      const auto tmp50 = 2 * tmp44;
      const auto tmp51 = 3 * tmp29;
      const auto tmp52 = 0.3333333333333333 * tmp51;
      const auto tmp53 = -1 * tmp52;
      const auto tmp54 = 1.0 + tmp53;
      const auto tmp55 = tmp54 * tmp9;
      const auto tmp56 = tmp55 * tmp32;
      const auto tmp57 = tmp56 + tmp56;
      const auto tmp58 = tmp57 / tmp50;
      const auto tmp59 = 4 * tmp58;
      const auto tmp60 = tmp59 / tmp49;
      const auto tmp61 = 2 / tmp48;
      const auto tmp62 = tmp61 * tmp60;
      const auto tmp63 = -1 * tmp62;
      const auto tmp64 = tmp63 / tmp48;
      const auto tmp65 = tmp32 * tmp64;
      const auto tmp66 = tmp55 * tmp61;
      const auto tmp67 = tmp66 + tmp65;
      const auto tmp68 = 3 * tmp38;
      const auto tmp69 = -1 * tmp68;
      const auto tmp70 = 0.3333333333333333 * tmp69;
      const auto tmp71 = -1 * tmp70;
      const auto tmp72 = 1.0 + tmp71;
      const auto tmp73 = tmp72 * tmp9;
      const auto tmp74 = tmp73 * tmp41;
      const auto tmp75 = tmp74 + tmp74;
      const auto tmp76 = tmp75 / tmp50;
      const auto tmp77 = 4 * tmp76;
      const auto tmp78 = tmp77 / tmp49;
      const auto tmp79 = tmp61 * tmp78;
      const auto tmp80 = -1 * tmp79;
      const auto tmp81 = tmp80 / tmp48;
      const auto tmp82 = tmp41 * tmp81;
      const auto tmp83 = tmp73 * tmp61;
      const auto tmp84 = tmp83 + tmp82;
      const auto tmp85 = tmp84 + tmp67;
      const auto tmp86 = -1 * tmp85;
      const auto tmp87 = tmp86 + tmp25;
      const auto tmp88 = -2 * tmp6;
      const auto tmp89 = std::exp( tmp88 );
      const auto tmp90 = -2 * tmp89;
      const auto tmp91 = tmp24 * tmp90;
      const auto tmp92 = tmp31 * tmp89;
      const auto tmp93 = tmp92 * tmp92;
      const auto tmp94 = tmp40 * tmp89;
      const auto tmp95 = tmp94 * tmp94;
      const auto tmp96 = tmp95 + tmp93;
      const auto tmp97 = std::sqrt( tmp96 );
      const auto tmp98 = 4 * tmp97;
      const auto tmp99 = 1 + tmp98;
      const auto tmp100 = std::sqrt( tmp99 );
      const auto tmp101 = 1 + tmp100;
      const auto tmp102 = 2 * tmp100;
      const auto tmp103 = 2 * tmp97;
      const auto tmp104 = tmp54 * tmp89;
      const auto tmp105 = tmp104 * tmp92;
      const auto tmp106 = tmp105 + tmp105;
      const auto tmp107 = tmp106 / tmp103;
      const auto tmp108 = 4 * tmp107;
      const auto tmp109 = tmp108 / tmp102;
      const auto tmp110 = 2 / tmp101;
      const auto tmp111 = tmp110 * tmp109;
      const auto tmp112 = -1 * tmp111;
      const auto tmp113 = tmp112 / tmp101;
      const auto tmp114 = tmp92 * tmp113;
      const auto tmp115 = tmp104 * tmp110;
      const auto tmp116 = tmp115 + tmp114;
      const auto tmp117 = tmp72 * tmp89;
      const auto tmp118 = tmp117 * tmp94;
      const auto tmp119 = tmp118 + tmp118;
      const auto tmp120 = tmp119 / tmp103;
      const auto tmp121 = 4 * tmp120;
      const auto tmp122 = tmp121 / tmp102;
      const auto tmp123 = tmp110 * tmp122;
      const auto tmp124 = -1 * tmp123;
      const auto tmp125 = tmp124 / tmp101;
      const auto tmp126 = tmp94 * tmp125;
      const auto tmp127 = tmp117 * tmp110;
      const auto tmp128 = tmp127 + tmp126;
      const auto tmp129 = tmp128 + tmp116;
      const auto tmp130 = -1 * tmp129;
      const auto tmp131 = tmp130 + tmp91;
      const auto tmp132 = tmp131 + tmp87;
      const auto tmp133 = 0.5 * tmp132;
      const auto tmp134 = -1 * tmp133;
      const auto tmp135 = tmp134 + tmp5;
      typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp136 = jacobianCoefficient< 0 >( x );
      const auto tmp137 = (tmp136[ 0 ])[ 1 ] * (tmp136[ 0 ])[ 1 ];
      const auto tmp138 = (tmp136[ 0 ])[ 0 ] * (tmp136[ 0 ])[ 0 ];
      const auto tmp139 = tmp138 + tmp137;
      const auto tmp140 = std::sqrt( tmp139 );
      const auto tmp141 = 4 * tmp140;
      const auto tmp142 = 1 + tmp141;
      const auto tmp143 = std::sqrt( tmp142 );
      const auto tmp144 = 1 + tmp143;
      const auto tmp145 = 2 / tmp144;
      const auto tmp146 = (tmp136[ 0 ])[ 0 ] * tmp145;
      const auto tmp147 = 0.5 * tmp146;
      const auto tmp148 = std::get< 1 >( u );
      const auto tmp149 = (tmp148[ 0 ])[ 1 ] * (tmp148[ 0 ])[ 1 ];
      const auto tmp150 = (tmp148[ 0 ])[ 0 ] * (tmp148[ 0 ])[ 0 ];
      const auto tmp151 = tmp150 + tmp149;
      const auto tmp152 = std::sqrt( tmp151 );
      const auto tmp153 = 4 * tmp152;
      const auto tmp154 = 1 + tmp153;
      const auto tmp155 = std::sqrt( tmp154 );
      const auto tmp156 = 1 + tmp155;
      const auto tmp157 = 2 / tmp156;
      const auto tmp158 = (tmp148[ 0 ])[ 0 ] * tmp157;
      const auto tmp159 = 0.5 * tmp158;
      const auto tmp160 = tmp159 + tmp147;
      const auto tmp161 = (tmp136[ 0 ])[ 1 ] * tmp145;
      const auto tmp162 = 0.5 * tmp161;
      const auto tmp163 = (tmp148[ 0 ])[ 1 ] * tmp157;
      const auto tmp164 = 0.5 * tmp163;
      const auto tmp165 = tmp164 + tmp162;
      return RangeValueType{ { tmp135 }, { { tmp160, tmp165 } } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      using std::sqrt;
      double tmp0 = constant< 0 >();
      const auto tmp1 = 1 / tmp0;
      const auto tmp2 = std::get< 1 >( u );
      const auto tmp3 = (tmp2[ 0 ])[ 1 ] * (tmp2[ 0 ])[ 1 ];
      const auto tmp4 = (tmp2[ 0 ])[ 0 ] * (tmp2[ 0 ])[ 0 ];
      const auto tmp5 = tmp4 + tmp3;
      const auto tmp6 = std::sqrt( tmp5 );
      const auto tmp7 = 4 * tmp6;
      const auto tmp8 = 1 + tmp7;
      const auto tmp9 = std::sqrt( tmp8 );
      const auto tmp10 = 1 + tmp9;
      const auto tmp11 = 2 / tmp10;
      const auto tmp12 = 2 * tmp9;
      const auto tmp13 = 2 * tmp6;
      const auto tmp14 = (tmp2[ 0 ])[ 0 ] + (tmp2[ 0 ])[ 0 ];
      const auto tmp15 = tmp14 / tmp13;
      const auto tmp16 = 4 * tmp15;
      const auto tmp17 = tmp16 / tmp12;
      const auto tmp18 = tmp11 * tmp17;
      const auto tmp19 = -1 * tmp18;
      const auto tmp20 = tmp19 / tmp10;
      const auto tmp21 = (tmp2[ 0 ])[ 0 ] * tmp20;
      const auto tmp22 = tmp21 + tmp11;
      const auto tmp23 = 0.5 * tmp22;
      const auto tmp24 = (tmp2[ 0 ])[ 1 ] + (tmp2[ 0 ])[ 1 ];
      const auto tmp25 = tmp24 / tmp13;
      const auto tmp26 = 4 * tmp25;
      const auto tmp27 = tmp26 / tmp12;
      const auto tmp28 = tmp11 * tmp27;
      const auto tmp29 = -1 * tmp28;
      const auto tmp30 = tmp29 / tmp10;
      const auto tmp31 = (tmp2[ 0 ])[ 0 ] * tmp30;
      const auto tmp32 = 0.5 * tmp31;
      const auto tmp33 = (tmp2[ 0 ])[ 1 ] * tmp20;
      const auto tmp34 = 0.5 * tmp33;
      const auto tmp35 = (tmp2[ 0 ])[ 1 ] * tmp30;
      const auto tmp36 = tmp35 + tmp11;
      const auto tmp37 = 0.5 * tmp36;
      return [ tmp1, tmp23, tmp32, tmp34, tmp37 ] ( const DomainValueType &phi ) {
          return RangeValueType{ { tmp1 * (std::get< 0 >( phi ))[ 0 ] }, { { tmp23 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ] + tmp32 * ((std::get< 1 >( phi ))[ 0 ])[ 1 ], tmp34 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ] + tmp37 * ((std::get< 1 >( phi ))[ 0 ])[ 1 ] } } };
        };
    }

    template< class Point >
    RangeValueType boundary ( const Point &x, const DomainValueType &u ) const
    {
      using std::exp;
      using std::pow;
      using std::sqrt;
      double tmp0 = constant< 1 >();
      double tmp1 = constant< 0 >();
      const auto tmp2 = tmp1 + tmp0;
      const auto tmp3 = -2 * tmp2;
      const auto tmp4 = std::exp( tmp3 );
      GlobalCoordinateType tmp5 = geometry().global( Dune::Fem::coordinate( x ) );
      const auto tmp6 = std::pow( tmp5[ 0 ], 2 );
      const auto tmp7 = 3 * tmp6;
      const auto tmp8 = 0.3333333333333333 * tmp7;
      const auto tmp9 = -1 * tmp8;
      const auto tmp10 = 2 * tmp5[ 0 ];
      const auto tmp11 = 0.5 * tmp10;
      const auto tmp12 = tmp11 + tmp9;
      const auto tmp13 = tmp12 * tmp4;
      const auto tmp14 = tmp13 * tmp13;
      const auto tmp15 = std::pow( tmp5[ 1 ], 2 );
      const auto tmp16 = 3 * tmp15;
      const auto tmp17 = -1 * tmp16;
      const auto tmp18 = 0.3333333333333333 * tmp17;
      const auto tmp19 = -1 * tmp18;
      const auto tmp20 = 2 * tmp5[ 1 ];
      const auto tmp21 = 0.5 * tmp20;
      const auto tmp22 = tmp21 + tmp19;
      const auto tmp23 = tmp22 * tmp4;
      const auto tmp24 = tmp23 * tmp23;
      const auto tmp25 = tmp24 + tmp14;
      const auto tmp26 = std::sqrt( tmp25 );
      const auto tmp27 = 4 * tmp26;
      const auto tmp28 = 1 + tmp27;
      const auto tmp29 = std::sqrt( tmp28 );
      const auto tmp30 = 1 + tmp29;
      const auto tmp31 = 2 / tmp30;
      const auto tmp32 = tmp13 * tmp31;
      const auto tmp33 = -2 * tmp0;
      const auto tmp34 = std::exp( tmp33 );
      const auto tmp35 = tmp12 * tmp34;
      const auto tmp36 = tmp35 * tmp35;
      const auto tmp37 = tmp22 * tmp34;
      const auto tmp38 = tmp37 * tmp37;
      const auto tmp39 = tmp38 + tmp36;
      const auto tmp40 = std::sqrt( tmp39 );
      const auto tmp41 = 4 * tmp40;
      const auto tmp42 = 1 + tmp41;
      const auto tmp43 = std::sqrt( tmp42 );
      const auto tmp44 = 1 + tmp43;
      const auto tmp45 = 2 / tmp44;
      const auto tmp46 = tmp35 * tmp45;
      const auto tmp47 = tmp46 + tmp32;
      GlobalCoordinateType tmp48 = intersection_.unitOuterNormal( x.localPosition() );
      const auto tmp49 = tmp48[ 0 ] * tmp47;
      const auto tmp50 = tmp23 * tmp31;
      const auto tmp51 = tmp37 * tmp45;
      const auto tmp52 = tmp51 + tmp50;
      const auto tmp53 = tmp48[ 1 ] * tmp52;
      const auto tmp54 = tmp53 + tmp49;
      const auto tmp55 = 0.5 * tmp54;
      const auto tmp56 = -1 * tmp55;
      return RangeValueType{ { tmp56 }, { { 0, 0 } } };
    }

    template< class Point >
    auto linearizedBoundary ( const Point &x, const DomainValueType &u ) const
    {
      return [] ( const DomainValueType &phi ) {
          return RangeValueType{ { 0 }, { { 0, 0 } } };
        };
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

    const Condt &condt () const
    {
      return *std::get< 0 >( constants_ );
    }

    Condt &condt ()
    {
      return *std::get< 0 >( constants_ );
    }

    const Cont &cont () const
    {
      return *std::get< 1 >( constants_ );
    }

    Cont &cont ()
    {
      return *std::get< 1 >( constants_ );
    }

    const EntityType &entity () const
    {
      return entity_;
    }

    const Geometry &geometry () const
    {
      return *geometry_;
    }

    template< std::size_t i, class Point >
    typename CoefficientFunctionSpaceType< i >::RangeType evaluateCoefficient ( const Point &x ) const
    {
      typename CoefficientFunctionSpaceType< i >::RangeType result;
      std::get< i >( coefficients_ ).evaluate( x, result );;
      return result;
    }

    template< std::size_t i, class Point >
    typename CoefficientFunctionSpaceType< i >::JacobianRangeType jacobianCoefficient ( const Point &x ) const
    {
      typename CoefficientFunctionSpaceType< i >::JacobianRangeType result;
      std::get< i >( coefficients_ ).jacobian( x, result );;
      return result;
    }

    template< std::size_t i, class Point >
    typename CoefficientFunctionSpaceType< i >::HessianRangeType hessianCoefficient ( const Point &x ) const
    {
      typename CoefficientFunctionSpaceType< i >::HessianRangeType result;
      std::get< i >( coefficients_ ).hessian( x, result );;
      return result;
    }
    EntityType entity_;
    IntersectionType intersection_;
    std::optional< Geometry > geometry_;
    ConstantTupleType constants_;
    std::tuple< Dune::Fem::ConstLocalFunction< Coeffprevious > > coefficients_;
  };

} // namespace Integrands_359275bea09cbad4b58d9f60b3efa24bv1_3

PYBIND11_MODULE( integrands_359275bea09cbad4b58d9f60b3efa24bv1_3_19659fe2f2f74bcfb5d532f1a0d38be9, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView > GridPart;
  typedef Integrands_359275bea09cbad4b58d9f60b3efa24bv1_3::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_359275bea09cbad4b58d9f60b3efa24bv1_3::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_359275bea09cbad4b58d9f60b3efa24bv1_3_19659fe2f2f74bcfb5d532f1a0d38be9.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] ( const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffprevious ) { return new Integrands( coeffprevious ); } ), pybind11::keep_alive< 1, 2 >() );
      cls.def_property( "dt", [] ( Integrands &self ) -> Integrands::Condt { return self.condt(); }, [] ( Integrands &self, const Integrands::Condt &v ) { self.condt() = v; } );
      cls.def_property( "t", [] ( Integrands &self ) -> Integrands::Cont { return self.cont(); }, [] ( Integrands &self, const Integrands::Cont &v ) { self.cont() = v; } );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return false;});
  }
}
#endif // GuardIntegrands_359275bea09cbad4b58d9f60b3efa24bv1_3
