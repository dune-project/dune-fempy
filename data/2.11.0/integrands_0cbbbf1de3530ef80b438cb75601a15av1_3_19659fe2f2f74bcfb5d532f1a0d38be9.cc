#ifndef GuardIntegrands_0cbbbf1de3530ef80b438cb75601a15av1_3
#define GuardIntegrands_0cbbbf1de3530ef80b438cb75601a15av1_3
#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/python/grid/hierarchical.hh>
#include <tuple>
#include <utility>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/extensions.h>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/geometry/edgelength.hh>
#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/function/simplegridfunction.hh>
#include <dune/fem/misc/gridfunctionview.hh>
#include <dune/fempy/py/integrands.hh>

namespace Integrands_0cbbbf1de3530ef80b438cb75601a15av1_3
{




  // Integrands
  // ----------

  template< class GridPart, class Coeffanalytical_solution_force_term, class Coeffanalytical_solution_x, class Coeffanalytical_solution_y >
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
    typedef double Coninterior_penalty_parameter;
    typedef double Conswitch_interior_penalty_method;
    typedef double Congrad_div_stab_param;
    typedef double Conmass_flux_stab_param;
    typedef double Constab_cburn;
    typedef double Conall_boundaries;
    typedef std::tuple< std::shared_ptr< Coninterior_penalty_parameter >, std::shared_ptr< Conswitch_interior_penalty_method >, std::shared_ptr< Congrad_div_stab_param >, std::shared_ptr< Conmass_flux_stab_param >, std::shared_ptr< Constab_cburn >, std::shared_ptr< Conall_boundaries > > ConstantTupleType;
    template< std::size_t i >
    using ConstantsRangeType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;
    typedef std::tuple< Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 2 > >, Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > >, Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > > > CoefficientFunctionSpaceTupleType;
    typedef std::tuple< Coeffanalytical_solution_force_term, Coeffanalytical_solution_x, Coeffanalytical_solution_y > CoefficientTupleType;
    template< std::size_t i >
    using CoefficientFunctionSpaceType = std::tuple_element_t< i, CoefficientFunctionSpaceTupleType >;
    template< std::size_t i >
    using CoefficientRangeType = typename CoefficientFunctionSpaceType< i >::RangeType;
    template< std::size_t i >
    using CoefficientJacobianRangeType = typename CoefficientFunctionSpaceType< i >::JacobianRangeType;
    static constexpr bool gridPartValid = Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffanalytical_solution_force_term>>() && Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffanalytical_solution_x>>() && Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffanalytical_solution_y>>();
    template< std::size_t i >
    using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
    template< std::size_t i >
    using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;

    Integrands ( const Coeffanalytical_solution_force_term &coeffanalytical_solution_force_term, const Coeffanalytical_solution_x &coeffanalytical_solution_x, const Coeffanalytical_solution_y &coeffanalytical_solution_y, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
      : coefficients_{{ std::tuple< Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_force_term >, Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_x >, Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_y > >( Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_force_term >( coeffanalytical_solution_force_term ), Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_x >( coeffanalytical_solution_x ), Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_y >( coeffanalytical_solution_y ) ), std::tuple< Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_force_term >, Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_x >, Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_y > >( Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_force_term >( coeffanalytical_solution_force_term ), Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_x >( coeffanalytical_solution_x ), Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_y >( coeffanalytical_solution_y ) ) }}
    {
      std::get< 0 >( constants_ ) = std::make_shared< Coninterior_penalty_parameter >( (Coninterior_penalty_parameter(0)) );
      std::get< 1 >( constants_ ) = std::make_shared< Conswitch_interior_penalty_method >( (Conswitch_interior_penalty_method(0)) );
      std::get< 2 >( constants_ ) = std::make_shared< Congrad_div_stab_param >( (Congrad_div_stab_param(0)) );
      std::get< 3 >( constants_ ) = std::make_shared< Conmass_flux_stab_param >( (Conmass_flux_stab_param(0)) );
      std::get< 4 >( constants_ ) = std::make_shared< Constab_cburn >( (Constab_cburn(0)) );
      std::get< 5 >( constants_ ) = std::make_shared< Conall_boundaries >( (Conall_boundaries(0)) );
    }

    bool init ( const EntityType &entity )
    {
      entity_[ (static_cast< std::size_t >( Side::in )) ] = entity;
      geometry_.emplace( this->entity().geometry() );
      cellVolume_[ (static_cast< std::size_t >( Side::in )) ] = (geometry().volume());
      std::get< 0 >( coefficients_[ static_cast< std::size_t >( Side::in ) ] ).bind( this->entity() );
      std::get< 0 >( coefficients_[ static_cast< std::size_t >( Side::out ) ] ).bind( this->entity() );
      std::get< 1 >( coefficients_[ static_cast< std::size_t >( Side::in ) ] ).bind( this->entity() );
      std::get< 1 >( coefficients_[ static_cast< std::size_t >( Side::out ) ] ).bind( this->entity() );
      std::get< 2 >( coefficients_[ static_cast< std::size_t >( Side::in ) ] ).bind( this->entity() );
      std::get< 2 >( coefficients_[ static_cast< std::size_t >( Side::out ) ] ).bind( this->entity() );
      return true;
    }

    void unbind ()
    {
      std::get< 0 >( coefficients_[ static_cast< std::size_t >( Side::in ) ] ).unbind( );
      std::get< 0 >( coefficients_[ static_cast< std::size_t >( Side::out ) ] ).unbind( );
      std::get< 1 >( coefficients_[ static_cast< std::size_t >( Side::in ) ] ).unbind( );
      std::get< 1 >( coefficients_[ static_cast< std::size_t >( Side::out ) ] ).unbind( );
      std::get< 2 >( coefficients_[ static_cast< std::size_t >( Side::in ) ] ).unbind( );
      std::get< 2 >( coefficients_[ static_cast< std::size_t >( Side::out ) ] ).unbind( );
    }

    bool init ( const IntersectionType &intersection )
    {
      intersection_ = intersection;
      const auto& intersectionGeom = intersection_.geometry();
      facetArea_ = (intersectionGeom.volume());
      entity_[ (static_cast< std::size_t >( Side::in )) ] = (intersection.inside());
      geometry_.emplace( this->entity().geometry() );
      cellVolume_[ (static_cast< std::size_t >( Side::in )) ] = (geometry().volume());
      std::get< 0 >( coefficients_[ static_cast< std::size_t >( Side::in ) ] ).bind( intersection_, Side::in  );
      std::get< 1 >( coefficients_[ static_cast< std::size_t >( Side::in ) ] ).bind( intersection_, Side::in  );
      std::get< 2 >( coefficients_[ static_cast< std::size_t >( Side::in ) ] ).bind( intersection_, Side::in  );
      if( intersection.neighbor() )
      {
        entity_[ static_cast< std::size_t >( Side::out ) ] = intersection.outside();
        const Geometry& outsideGeom = entity_[ static_cast< std::size_t >( Side::out ) ].geometry();
        cellVolume_[ static_cast< std::size_t >( Side::out ) ] = outsideGeom.volume();
      std::get< 0 >( coefficients_[ static_cast< std::size_t >( Side::out ) ] ).bind( intersection_, Side::out );
      std::get< 1 >( coefficients_[ static_cast< std::size_t >( Side::out ) ] ).bind( intersection_, Side::out );
      std::get< 2 >( coefficients_[ static_cast< std::size_t >( Side::out ) ] ).bind( intersection_, Side::out );
      }
      return true;
    }
    typedef std::tuple< Dune::FieldVector< double, 3 >, Dune::FieldMatrix< double, 3, 2 > > DomainValueType;
    typedef std::tuple< Dune::FieldVector< double, 3 >, Dune::FieldMatrix< double, 3, 2 > > RangeValueType;
    static constexpr bool _nonlinear = false;

    bool nonlinear () const
    {
      return _nonlinear;
    }

    template< class Point >
    RangeValueType interior ( const Point &x, const DomainValueType &u ) const
    {
      typename CoefficientFunctionSpaceType< 0 >::RangeType tmp0 = evaluateCoefficient< 0 >( x );
      const auto tmp1 = -1 * tmp0[ 0 ];
      const auto tmp2 = -1 * tmp0[ 1 ];
      const auto tmp3 = std::get< 1 >( u );
      const auto tmp4 = (tmp3[ 0 ])[ 0 ] + (tmp3[ 1 ])[ 1 ];
      const auto tmp5 = -1 * tmp4;
      double tmp6 = constant< 2 >();
      const auto tmp7 = tmp6 * tmp4;
      const auto tmp8 = std::get< 0 >( u );
      const auto tmp9 = -1 * tmp8[ 2 ];
      const auto tmp10 = (tmp3[ 0 ])[ 0 ] + tmp9;
      const auto tmp11 = tmp10 + tmp7;
      const auto tmp12 = (tmp3[ 1 ])[ 1 ] + tmp9;
      const auto tmp13 = tmp12 + tmp7;
      return RangeValueType{ { tmp1, tmp2, tmp5 }, { { tmp11, (tmp3[ 0 ])[ 1 ] }, { (tmp3[ 1 ])[ 0 ], tmp13 }, { 0.0, 0.0 } } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      double tmp0 = constant< 2 >();
      const auto tmp1 = 1 + tmp0;
      return [ tmp0, tmp1 ] ( const DomainValueType &phi ) {
          return RangeValueType{ { 0, 0, -1 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ] + -1 * ((std::get< 1 >( phi ))[ 1 ])[ 1 ] }, { { (-1 * (std::get< 0 >( phi ))[ 2 ] + tmp1 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ]) + tmp0 * ((std::get< 1 >( phi ))[ 1 ])[ 1 ], ((std::get< 1 >( phi ))[ 0 ])[ 1 ] }, { ((std::get< 1 >( phi ))[ 1 ])[ 0 ], (-1 * (std::get< 0 >( phi ))[ 2 ] + tmp0 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ]) + tmp1 * ((std::get< 1 >( phi ))[ 1 ])[ 1 ] }, { 0, 0 } } };
        };
    }

    template< class Point >
    RangeValueType boundary ( const Point &x, const DomainValueType &u ) const
    {
      auto tmp0 = facetArea_;
      auto tmp1 = cellVolume();
      const auto tmp2 = tmp1 / tmp0;
      const auto tmp3 = 1 / tmp2;
      typename CoefficientFunctionSpaceType< 2 >::RangeType tmp4 = evaluateCoefficient< 2 >( x );
      const auto tmp5 = -1 * tmp4[ 0 ];
      const auto tmp6 = std::get< 0 >( u );
      const auto tmp7 = tmp6[ 1 ] + tmp5;
      GlobalCoordinateType tmp8 = intersection_.unitOuterNormal( x.localPosition() );
      const auto tmp9 = tmp8[ 1 ] * tmp7;
      typename CoefficientFunctionSpaceType< 1 >::RangeType tmp10 = evaluateCoefficient< 1 >( x );
      const auto tmp11 = -1 * tmp10[ 0 ];
      const auto tmp12 = tmp6[ 0 ] + tmp11;
      const auto tmp13 = tmp8[ 0 ] * tmp12;
      const auto tmp14 = tmp13 + tmp9;
      const auto tmp15 = tmp14 * tmp3;
      const auto tmp16 = tmp8[ 0 ] * tmp15;
      double tmp17 = constant< 3 >();
      const auto tmp18 = tmp17 * tmp16;
      const auto tmp19 = tmp8[ 0 ] * tmp6[ 2 ];
      double tmp20 = constant< 0 >();
      const auto tmp21 = tmp20 / tmp2;
      const auto tmp22 = tmp12 * tmp21;
      double tmp23 = constant< 5 >();
      const auto tmp24 = tmp23 * tmp22;
      const auto tmp25 = std::get< 1 >( u );
      const auto tmp26 = (tmp25[ 0 ])[ 1 ] * tmp8[ 1 ];
      const auto tmp27 = tmp8[ 0 ] * (tmp25[ 0 ])[ 0 ];
      const auto tmp28 = tmp27 + tmp26;
      const auto tmp29 = tmp23 * tmp28;
      const auto tmp30 = -1 * tmp29;
      const auto tmp31 = tmp30 + tmp24;
      const auto tmp32 = tmp31 + tmp19;
      const auto tmp33 = tmp32 + tmp18;
      const auto tmp34 = tmp8[ 1 ] * tmp15;
      const auto tmp35 = tmp17 * tmp34;
      const auto tmp36 = tmp8[ 1 ] * tmp6[ 2 ];
      const auto tmp37 = tmp7 * tmp21;
      const auto tmp38 = tmp23 * tmp37;
      const auto tmp39 = tmp8[ 1 ] * (tmp25[ 1 ])[ 1 ];
      const auto tmp40 = tmp8[ 0 ] * (tmp25[ 1 ])[ 0 ];
      const auto tmp41 = tmp40 + tmp39;
      const auto tmp42 = tmp23 * tmp41;
      const auto tmp43 = -1 * tmp42;
      const auto tmp44 = tmp43 + tmp38;
      const auto tmp45 = tmp44 + tmp36;
      const auto tmp46 = tmp45 + tmp35;
      double tmp47 = constant< 1 >();
      const auto tmp48 = tmp47 * tmp13;
      const auto tmp49 = tmp23 * tmp48;
      const auto tmp50 = tmp8[ 1 ] * tmp12;
      const auto tmp51 = tmp47 * tmp50;
      const auto tmp52 = tmp23 * tmp51;
      const auto tmp53 = tmp8[ 0 ] * tmp7;
      const auto tmp54 = tmp47 * tmp53;
      const auto tmp55 = tmp23 * tmp54;
      const auto tmp56 = tmp47 * tmp9;
      const auto tmp57 = tmp23 * tmp56;
      return RangeValueType{ { tmp33, tmp46, tmp14 }, { { tmp49, tmp52 }, { tmp55, tmp57 }, { 0.0, 0.0 } } };
    }

    template< class Point >
    auto linearizedBoundary ( const Point &x, const DomainValueType &u ) const
    {
      auto tmp0 = facetArea_;
      auto tmp1 = cellVolume();
      const auto tmp2 = tmp1 / tmp0;
      double tmp3 = constant< 0 >();
      const auto tmp4 = tmp3 / tmp2;
      double tmp5 = constant< 5 >();
      const auto tmp6 = tmp5 * tmp4;
      const auto tmp7 = 1 / tmp2;
      GlobalCoordinateType tmp8 = intersection_.unitOuterNormal( x.localPosition() );
      const auto tmp9 = tmp8[ 0 ] * tmp7;
      const auto tmp10 = tmp8[ 0 ] * tmp9;
      double tmp11 = constant< 3 >();
      const auto tmp12 = tmp11 * tmp10;
      const auto tmp13 = tmp12 + tmp6;
      const auto tmp14 = tmp8[ 1 ] * tmp7;
      const auto tmp15 = tmp8[ 0 ] * tmp14;
      const auto tmp16 = tmp11 * tmp15;
      const auto tmp17 = tmp8[ 1 ] * tmp9;
      const auto tmp18 = tmp11 * tmp17;
      const auto tmp19 = tmp8[ 1 ] * tmp14;
      const auto tmp20 = tmp11 * tmp19;
      const auto tmp21 = tmp20 + tmp6;
      const auto tmp22 = tmp8[ 0 ] * tmp5;
      const auto tmp23 = -1 * tmp22;
      const auto tmp24 = tmp8[ 1 ] * tmp5;
      const auto tmp25 = -1 * tmp24;
      double tmp26 = constant< 1 >();
      const auto tmp27 = tmp8[ 0 ] * tmp26;
      const auto tmp28 = tmp5 * tmp27;
      const auto tmp29 = tmp8[ 1 ] * tmp26;
      const auto tmp30 = tmp5 * tmp29;
      return [ tmp13, tmp16, tmp18, tmp21, tmp23, tmp25, tmp28, tmp30, tmp8 ] ( const DomainValueType &phi ) {
          return RangeValueType{ { (((tmp13 * (std::get< 0 >( phi ))[ 0 ] + tmp16 * (std::get< 0 >( phi ))[ 1 ]) + tmp8[ 0 ] * (std::get< 0 >( phi ))[ 2 ]) + tmp23 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ]) + tmp25 * ((std::get< 1 >( phi ))[ 0 ])[ 1 ], (((tmp18 * (std::get< 0 >( phi ))[ 0 ] + tmp21 * (std::get< 0 >( phi ))[ 1 ]) + tmp8[ 1 ] * (std::get< 0 >( phi ))[ 2 ]) + tmp23 * ((std::get< 1 >( phi ))[ 1 ])[ 0 ]) + tmp25 * ((std::get< 1 >( phi ))[ 1 ])[ 1 ], tmp8[ 0 ] * (std::get< 0 >( phi ))[ 0 ] + tmp8[ 1 ] * (std::get< 0 >( phi ))[ 1 ] }, { { tmp28 * (std::get< 0 >( phi ))[ 0 ], tmp30 * (std::get< 0 >( phi ))[ 0 ] }, { tmp28 * (std::get< 0 >( phi ))[ 1 ], tmp30 * (std::get< 0 >( phi ))[ 1 ] }, { 0, 0 } } };
        };
    }

    template< class Point >
    std::pair< RangeValueType, RangeValueType > skeleton ( const Point &xIn, const DomainValueType &uIn, const Point &xOut, const DomainValueType &uOut ) const
    {
      auto tmp0 = facetArea_;
      auto tmp1 = cellVolume_[ static_cast< std::size_t >( Side::out ) ];
      auto tmp2 = cellVolume_[ static_cast< std::size_t >( Side::in ) ];
      const auto tmp3 = tmp2 + tmp1;
      const auto tmp4 = 0.5 * tmp3;
      const auto tmp5 = tmp4 / tmp0;
      const auto tmp6 = 1 / tmp5;
      const auto tmp7 = std::get< 0 >( uOut );
      const auto tmp8 = -1 * tmp7[ 1 ];
      const auto tmp9 = std::get< 0 >( uIn );
      const auto tmp10 = tmp9[ 1 ] + tmp8;
      GlobalCoordinateType tmp11 = intersection_.unitOuterNormal( xIn.localPosition() );
      const auto tmp12 = tmp11[ 1 ] * tmp10;
      const auto tmp13 = -1 * tmp7[ 0 ];
      const auto tmp14 = tmp9[ 0 ] + tmp13;
      const auto tmp15 = tmp11[ 0 ] * tmp14;
      const auto tmp16 = tmp15 + tmp12;
      const auto tmp17 = tmp16 * tmp6;
      const auto tmp18 = tmp11[ 0 ] * tmp17;
      double tmp19 = constant< 3 >();
      const auto tmp20 = tmp19 * tmp18;
      const auto tmp21 = tmp9[ 2 ] + tmp7[ 2 ];
      const auto tmp22 = 0.5 * tmp21;
      const auto tmp23 = tmp11[ 0 ] * tmp22;
      double tmp24 = constant< 0 >();
      const auto tmp25 = tmp24 / tmp5;
      const auto tmp26 = tmp14 * tmp25;
      const auto tmp27 = std::get< 1 >( uOut );
      const auto tmp28 = std::get< 1 >( uIn );
      const auto tmp29 = (tmp28[ 0 ])[ 1 ] + (tmp27[ 0 ])[ 1 ];
      const auto tmp30 = 0.5 * tmp29;
      const auto tmp31 = tmp11[ 1 ] * tmp30;
      const auto tmp32 = (tmp28[ 0 ])[ 0 ] + (tmp27[ 0 ])[ 0 ];
      const auto tmp33 = 0.5 * tmp32;
      const auto tmp34 = tmp11[ 0 ] * tmp33;
      const auto tmp35 = tmp34 + tmp31;
      const auto tmp36 = -1 * tmp35;
      const auto tmp37 = tmp36 + tmp26;
      const auto tmp38 = tmp37 + tmp23;
      const auto tmp39 = tmp38 + tmp20;
      const auto tmp40 = tmp11[ 1 ] * tmp17;
      const auto tmp41 = tmp19 * tmp40;
      const auto tmp42 = tmp11[ 1 ] * tmp22;
      const auto tmp43 = tmp10 * tmp25;
      const auto tmp44 = (tmp28[ 1 ])[ 1 ] + (tmp27[ 1 ])[ 1 ];
      const auto tmp45 = 0.5 * tmp44;
      const auto tmp46 = tmp11[ 1 ] * tmp45;
      const auto tmp47 = (tmp28[ 1 ])[ 0 ] + (tmp27[ 1 ])[ 0 ];
      const auto tmp48 = 0.5 * tmp47;
      const auto tmp49 = tmp11[ 0 ] * tmp48;
      const auto tmp50 = tmp49 + tmp46;
      const auto tmp51 = -1 * tmp50;
      const auto tmp52 = tmp51 + tmp43;
      const auto tmp53 = tmp52 + tmp42;
      const auto tmp54 = tmp53 + tmp41;
      double tmp55 = constant< 4 >();
      const auto tmp56 = tmp55 * tmp5;
      const auto tmp57 = -1 * tmp7[ 2 ];
      const auto tmp58 = tmp9[ 2 ] + tmp57;
      const auto tmp59 = tmp11[ 1 ] * tmp58;
      const auto tmp60 = tmp11[ 1 ] * tmp59;
      const auto tmp61 = tmp11[ 0 ] * tmp58;
      const auto tmp62 = tmp11[ 0 ] * tmp61;
      const auto tmp63 = tmp62 + tmp60;
      const auto tmp64 = tmp63 * tmp56;
      const auto tmp65 = 0.5 * tmp16;
      const auto tmp66 = tmp65 + tmp64;
      const auto tmp67 = 0.5 * tmp11[ 0 ];
      const auto tmp68 = tmp14 * tmp67;
      double tmp69 = constant< 1 >();
      const auto tmp70 = tmp69 * tmp68;
      const auto tmp71 = 0.5 * tmp11[ 1 ];
      const auto tmp72 = tmp14 * tmp71;
      const auto tmp73 = tmp69 * tmp72;
      const auto tmp74 = tmp10 * tmp67;
      const auto tmp75 = tmp69 * tmp74;
      const auto tmp76 = tmp10 * tmp71;
      const auto tmp77 = tmp69 * tmp76;
      const auto tmp78 = -1 * tmp11[ 0 ];
      const auto tmp79 = tmp78 * tmp17;
      const auto tmp80 = tmp19 * tmp79;
      const auto tmp81 = tmp78 * tmp22;
      const auto tmp82 = -1 * tmp14;
      const auto tmp83 = tmp82 * tmp25;
      const auto tmp84 = -1 * tmp36;
      const auto tmp85 = tmp84 + tmp83;
      const auto tmp86 = tmp85 + tmp81;
      const auto tmp87 = tmp86 + tmp80;
      const auto tmp88 = -1 * tmp11[ 1 ];
      const auto tmp89 = tmp88 * tmp17;
      const auto tmp90 = tmp19 * tmp89;
      const auto tmp91 = tmp88 * tmp22;
      const auto tmp92 = -1 * tmp10;
      const auto tmp93 = tmp92 * tmp25;
      const auto tmp94 = -1 * tmp51;
      const auto tmp95 = tmp94 + tmp93;
      const auto tmp96 = tmp95 + tmp91;
      const auto tmp97 = tmp96 + tmp90;
      const auto tmp98 = tmp88 * tmp59;
      const auto tmp99 = tmp78 * tmp61;
      const auto tmp100 = tmp99 + tmp98;
      const auto tmp101 = tmp100 * tmp56;
      const auto tmp102 = tmp65 + tmp101;
      return std::make_pair( RangeValueType{ { tmp39, tmp54, tmp66 }, { { tmp70, tmp73 }, { tmp75, tmp77 }, { 0.0, 0.0 } } }, RangeValueType{ { tmp87, tmp97, tmp102 }, { { tmp70, tmp73 }, { tmp75, tmp77 }, { 0.0, 0.0 } } } );
    }

    template< class Point >
    auto linearizedSkeleton ( const Point &xIn, const DomainValueType &uIn, const Point &xOut, const DomainValueType &uOut ) const
    {
      auto tmp0 = facetArea_;
      auto tmp1 = cellVolume_[ static_cast< std::size_t >( Side::out ) ];
      auto tmp2 = cellVolume_[ static_cast< std::size_t >( Side::in ) ];
      const auto tmp3 = tmp2 + tmp1;
      const auto tmp4 = 0.5 * tmp3;
      const auto tmp5 = tmp4 / tmp0;
      double tmp6 = constant< 0 >();
      const auto tmp7 = tmp6 / tmp5;
      const auto tmp8 = 1 / tmp5;
      GlobalCoordinateType tmp9 = intersection_.unitOuterNormal( xIn.localPosition() );
      const auto tmp10 = tmp9[ 0 ] * tmp8;
      const auto tmp11 = tmp9[ 0 ] * tmp10;
      double tmp12 = constant< 3 >();
      const auto tmp13 = tmp12 * tmp11;
      const auto tmp14 = tmp13 + tmp7;
      const auto tmp15 = tmp9[ 1 ] * tmp8;
      const auto tmp16 = tmp9[ 0 ] * tmp15;
      const auto tmp17 = tmp12 * tmp16;
      const auto tmp18 = 0.5 * tmp9[ 0 ];
      const auto tmp19 = tmp9[ 1 ] * tmp10;
      const auto tmp20 = tmp12 * tmp19;
      const auto tmp21 = tmp9[ 1 ] * tmp15;
      const auto tmp22 = tmp12 * tmp21;
      const auto tmp23 = tmp22 + tmp7;
      const auto tmp24 = 0.5 * tmp9[ 1 ];
      double tmp25 = constant< 4 >();
      const auto tmp26 = tmp25 * tmp5;
      const auto tmp27 = tmp9[ 1 ] * tmp9[ 1 ];
      const auto tmp28 = tmp9[ 0 ] * tmp9[ 0 ];
      const auto tmp29 = tmp28 + tmp27;
      const auto tmp30 = tmp29 * tmp26;
      const auto tmp31 = -1 * tmp18;
      const auto tmp32 = -1 * tmp24;
      double tmp33 = constant< 1 >();
      const auto tmp34 = tmp33 * tmp18;
      const auto tmp35 = tmp33 * tmp24;
      const auto tmp36 = -1 * tmp7;
      const auto tmp37 = -1 * tmp9[ 0 ];
      const auto tmp38 = tmp37 * tmp10;
      const auto tmp39 = tmp12 * tmp38;
      const auto tmp40 = tmp39 + tmp36;
      const auto tmp41 = tmp37 * tmp15;
      const auto tmp42 = tmp12 * tmp41;
      const auto tmp43 = 0.5 * tmp37;
      const auto tmp44 = -1 * tmp9[ 1 ];
      const auto tmp45 = tmp44 * tmp10;
      const auto tmp46 = tmp12 * tmp45;
      const auto tmp47 = tmp44 * tmp15;
      const auto tmp48 = tmp12 * tmp47;
      const auto tmp49 = tmp48 + tmp36;
      const auto tmp50 = 0.5 * tmp44;
      const auto tmp51 = tmp9[ 1 ] * tmp44;
      const auto tmp52 = tmp9[ 0 ] * tmp37;
      const auto tmp53 = tmp52 + tmp51;
      const auto tmp54 = tmp53 * tmp26;
      const auto tmp55 = -1 * tmp31;
      const auto tmp56 = -1 * tmp32;
      const auto tmp57 = tmp37 * tmp8;
      const auto tmp58 = tmp9[ 0 ] * tmp57;
      const auto tmp59 = tmp12 * tmp58;
      const auto tmp60 = tmp59 + tmp36;
      const auto tmp61 = tmp44 * tmp8;
      const auto tmp62 = tmp9[ 0 ] * tmp61;
      const auto tmp63 = tmp12 * tmp62;
      const auto tmp64 = tmp9[ 1 ] * tmp57;
      const auto tmp65 = tmp12 * tmp64;
      const auto tmp66 = tmp9[ 1 ] * tmp61;
      const auto tmp67 = tmp12 * tmp66;
      const auto tmp68 = tmp67 + tmp36;
      const auto tmp69 = tmp33 * tmp31;
      const auto tmp70 = tmp33 * tmp32;
      const auto tmp71 = tmp37 * tmp57;
      const auto tmp72 = tmp12 * tmp71;
      const auto tmp73 = tmp72 + tmp7;
      const auto tmp74 = tmp37 * tmp61;
      const auto tmp75 = tmp12 * tmp74;
      const auto tmp76 = tmp44 * tmp57;
      const auto tmp77 = tmp12 * tmp76;
      const auto tmp78 = tmp44 * tmp61;
      const auto tmp79 = tmp12 * tmp78;
      const auto tmp80 = tmp79 + tmp7;
      const auto tmp81 = tmp44 * tmp44;
      const auto tmp82 = tmp37 * tmp37;
      const auto tmp83 = tmp82 + tmp81;
      const auto tmp84 = tmp83 * tmp26;
      return std::make_pair( [ tmp14, tmp17, tmp18, tmp20, tmp23, tmp24, tmp30, tmp31, tmp32, tmp34, tmp35, tmp40, tmp42, tmp43, tmp46, tmp49, tmp50, tmp54, tmp55, tmp56 ] ( const DomainValueType &phiIn ) {
          return std::make_pair( RangeValueType{ { (((tmp14 * (std::get< 0 >( phiIn ))[ 0 ] + tmp17 * (std::get< 0 >( phiIn ))[ 1 ]) + tmp18 * (std::get< 0 >( phiIn ))[ 2 ]) + tmp31 * ((std::get< 1 >( phiIn ))[ 0 ])[ 0 ]) + tmp32 * ((std::get< 1 >( phiIn ))[ 0 ])[ 1 ], (((tmp20 * (std::get< 0 >( phiIn ))[ 0 ] + tmp23 * (std::get< 0 >( phiIn ))[ 1 ]) + tmp24 * (std::get< 0 >( phiIn ))[ 2 ]) + tmp31 * ((std::get< 1 >( phiIn ))[ 1 ])[ 0 ]) + tmp32 * ((std::get< 1 >( phiIn ))[ 1 ])[ 1 ], (tmp18 * (std::get< 0 >( phiIn ))[ 0 ] + tmp24 * (std::get< 0 >( phiIn ))[ 1 ]) + tmp30 * (std::get< 0 >( phiIn ))[ 2 ] }, { { tmp34 * (std::get< 0 >( phiIn ))[ 0 ], tmp35 * (std::get< 0 >( phiIn ))[ 0 ] }, { tmp34 * (std::get< 0 >( phiIn ))[ 1 ], tmp35 * (std::get< 0 >( phiIn ))[ 1 ] }, { 0, 0 } } }, RangeValueType{ { (((tmp40 * (std::get< 0 >( phiIn ))[ 0 ] + tmp42 * (std::get< 0 >( phiIn ))[ 1 ]) + tmp43 * (std::get< 0 >( phiIn ))[ 2 ]) + tmp55 * ((std::get< 1 >( phiIn ))[ 0 ])[ 0 ]) + tmp56 * ((std::get< 1 >( phiIn ))[ 0 ])[ 1 ], (((tmp46 * (std::get< 0 >( phiIn ))[ 0 ] + tmp49 * (std::get< 0 >( phiIn ))[ 1 ]) + tmp50 * (std::get< 0 >( phiIn ))[ 2 ]) + tmp55 * ((std::get< 1 >( phiIn ))[ 1 ])[ 0 ]) + tmp56 * ((std::get< 1 >( phiIn ))[ 1 ])[ 1 ], (tmp18 * (std::get< 0 >( phiIn ))[ 0 ] + tmp24 * (std::get< 0 >( phiIn ))[ 1 ]) + tmp54 * (std::get< 0 >( phiIn ))[ 2 ] }, { { tmp34 * (std::get< 0 >( phiIn ))[ 0 ], tmp35 * (std::get< 0 >( phiIn ))[ 0 ] }, { tmp34 * (std::get< 0 >( phiIn ))[ 1 ], tmp35 * (std::get< 0 >( phiIn ))[ 1 ] }, { 0, 0 } } } );
        }, [ tmp18, tmp24, tmp31, tmp32, tmp43, tmp50, tmp54, tmp55, tmp56, tmp60, tmp63, tmp65, tmp68, tmp69, tmp70, tmp73, tmp75, tmp77, tmp80, tmp84 ] ( const DomainValueType &phiOut ) {
          return std::make_pair( RangeValueType{ { (((tmp60 * (std::get< 0 >( phiOut ))[ 0 ] + tmp63 * (std::get< 0 >( phiOut ))[ 1 ]) + tmp18 * (std::get< 0 >( phiOut ))[ 2 ]) + tmp31 * ((std::get< 1 >( phiOut ))[ 0 ])[ 0 ]) + tmp32 * ((std::get< 1 >( phiOut ))[ 0 ])[ 1 ], (((tmp65 * (std::get< 0 >( phiOut ))[ 0 ] + tmp68 * (std::get< 0 >( phiOut ))[ 1 ]) + tmp24 * (std::get< 0 >( phiOut ))[ 2 ]) + tmp31 * ((std::get< 1 >( phiOut ))[ 1 ])[ 0 ]) + tmp32 * ((std::get< 1 >( phiOut ))[ 1 ])[ 1 ], (tmp43 * (std::get< 0 >( phiOut ))[ 0 ] + tmp50 * (std::get< 0 >( phiOut ))[ 1 ]) + tmp54 * (std::get< 0 >( phiOut ))[ 2 ] }, { { tmp69 * (std::get< 0 >( phiOut ))[ 0 ], tmp70 * (std::get< 0 >( phiOut ))[ 0 ] }, { tmp69 * (std::get< 0 >( phiOut ))[ 1 ], tmp70 * (std::get< 0 >( phiOut ))[ 1 ] }, { 0, 0 } } }, RangeValueType{ { (((tmp73 * (std::get< 0 >( phiOut ))[ 0 ] + tmp75 * (std::get< 0 >( phiOut ))[ 1 ]) + tmp43 * (std::get< 0 >( phiOut ))[ 2 ]) + tmp55 * ((std::get< 1 >( phiOut ))[ 0 ])[ 0 ]) + tmp56 * ((std::get< 1 >( phiOut ))[ 0 ])[ 1 ], (((tmp77 * (std::get< 0 >( phiOut ))[ 0 ] + tmp80 * (std::get< 0 >( phiOut ))[ 1 ]) + tmp50 * (std::get< 0 >( phiOut ))[ 2 ]) + tmp55 * ((std::get< 1 >( phiOut ))[ 1 ])[ 0 ]) + tmp56 * ((std::get< 1 >( phiOut ))[ 1 ])[ 1 ], (tmp43 * (std::get< 0 >( phiOut ))[ 0 ] + tmp50 * (std::get< 0 >( phiOut ))[ 1 ]) + tmp84 * (std::get< 0 >( phiOut ))[ 2 ] }, { { tmp69 * (std::get< 0 >( phiOut ))[ 0 ], tmp70 * (std::get< 0 >( phiOut ))[ 0 ] }, { tmp69 * (std::get< 0 >( phiOut ))[ 1 ], tmp70 * (std::get< 0 >( phiOut ))[ 1 ] }, { 0, 0 } } } );
        } );
    }
    typedef Dune::FieldVector< double, 3 >  RRangeType;
    typedef Dune::FieldMatrix< double, 3, GridPartType::dimension >  RJacobianRangeType;
    typedef Dune::Fem::BoundaryIdProvider< typename GridPartType::GridType > BoundaryIdProviderType;
    typedef std::array<int,3> DirichletComponentType;

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

    const Coninterior_penalty_parameter &coninterior_penalty_parameter () const
    {
      return *std::get< 0 >( constants_ );
    }

    Coninterior_penalty_parameter &coninterior_penalty_parameter ()
    {
      return *std::get< 0 >( constants_ );
    }

    const Conswitch_interior_penalty_method &conswitch_interior_penalty_method () const
    {
      return *std::get< 1 >( constants_ );
    }

    Conswitch_interior_penalty_method &conswitch_interior_penalty_method ()
    {
      return *std::get< 1 >( constants_ );
    }

    const Congrad_div_stab_param &congrad_div_stab_param () const
    {
      return *std::get< 2 >( constants_ );
    }

    Congrad_div_stab_param &congrad_div_stab_param ()
    {
      return *std::get< 2 >( constants_ );
    }

    const Conmass_flux_stab_param &conmass_flux_stab_param () const
    {
      return *std::get< 3 >( constants_ );
    }

    Conmass_flux_stab_param &conmass_flux_stab_param ()
    {
      return *std::get< 3 >( constants_ );
    }

    const Constab_cburn &constab_cburn () const
    {
      return *std::get< 4 >( constants_ );
    }

    Constab_cburn &constab_cburn ()
    {
      return *std::get< 4 >( constants_ );
    }

    const Conall_boundaries &conall_boundaries () const
    {
      return *std::get< 5 >( constants_ );
    }

    Conall_boundaries &conall_boundaries ()
    {
      return *std::get< 5 >( constants_ );
    }

    const EntityType &entity () const
    {
      return entity_[ (static_cast< std::size_t >( Side::in )) ];
    }

    const Geometry &geometry () const
    {
      return *geometry_;
    }

    const ctype cellVolume () const
    {
      return cellVolume_[ (static_cast< std::size_t >( Side::in )) ];
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
    std::optional< Geometry > geometry_;
    std::array< ctype, 2 > cellVolume_;
    ctype facetArea_;
    ConstantTupleType constants_;
    std::array< std::tuple< Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_force_term >, Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_x >, Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_y > >, 2 > coefficients_;
  };

} // namespace Integrands_0cbbbf1de3530ef80b438cb75601a15av1_3

PYBIND11_MODULE( integrands_0cbbbf1de3530ef80b438cb75601a15av1_3_19659fe2f2f74bcfb5d532f1a0d38be9, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView > GridPart;
  typedef Integrands_0cbbbf1de3530ef80b438cb75601a15av1_3::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 2 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_0cbbbf1de3530ef80b438cb75601a15av1_3::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 2 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_0cbbbf1de3530ef80b438cb75601a15av1_3_19659fe2f2f74bcfb5d532f1a0d38be9.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] ( const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 2 > > &coeffanalytical_solution_force_term, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffanalytical_solution_x, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffanalytical_solution_y ) { return new Integrands( coeffanalytical_solution_force_term, coeffanalytical_solution_x, coeffanalytical_solution_y ); } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 4 >() );
      cls.def_property( "interior_penalty_parameter", [] ( Integrands &self ) -> Integrands::Coninterior_penalty_parameter { return self.coninterior_penalty_parameter(); }, [] ( Integrands &self, const Integrands::Coninterior_penalty_parameter &v ) { self.coninterior_penalty_parameter() = v; } );
      cls.def_property( "switch_interior_penalty_method", [] ( Integrands &self ) -> Integrands::Conswitch_interior_penalty_method { return self.conswitch_interior_penalty_method(); }, [] ( Integrands &self, const Integrands::Conswitch_interior_penalty_method &v ) { self.conswitch_interior_penalty_method() = v; } );
      cls.def_property( "grad_div_stab_param", [] ( Integrands &self ) -> Integrands::Congrad_div_stab_param { return self.congrad_div_stab_param(); }, [] ( Integrands &self, const Integrands::Congrad_div_stab_param &v ) { self.congrad_div_stab_param() = v; } );
      cls.def_property( "mass_flux_stab_param", [] ( Integrands &self ) -> Integrands::Conmass_flux_stab_param { return self.conmass_flux_stab_param(); }, [] ( Integrands &self, const Integrands::Conmass_flux_stab_param &v ) { self.conmass_flux_stab_param() = v; } );
      cls.def_property( "stab_cburn", [] ( Integrands &self ) -> Integrands::Constab_cburn { return self.constab_cburn(); }, [] ( Integrands &self, const Integrands::Constab_cburn &v ) { self.constab_cburn() = v; } );
      cls.def_property( "all_boundaries", [] ( Integrands &self ) -> Integrands::Conall_boundaries { return self.conall_boundaries(); }, [] ( Integrands &self, const Integrands::Conall_boundaries &v ) { self.conall_boundaries() = v; } );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return false;});
  }
}
#endif // GuardIntegrands_0cbbbf1de3530ef80b438cb75601a15av1_3
