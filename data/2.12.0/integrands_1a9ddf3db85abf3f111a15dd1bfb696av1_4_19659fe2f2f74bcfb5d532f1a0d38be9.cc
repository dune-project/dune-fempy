#ifndef GuardIntegrands_1a9ddf3db85abf3f111a15dd1bfb696av1_4
#define GuardIntegrands_1a9ddf3db85abf3f111a15dd1bfb696av1_4
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
#include <dune/fempy/geometry/celldiameter.hh>
#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/function/simplegridfunction.hh>
#include <dune/fem/misc/gridfunctionview.hh>
#include <dune/fempy/py/integrands.hh>

namespace Integrands_1a9ddf3db85abf3f111a15dd1bfb696av1_4
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

    void decouple ()
    {
      std::get< 0 >( constants_ ) = std::make_shared< Coninterior_penalty_parameter >( (constant<0>()) );
      std::get< 1 >( constants_ ) = std::make_shared< Conswitch_interior_penalty_method >( (constant<1>()) );
      std::get< 2 >( constants_ ) = std::make_shared< Congrad_div_stab_param >( (constant<2>()) );
      std::get< 3 >( constants_ ) = std::make_shared< Conmass_flux_stab_param >( (constant<3>()) );
      std::get< 4 >( constants_ ) = std::make_shared< Constab_cburn >( (constant<4>()) );
      std::get< 5 >( constants_ ) = std::make_shared< Conall_boundaries >( (constant<5>()) );
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
        assert(tmp0 == tmp0);
      const auto tmp1 = -1 * tmp0[ 0 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = -1 * tmp0[ 1 ];
        assert(tmp2 == tmp2);
      const auto tmp3 = std::get< 1 >( u );
        assert(tmp3 == tmp3);
      const auto tmp4 = (tmp3[ 0 ])[ 0 ] + (tmp3[ 1 ])[ 1 ];
        assert(tmp4 == tmp4);
      const auto tmp5 = -1 * tmp4;
        assert(tmp5 == tmp5);
      double tmp6 = constant< 2 >();
        assert(tmp6 == tmp6);
      const auto tmp7 = tmp4 * tmp6;
        assert(tmp7 == tmp7);
      const auto tmp8 = std::get< 0 >( u );
        assert(tmp8 == tmp8);
      const auto tmp9 = -1 * tmp8[ 2 ];
        assert(tmp9 == tmp9);
      const auto tmp10 = (tmp3[ 0 ])[ 0 ] + tmp9;
        assert(tmp10 == tmp10);
      const auto tmp11 = tmp10 + tmp7;
        assert(tmp11 == tmp11);
      const auto tmp12 = (tmp3[ 1 ])[ 1 ] + tmp9;
        assert(tmp12 == tmp12);
      const auto tmp13 = tmp12 + tmp7;
        assert(tmp13 == tmp13);
      return RangeValueType{ { tmp1, tmp2, tmp5 }, { { tmp11, (tmp3[ 0 ])[ 1 ] }, { (tmp3[ 1 ])[ 0 ], tmp13 }, { 0.0, 0.0 } } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      double tmp0 = constant< 2 >();
        assert(tmp0 == tmp0);
      const auto tmp1 = 1 + tmp0;
        assert(tmp1 == tmp1);
      return [ tmp0, tmp1 ] ( const DomainValueType &phi ) {
          return RangeValueType{ { 0, 0, -1 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ] + -1 * ((std::get< 1 >( phi ))[ 1 ])[ 1 ] }, { { (-1 * (std::get< 0 >( phi ))[ 2 ] + tmp1 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ]) + tmp0 * ((std::get< 1 >( phi ))[ 1 ])[ 1 ], ((std::get< 1 >( phi ))[ 0 ])[ 1 ] }, { ((std::get< 1 >( phi ))[ 1 ])[ 0 ], (-1 * (std::get< 0 >( phi ))[ 2 ] + tmp0 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ]) + tmp1 * ((std::get< 1 >( phi ))[ 1 ])[ 1 ] }, { 0, 0 } } };
        };
    }

    template< class Point >
    RangeValueType boundary ( const Point &x, const DomainValueType &u ) const
    {
      double tmp0 = constant< 3 >();
        assert(tmp0 == tmp0);
      auto tmp1 = facetArea_;
        assert(tmp1 == tmp1);
      auto tmp2 = cellVolume();
        assert(tmp2 == tmp2);
      const auto tmp3 = tmp2 / tmp1;
        assert(tmp3 == tmp3);
      const auto tmp4 = 1 / tmp3;
        assert(tmp4 == tmp4);
      typename CoefficientFunctionSpaceType< 2 >::RangeType tmp5 = evaluateCoefficient< 2 >( x );
        assert(tmp5 == tmp5);
      const auto tmp6 = -1 * tmp5[ 0 ];
        assert(tmp6 == tmp6);
      const auto tmp7 = std::get< 0 >( u );
        assert(tmp7 == tmp7);
      const auto tmp8 = tmp7[ 1 ] + tmp6;
        assert(tmp8 == tmp8);
      GlobalCoordinateType tmp9 = intersection_.unitOuterNormal( x.localPosition() );
        assert(tmp9 == tmp9);
      const auto tmp10 = tmp9[ 1 ] * tmp8;
        assert(tmp10 == tmp10);
      typename CoefficientFunctionSpaceType< 1 >::RangeType tmp11 = evaluateCoefficient< 1 >( x );
        assert(tmp11 == tmp11);
      const auto tmp12 = -1 * tmp11[ 0 ];
        assert(tmp12 == tmp12);
      const auto tmp13 = tmp7[ 0 ] + tmp12;
        assert(tmp13 == tmp13);
      const auto tmp14 = tmp9[ 0 ] * tmp13;
        assert(tmp14 == tmp14);
      const auto tmp15 = tmp14 + tmp10;
        assert(tmp15 == tmp15);
      const auto tmp16 = tmp15 * tmp4;
        assert(tmp16 == tmp16);
      const auto tmp17 = tmp9[ 0 ] * tmp16;
        assert(tmp17 == tmp17);
      const auto tmp18 = tmp17 * tmp0;
        assert(tmp18 == tmp18);
      const auto tmp19 = tmp9[ 0 ] * tmp7[ 2 ];
        assert(tmp19 == tmp19);
      double tmp20 = constant< 5 >();
        assert(tmp20 == tmp20);
      double tmp21 = constant< 0 >();
        assert(tmp21 == tmp21);
      const auto tmp22 = tmp21 / tmp3;
        assert(tmp22 == tmp22);
      const auto tmp23 = tmp13 * tmp22;
        assert(tmp23 == tmp23);
      const auto tmp24 = tmp23 * tmp20;
        assert(tmp24 == tmp24);
      const auto tmp25 = std::get< 1 >( u );
        assert(tmp25 == tmp25);
      const auto tmp26 = (tmp25[ 0 ])[ 1 ] * tmp9[ 1 ];
        assert(tmp26 == tmp26);
      const auto tmp27 = tmp9[ 0 ] * (tmp25[ 0 ])[ 0 ];
        assert(tmp27 == tmp27);
      const auto tmp28 = tmp27 + tmp26;
        assert(tmp28 == tmp28);
      const auto tmp29 = tmp28 * tmp20;
        assert(tmp29 == tmp29);
      const auto tmp30 = -1 * tmp29;
        assert(tmp30 == tmp30);
      const auto tmp31 = tmp30 + tmp24;
        assert(tmp31 == tmp31);
      const auto tmp32 = tmp31 + tmp19;
        assert(tmp32 == tmp32);
      const auto tmp33 = tmp32 + tmp18;
        assert(tmp33 == tmp33);
      const auto tmp34 = tmp9[ 1 ] * tmp16;
        assert(tmp34 == tmp34);
      const auto tmp35 = tmp34 * tmp0;
        assert(tmp35 == tmp35);
      const auto tmp36 = tmp9[ 1 ] * tmp7[ 2 ];
        assert(tmp36 == tmp36);
      const auto tmp37 = tmp8 * tmp22;
        assert(tmp37 == tmp37);
      const auto tmp38 = tmp37 * tmp20;
        assert(tmp38 == tmp38);
      const auto tmp39 = tmp9[ 1 ] * (tmp25[ 1 ])[ 1 ];
        assert(tmp39 == tmp39);
      const auto tmp40 = tmp9[ 0 ] * (tmp25[ 1 ])[ 0 ];
        assert(tmp40 == tmp40);
      const auto tmp41 = tmp40 + tmp39;
        assert(tmp41 == tmp41);
      const auto tmp42 = tmp41 * tmp20;
        assert(tmp42 == tmp42);
      const auto tmp43 = -1 * tmp42;
        assert(tmp43 == tmp43);
      const auto tmp44 = tmp43 + tmp38;
        assert(tmp44 == tmp44);
      const auto tmp45 = tmp44 + tmp36;
        assert(tmp45 == tmp45);
      const auto tmp46 = tmp45 + tmp35;
        assert(tmp46 == tmp46);
      double tmp47 = constant< 1 >();
        assert(tmp47 == tmp47);
      const auto tmp48 = tmp14 * tmp47;
        assert(tmp48 == tmp48);
      const auto tmp49 = tmp48 * tmp20;
        assert(tmp49 == tmp49);
      const auto tmp50 = tmp9[ 1 ] * tmp13;
        assert(tmp50 == tmp50);
      const auto tmp51 = tmp50 * tmp47;
        assert(tmp51 == tmp51);
      const auto tmp52 = tmp51 * tmp20;
        assert(tmp52 == tmp52);
      const auto tmp53 = tmp9[ 0 ] * tmp8;
        assert(tmp53 == tmp53);
      const auto tmp54 = tmp53 * tmp47;
        assert(tmp54 == tmp54);
      const auto tmp55 = tmp54 * tmp20;
        assert(tmp55 == tmp55);
      const auto tmp56 = tmp10 * tmp47;
        assert(tmp56 == tmp56);
      const auto tmp57 = tmp56 * tmp20;
        assert(tmp57 == tmp57);
      return RangeValueType{ { tmp33, tmp46, tmp15 }, { { tmp49, tmp52 }, { tmp55, tmp57 }, { 0.0, 0.0 } } };
    }

    template< class Point >
    auto linearizedBoundary ( const Point &x, const DomainValueType &u ) const
    {
      double tmp0 = constant< 5 >();
        assert(tmp0 == tmp0);
      auto tmp1 = facetArea_;
        assert(tmp1 == tmp1);
      auto tmp2 = cellVolume();
        assert(tmp2 == tmp2);
      const auto tmp3 = tmp2 / tmp1;
        assert(tmp3 == tmp3);
      double tmp4 = constant< 0 >();
        assert(tmp4 == tmp4);
      const auto tmp5 = tmp4 / tmp3;
        assert(tmp5 == tmp5);
      const auto tmp6 = tmp5 * tmp0;
        assert(tmp6 == tmp6);
      double tmp7 = constant< 3 >();
        assert(tmp7 == tmp7);
      const auto tmp8 = 1 / tmp3;
        assert(tmp8 == tmp8);
      GlobalCoordinateType tmp9 = intersection_.unitOuterNormal( x.localPosition() );
        assert(tmp9 == tmp9);
      const auto tmp10 = tmp9[ 0 ] * tmp8;
        assert(tmp10 == tmp10);
      const auto tmp11 = tmp9[ 0 ] * tmp10;
        assert(tmp11 == tmp11);
      const auto tmp12 = tmp11 * tmp7;
        assert(tmp12 == tmp12);
      const auto tmp13 = tmp12 + tmp6;
        assert(tmp13 == tmp13);
      const auto tmp14 = tmp9[ 1 ] * tmp8;
        assert(tmp14 == tmp14);
      const auto tmp15 = tmp9[ 0 ] * tmp14;
        assert(tmp15 == tmp15);
      const auto tmp16 = tmp15 * tmp7;
        assert(tmp16 == tmp16);
      const auto tmp17 = tmp9[ 1 ] * tmp10;
        assert(tmp17 == tmp17);
      const auto tmp18 = tmp17 * tmp7;
        assert(tmp18 == tmp18);
      const auto tmp19 = tmp9[ 1 ] * tmp14;
        assert(tmp19 == tmp19);
      const auto tmp20 = tmp19 * tmp7;
        assert(tmp20 == tmp20);
      const auto tmp21 = tmp20 + tmp6;
        assert(tmp21 == tmp21);
      const auto tmp22 = tmp9[ 0 ] * tmp0;
        assert(tmp22 == tmp22);
      const auto tmp23 = -1 * tmp22;
        assert(tmp23 == tmp23);
      const auto tmp24 = tmp9[ 1 ] * tmp0;
        assert(tmp24 == tmp24);
      const auto tmp25 = -1 * tmp24;
        assert(tmp25 == tmp25);
      double tmp26 = constant< 1 >();
        assert(tmp26 == tmp26);
      const auto tmp27 = tmp9[ 0 ] * tmp26;
        assert(tmp27 == tmp27);
      const auto tmp28 = tmp27 * tmp0;
        assert(tmp28 == tmp28);
      const auto tmp29 = tmp9[ 1 ] * tmp26;
        assert(tmp29 == tmp29);
      const auto tmp30 = tmp29 * tmp0;
        assert(tmp30 == tmp30);
      return [ tmp13, tmp16, tmp18, tmp21, tmp23, tmp25, tmp28, tmp30, tmp9 ] ( const DomainValueType &phi ) {
          return RangeValueType{ { (((tmp13 * (std::get< 0 >( phi ))[ 0 ] + tmp16 * (std::get< 0 >( phi ))[ 1 ]) + tmp9[ 0 ] * (std::get< 0 >( phi ))[ 2 ]) + tmp23 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ]) + tmp25 * ((std::get< 1 >( phi ))[ 0 ])[ 1 ], (((tmp18 * (std::get< 0 >( phi ))[ 0 ] + tmp21 * (std::get< 0 >( phi ))[ 1 ]) + tmp9[ 1 ] * (std::get< 0 >( phi ))[ 2 ]) + tmp23 * ((std::get< 1 >( phi ))[ 1 ])[ 0 ]) + tmp25 * ((std::get< 1 >( phi ))[ 1 ])[ 1 ], tmp9[ 0 ] * (std::get< 0 >( phi ))[ 0 ] + tmp9[ 1 ] * (std::get< 0 >( phi ))[ 1 ] }, { { tmp28 * (std::get< 0 >( phi ))[ 0 ], tmp30 * (std::get< 0 >( phi ))[ 0 ] }, { tmp28 * (std::get< 0 >( phi ))[ 1 ], tmp30 * (std::get< 0 >( phi ))[ 1 ] }, { 0, 0 } } };
        };
    }

    template< class Point >
    std::pair< RangeValueType, RangeValueType > skeleton ( const Point &xIn, const DomainValueType &uIn, const Point &xOut, const DomainValueType &uOut ) const
    {
      double tmp0 = constant< 3 >();
        assert(tmp0 == tmp0);
      auto tmp1 = facetArea_;
        assert(tmp1 == tmp1);
      auto tmp2 = cellVolume_[ static_cast< std::size_t >( Side::out ) ];
        assert(tmp2 == tmp2);
      auto tmp3 = cellVolume_[ static_cast< std::size_t >( Side::in ) ];
        assert(tmp3 == tmp3);
      const auto tmp4 = tmp3 + tmp2;
        assert(tmp4 == tmp4);
      const auto tmp5 = 0.5 * tmp4;
        assert(tmp5 == tmp5);
      const auto tmp6 = tmp5 / tmp1;
        assert(tmp6 == tmp6);
      const auto tmp7 = 1 / tmp6;
        assert(tmp7 == tmp7);
      const auto tmp8 = std::get< 0 >( uOut );
        assert(tmp8 == tmp8);
      const auto tmp9 = -1 * tmp8[ 1 ];
        assert(tmp9 == tmp9);
      const auto tmp10 = std::get< 0 >( uIn );
        assert(tmp10 == tmp10);
      const auto tmp11 = tmp10[ 1 ] + tmp9;
        assert(tmp11 == tmp11);
      GlobalCoordinateType tmp12 = intersection_.unitOuterNormal( xIn.localPosition() );
        assert(tmp12 == tmp12);
      const auto tmp13 = tmp12[ 1 ] * tmp11;
        assert(tmp13 == tmp13);
      const auto tmp14 = -1 * tmp8[ 0 ];
        assert(tmp14 == tmp14);
      const auto tmp15 = tmp10[ 0 ] + tmp14;
        assert(tmp15 == tmp15);
      const auto tmp16 = tmp12[ 0 ] * tmp15;
        assert(tmp16 == tmp16);
      const auto tmp17 = tmp16 + tmp13;
        assert(tmp17 == tmp17);
      const auto tmp18 = tmp17 * tmp7;
        assert(tmp18 == tmp18);
      const auto tmp19 = tmp12[ 0 ] * tmp18;
        assert(tmp19 == tmp19);
      const auto tmp20 = tmp19 * tmp0;
        assert(tmp20 == tmp20);
      const auto tmp21 = tmp10[ 2 ] + tmp8[ 2 ];
        assert(tmp21 == tmp21);
      const auto tmp22 = 0.5 * tmp21;
        assert(tmp22 == tmp22);
      const auto tmp23 = tmp12[ 0 ] * tmp22;
        assert(tmp23 == tmp23);
      double tmp24 = constant< 0 >();
        assert(tmp24 == tmp24);
      const auto tmp25 = tmp24 / tmp6;
        assert(tmp25 == tmp25);
      const auto tmp26 = tmp15 * tmp25;
        assert(tmp26 == tmp26);
      const auto tmp27 = std::get< 1 >( uOut );
        assert(tmp27 == tmp27);
      const auto tmp28 = std::get< 1 >( uIn );
        assert(tmp28 == tmp28);
      const auto tmp29 = (tmp28[ 0 ])[ 1 ] + (tmp27[ 0 ])[ 1 ];
        assert(tmp29 == tmp29);
      const auto tmp30 = 0.5 * tmp29;
        assert(tmp30 == tmp30);
      const auto tmp31 = tmp12[ 1 ] * tmp30;
        assert(tmp31 == tmp31);
      const auto tmp32 = (tmp28[ 0 ])[ 0 ] + (tmp27[ 0 ])[ 0 ];
        assert(tmp32 == tmp32);
      const auto tmp33 = 0.5 * tmp32;
        assert(tmp33 == tmp33);
      const auto tmp34 = tmp12[ 0 ] * tmp33;
        assert(tmp34 == tmp34);
      const auto tmp35 = tmp34 + tmp31;
        assert(tmp35 == tmp35);
      const auto tmp36 = -1 * tmp35;
        assert(tmp36 == tmp36);
      const auto tmp37 = tmp36 + tmp26;
        assert(tmp37 == tmp37);
      const auto tmp38 = tmp37 + tmp23;
        assert(tmp38 == tmp38);
      const auto tmp39 = tmp38 + tmp20;
        assert(tmp39 == tmp39);
      const auto tmp40 = tmp12[ 1 ] * tmp18;
        assert(tmp40 == tmp40);
      const auto tmp41 = tmp40 * tmp0;
        assert(tmp41 == tmp41);
      const auto tmp42 = tmp12[ 1 ] * tmp22;
        assert(tmp42 == tmp42);
      const auto tmp43 = tmp11 * tmp25;
        assert(tmp43 == tmp43);
      const auto tmp44 = (tmp28[ 1 ])[ 1 ] + (tmp27[ 1 ])[ 1 ];
        assert(tmp44 == tmp44);
      const auto tmp45 = 0.5 * tmp44;
        assert(tmp45 == tmp45);
      const auto tmp46 = tmp12[ 1 ] * tmp45;
        assert(tmp46 == tmp46);
      const auto tmp47 = (tmp28[ 1 ])[ 0 ] + (tmp27[ 1 ])[ 0 ];
        assert(tmp47 == tmp47);
      const auto tmp48 = 0.5 * tmp47;
        assert(tmp48 == tmp48);
      const auto tmp49 = tmp12[ 0 ] * tmp48;
        assert(tmp49 == tmp49);
      const auto tmp50 = tmp49 + tmp46;
        assert(tmp50 == tmp50);
      const auto tmp51 = -1 * tmp50;
        assert(tmp51 == tmp51);
      const auto tmp52 = tmp51 + tmp43;
        assert(tmp52 == tmp52);
      const auto tmp53 = tmp52 + tmp42;
        assert(tmp53 == tmp53);
      const auto tmp54 = tmp53 + tmp41;
        assert(tmp54 == tmp54);
      double tmp55 = constant< 4 >();
        assert(tmp55 == tmp55);
      const auto tmp56 = tmp6 * tmp55;
        assert(tmp56 == tmp56);
      const auto tmp57 = -1 * tmp8[ 2 ];
        assert(tmp57 == tmp57);
      const auto tmp58 = tmp10[ 2 ] + tmp57;
        assert(tmp58 == tmp58);
      const auto tmp59 = tmp12[ 1 ] * tmp58;
        assert(tmp59 == tmp59);
      const auto tmp60 = tmp12[ 1 ] * tmp59;
        assert(tmp60 == tmp60);
      const auto tmp61 = tmp12[ 0 ] * tmp58;
        assert(tmp61 == tmp61);
      const auto tmp62 = tmp12[ 0 ] * tmp61;
        assert(tmp62 == tmp62);
      const auto tmp63 = tmp62 + tmp60;
        assert(tmp63 == tmp63);
      const auto tmp64 = tmp63 * tmp56;
        assert(tmp64 == tmp64);
      const auto tmp65 = 0.5 * tmp17;
        assert(tmp65 == tmp65);
      const auto tmp66 = tmp65 + tmp64;
        assert(tmp66 == tmp66);
      double tmp67 = constant< 1 >();
        assert(tmp67 == tmp67);
      const auto tmp68 = 0.5 * tmp12[ 0 ];
        assert(tmp68 == tmp68);
      const auto tmp69 = tmp15 * tmp68;
        assert(tmp69 == tmp69);
      const auto tmp70 = tmp69 * tmp67;
        assert(tmp70 == tmp70);
      const auto tmp71 = 0.5 * tmp12[ 1 ];
        assert(tmp71 == tmp71);
      const auto tmp72 = tmp15 * tmp71;
        assert(tmp72 == tmp72);
      const auto tmp73 = tmp72 * tmp67;
        assert(tmp73 == tmp73);
      const auto tmp74 = tmp11 * tmp68;
        assert(tmp74 == tmp74);
      const auto tmp75 = tmp74 * tmp67;
        assert(tmp75 == tmp75);
      const auto tmp76 = tmp11 * tmp71;
        assert(tmp76 == tmp76);
      const auto tmp77 = tmp76 * tmp67;
        assert(tmp77 == tmp77);
      const auto tmp78 = -1 * tmp12[ 0 ];
        assert(tmp78 == tmp78);
      const auto tmp79 = tmp78 * tmp18;
        assert(tmp79 == tmp79);
      const auto tmp80 = tmp79 * tmp0;
        assert(tmp80 == tmp80);
      const auto tmp81 = tmp78 * tmp22;
        assert(tmp81 == tmp81);
      const auto tmp82 = -1 * tmp15;
        assert(tmp82 == tmp82);
      const auto tmp83 = tmp82 * tmp25;
        assert(tmp83 == tmp83);
      const auto tmp84 = -1 * tmp36;
        assert(tmp84 == tmp84);
      const auto tmp85 = tmp84 + tmp83;
        assert(tmp85 == tmp85);
      const auto tmp86 = tmp85 + tmp81;
        assert(tmp86 == tmp86);
      const auto tmp87 = tmp86 + tmp80;
        assert(tmp87 == tmp87);
      const auto tmp88 = -1 * tmp12[ 1 ];
        assert(tmp88 == tmp88);
      const auto tmp89 = tmp88 * tmp18;
        assert(tmp89 == tmp89);
      const auto tmp90 = tmp89 * tmp0;
        assert(tmp90 == tmp90);
      const auto tmp91 = tmp88 * tmp22;
        assert(tmp91 == tmp91);
      const auto tmp92 = -1 * tmp11;
        assert(tmp92 == tmp92);
      const auto tmp93 = tmp92 * tmp25;
        assert(tmp93 == tmp93);
      const auto tmp94 = -1 * tmp51;
        assert(tmp94 == tmp94);
      const auto tmp95 = tmp94 + tmp93;
        assert(tmp95 == tmp95);
      const auto tmp96 = tmp95 + tmp91;
        assert(tmp96 == tmp96);
      const auto tmp97 = tmp96 + tmp90;
        assert(tmp97 == tmp97);
      const auto tmp98 = tmp88 * tmp59;
        assert(tmp98 == tmp98);
      const auto tmp99 = tmp78 * tmp61;
        assert(tmp99 == tmp99);
      const auto tmp100 = tmp99 + tmp98;
        assert(tmp100 == tmp100);
      const auto tmp101 = tmp100 * tmp56;
        assert(tmp101 == tmp101);
      const auto tmp102 = tmp65 + tmp101;
        assert(tmp102 == tmp102);
      return std::make_pair( RangeValueType{ { tmp39, tmp54, tmp66 }, { { tmp70, tmp73 }, { tmp75, tmp77 }, { 0.0, 0.0 } } }, RangeValueType{ { tmp87, tmp97, tmp102 }, { { tmp70, tmp73 }, { tmp75, tmp77 }, { 0.0, 0.0 } } } );
    }

    template< class Point >
    auto linearizedSkeleton ( const Point &xIn, const DomainValueType &uIn, const Point &xOut, const DomainValueType &uOut ) const
    {
      auto tmp0 = facetArea_;
        assert(tmp0 == tmp0);
      auto tmp1 = cellVolume_[ static_cast< std::size_t >( Side::out ) ];
        assert(tmp1 == tmp1);
      auto tmp2 = cellVolume_[ static_cast< std::size_t >( Side::in ) ];
        assert(tmp2 == tmp2);
      const auto tmp3 = tmp2 + tmp1;
        assert(tmp3 == tmp3);
      const auto tmp4 = 0.5 * tmp3;
        assert(tmp4 == tmp4);
      const auto tmp5 = tmp4 / tmp0;
        assert(tmp5 == tmp5);
      double tmp6 = constant< 0 >();
        assert(tmp6 == tmp6);
      const auto tmp7 = tmp6 / tmp5;
        assert(tmp7 == tmp7);
      double tmp8 = constant< 3 >();
        assert(tmp8 == tmp8);
      const auto tmp9 = 1 / tmp5;
        assert(tmp9 == tmp9);
      GlobalCoordinateType tmp10 = intersection_.unitOuterNormal( xIn.localPosition() );
        assert(tmp10 == tmp10);
      const auto tmp11 = tmp10[ 0 ] * tmp9;
        assert(tmp11 == tmp11);
      const auto tmp12 = tmp10[ 0 ] * tmp11;
        assert(tmp12 == tmp12);
      const auto tmp13 = tmp12 * tmp8;
        assert(tmp13 == tmp13);
      const auto tmp14 = tmp13 + tmp7;
        assert(tmp14 == tmp14);
      const auto tmp15 = tmp10[ 1 ] * tmp9;
        assert(tmp15 == tmp15);
      const auto tmp16 = tmp10[ 0 ] * tmp15;
        assert(tmp16 == tmp16);
      const auto tmp17 = tmp16 * tmp8;
        assert(tmp17 == tmp17);
      const auto tmp18 = 0.5 * tmp10[ 0 ];
        assert(tmp18 == tmp18);
      const auto tmp19 = tmp10[ 1 ] * tmp11;
        assert(tmp19 == tmp19);
      const auto tmp20 = tmp19 * tmp8;
        assert(tmp20 == tmp20);
      const auto tmp21 = tmp10[ 1 ] * tmp15;
        assert(tmp21 == tmp21);
      const auto tmp22 = tmp21 * tmp8;
        assert(tmp22 == tmp22);
      const auto tmp23 = tmp22 + tmp7;
        assert(tmp23 == tmp23);
      const auto tmp24 = 0.5 * tmp10[ 1 ];
        assert(tmp24 == tmp24);
      double tmp25 = constant< 4 >();
        assert(tmp25 == tmp25);
      const auto tmp26 = tmp5 * tmp25;
        assert(tmp26 == tmp26);
      const auto tmp27 = tmp10[ 1 ] * tmp10[ 1 ];
        assert(tmp27 == tmp27);
      const auto tmp28 = tmp10[ 0 ] * tmp10[ 0 ];
        assert(tmp28 == tmp28);
      const auto tmp29 = tmp28 + tmp27;
        assert(tmp29 == tmp29);
      const auto tmp30 = tmp29 * tmp26;
        assert(tmp30 == tmp30);
      const auto tmp31 = -1 * tmp18;
        assert(tmp31 == tmp31);
      const auto tmp32 = -1 * tmp24;
        assert(tmp32 == tmp32);
      double tmp33 = constant< 1 >();
        assert(tmp33 == tmp33);
      const auto tmp34 = tmp18 * tmp33;
        assert(tmp34 == tmp34);
      const auto tmp35 = tmp24 * tmp33;
        assert(tmp35 == tmp35);
      const auto tmp36 = -1 * tmp10[ 0 ];
        assert(tmp36 == tmp36);
      const auto tmp37 = tmp36 * tmp11;
        assert(tmp37 == tmp37);
      const auto tmp38 = tmp37 * tmp8;
        assert(tmp38 == tmp38);
      const auto tmp39 = -1 * tmp7;
        assert(tmp39 == tmp39);
      const auto tmp40 = tmp39 + tmp38;
        assert(tmp40 == tmp40);
      const auto tmp41 = tmp36 * tmp15;
        assert(tmp41 == tmp41);
      const auto tmp42 = tmp41 * tmp8;
        assert(tmp42 == tmp42);
      const auto tmp43 = 0.5 * tmp36;
        assert(tmp43 == tmp43);
      const auto tmp44 = -1 * tmp10[ 1 ];
        assert(tmp44 == tmp44);
      const auto tmp45 = tmp44 * tmp11;
        assert(tmp45 == tmp45);
      const auto tmp46 = tmp45 * tmp8;
        assert(tmp46 == tmp46);
      const auto tmp47 = tmp44 * tmp15;
        assert(tmp47 == tmp47);
      const auto tmp48 = tmp47 * tmp8;
        assert(tmp48 == tmp48);
      const auto tmp49 = tmp39 + tmp48;
        assert(tmp49 == tmp49);
      const auto tmp50 = 0.5 * tmp44;
        assert(tmp50 == tmp50);
      const auto tmp51 = tmp10[ 1 ] * tmp44;
        assert(tmp51 == tmp51);
      const auto tmp52 = tmp10[ 0 ] * tmp36;
        assert(tmp52 == tmp52);
      const auto tmp53 = tmp52 + tmp51;
        assert(tmp53 == tmp53);
      const auto tmp54 = tmp53 * tmp26;
        assert(tmp54 == tmp54);
      const auto tmp55 = -1 * tmp31;
        assert(tmp55 == tmp55);
      const auto tmp56 = -1 * tmp32;
        assert(tmp56 == tmp56);
      const auto tmp57 = tmp36 * tmp9;
        assert(tmp57 == tmp57);
      const auto tmp58 = tmp10[ 0 ] * tmp57;
        assert(tmp58 == tmp58);
      const auto tmp59 = tmp58 * tmp8;
        assert(tmp59 == tmp59);
      const auto tmp60 = tmp39 + tmp59;
        assert(tmp60 == tmp60);
      const auto tmp61 = tmp44 * tmp9;
        assert(tmp61 == tmp61);
      const auto tmp62 = tmp10[ 0 ] * tmp61;
        assert(tmp62 == tmp62);
      const auto tmp63 = tmp62 * tmp8;
        assert(tmp63 == tmp63);
      const auto tmp64 = tmp10[ 1 ] * tmp57;
        assert(tmp64 == tmp64);
      const auto tmp65 = tmp64 * tmp8;
        assert(tmp65 == tmp65);
      const auto tmp66 = tmp10[ 1 ] * tmp61;
        assert(tmp66 == tmp66);
      const auto tmp67 = tmp66 * tmp8;
        assert(tmp67 == tmp67);
      const auto tmp68 = tmp39 + tmp67;
        assert(tmp68 == tmp68);
      const auto tmp69 = tmp31 * tmp33;
        assert(tmp69 == tmp69);
      const auto tmp70 = tmp32 * tmp33;
        assert(tmp70 == tmp70);
      const auto tmp71 = tmp36 * tmp57;
        assert(tmp71 == tmp71);
      const auto tmp72 = tmp71 * tmp8;
        assert(tmp72 == tmp72);
      const auto tmp73 = tmp72 + tmp7;
        assert(tmp73 == tmp73);
      const auto tmp74 = tmp36 * tmp61;
        assert(tmp74 == tmp74);
      const auto tmp75 = tmp74 * tmp8;
        assert(tmp75 == tmp75);
      const auto tmp76 = tmp44 * tmp57;
        assert(tmp76 == tmp76);
      const auto tmp77 = tmp76 * tmp8;
        assert(tmp77 == tmp77);
      const auto tmp78 = tmp44 * tmp61;
        assert(tmp78 == tmp78);
      const auto tmp79 = tmp78 * tmp8;
        assert(tmp79 == tmp79);
      const auto tmp80 = tmp79 + tmp7;
        assert(tmp80 == tmp80);
      const auto tmp81 = tmp44 * tmp44;
        assert(tmp81 == tmp81);
      const auto tmp82 = tmp36 * tmp36;
        assert(tmp82 == tmp82);
      const auto tmp83 = tmp82 + tmp81;
        assert(tmp83 == tmp83);
      const auto tmp84 = tmp83 * tmp26;
        assert(tmp84 == tmp84);
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

} // namespace Integrands_1a9ddf3db85abf3f111a15dd1bfb696av1_4

PYBIND11_MODULE( integrands_1a9ddf3db85abf3f111a15dd1bfb696av1_4_19659fe2f2f74bcfb5d532f1a0d38be9, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView > GridPart;
  typedef Integrands_1a9ddf3db85abf3f111a15dd1bfb696av1_4::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 2 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_1a9ddf3db85abf3f111a15dd1bfb696av1_4::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 2 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_1a9ddf3db85abf3f111a15dd1bfb696av1_4_19659fe2f2f74bcfb5d532f1a0d38be9.cc"})).first;
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
#endif // GuardIntegrands_1a9ddf3db85abf3f111a15dd1bfb696av1_4
