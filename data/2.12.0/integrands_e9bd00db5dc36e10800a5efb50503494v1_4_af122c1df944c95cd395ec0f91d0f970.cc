#ifndef GuardIntegrands_e9bd00db5dc36e10800a5efb50503494v1_4
#define GuardIntegrands_e9bd00db5dc36e10800a5efb50503494v1_4
#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fempy/py/gridview.hh>
#include <dune/python/alugrid/hierarchical.hh>
#include <dune/python/grid/gridview.hh>
#include <dune/python/grid/hierarchical.hh>
#include <cmath>
#include <tuple>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/extensions.h>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/geometry/edgelength.hh>
#include <dune/fempy/geometry/celldiameter.hh>
#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/function/simplegridfunction.hh>
#include <dune/fem/misc/gridfunctionview.hh>
#include <dune/fempy/py/integrands.hh>

namespace Integrands_e9bd00db5dc36e10800a5efb50503494v1_4
{




  // Integrands
  // ----------

  template< class GridPart, class Coeffcopy_of_solution >
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
    typedef std::tuple< std::shared_ptr< Condt > > ConstantTupleType;
    template< std::size_t i >
    using ConstantsRangeType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;
    typedef std::tuple< Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 2 > > > CoefficientFunctionSpaceTupleType;
    typedef std::tuple< Coeffcopy_of_solution > CoefficientTupleType;
    template< std::size_t i >
    using CoefficientFunctionSpaceType = std::tuple_element_t< i, CoefficientFunctionSpaceTupleType >;
    template< std::size_t i >
    using CoefficientRangeType = typename CoefficientFunctionSpaceType< i >::RangeType;
    template< std::size_t i >
    using CoefficientJacobianRangeType = typename CoefficientFunctionSpaceType< i >::JacobianRangeType;
    static constexpr bool gridPartValid = Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffcopy_of_solution>>();
    template< std::size_t i >
    using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
    template< std::size_t i >
    using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;

    Integrands ( const Coeffcopy_of_solution &coeffcopy_of_solution, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
      : coefficients_( Dune::Fem::ConstLocalFunction< Coeffcopy_of_solution >( coeffcopy_of_solution ) )
    {
      std::get< 0 >( constants_ ) = std::make_shared< Condt >( (Condt(0)) );
    }

    void decouple ()
    {
      std::get< 0 >( constants_ ) = std::make_shared< Condt >( (constant<0>()) );
    }

    bool init ( const EntityType &entity )
    {
      entity_ = entity;
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
    typedef std::tuple< Dune::FieldVector< double, 2 >, Dune::FieldMatrix< double, 2, 2 > > DomainValueType;
    typedef std::tuple< Dune::FieldVector< double, 2 >, Dune::FieldMatrix< double, 2, 2 > > RangeValueType;
    static constexpr bool _nonlinear = true;

    bool nonlinear () const
    {
      return _nonlinear;
    }

    template< class Point >
    RangeValueType interior ( const Point &x, const DomainValueType &u ) const
    {
      using std::atan;
      using std::atan2;
      using std::tan;
      typename CoefficientFunctionSpaceType< 0 >::RangeType tmp0 = evaluateCoefficient< 0 >( x );
        assert(tmp0 == tmp0);
      const auto tmp1 = -1 * tmp0[ 0 ];
        assert(tmp1 == tmp1);
      double tmp2 = constant< 0 >();
        assert(tmp2 == tmp2);
      const auto tmp3 = tmp2 / 0.0003;
        assert(tmp3 == tmp3);
      const auto tmp4 = std::get< 0 >( u );
        assert(tmp4 == tmp4);
      const auto tmp5 = tmp4[ 0 ] * tmp3;
        assert(tmp5 == tmp5);
      const auto tmp6 = -1 * tmp4[ 0 ];
        assert(tmp6 == tmp6);
      const auto tmp7 = 1.0 + tmp6;
        assert(tmp7 == tmp7);
      const auto tmp8 = tmp7 * tmp5;
        assert(tmp8 == tmp8);
      const auto tmp9 = 20.0 * tmp4[ 1 ];
        assert(tmp9 == tmp9);
      const auto tmp10 = std::atan( tmp9 );
        assert(tmp10 == tmp10);
      const auto tmp11 = 0.2864788975654116 * tmp10;
        assert(tmp11 == tmp11);
      const auto tmp12 = -1 * tmp11;
        assert(tmp12 == tmp12);
      const auto tmp13 = -0.5 + tmp4[ 0 ];
        assert(tmp13 == tmp13);
      const auto tmp14 = tmp13 + tmp12;
        assert(tmp14 == tmp14);
      const auto tmp15 = tmp14 * tmp8;
        assert(tmp15 == tmp15);
      const auto tmp16 = -1 * tmp15;
        assert(tmp16 == tmp16);
      const auto tmp17 = tmp4[ 0 ] + tmp16;
        assert(tmp17 == tmp17);
      const auto tmp18 = tmp17 + tmp1;
        assert(tmp18 == tmp18);
      const auto tmp19 = tmp0[ 1 ] + tmp1;
        assert(tmp19 == tmp19);
      const auto tmp20 = -1 * tmp19;
        assert(tmp20 == tmp20);
      const auto tmp21 = tmp4[ 1 ] + tmp6;
        assert(tmp21 == tmp21);
      const auto tmp22 = tmp21 + tmp20;
        assert(tmp22 == tmp22);
      const auto tmp23 = 0.000225 * tmp2;
        assert(tmp23 == tmp23);
      const auto tmp24 = tmp23 / 0.0003;
        assert(tmp24 == tmp24);
      typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp25 = jacobianCoefficient< 0 >( x );
        assert(tmp25 == tmp25);
      const auto tmp26 = std::atan2( (tmp25[ 0 ])[ 1 ], (tmp25[ 0 ])[ 0 ] );
        assert(tmp26 == tmp26);
      const auto tmp27 = 0.39269908169872414 + tmp26;
        assert(tmp27 == tmp27);
      const auto tmp28 = 3.0 * tmp27;
        assert(tmp28 == tmp28);
      const auto tmp29 = std::tan( tmp28 );
        assert(tmp29 == tmp29);
      const auto tmp30 = tmp29 * tmp29;
        assert(tmp30 == tmp30);
      const auto tmp31 = 1.0 + tmp30;
        assert(tmp31 == tmp31);
      const auto tmp32 = -12.0 * tmp29;
        assert(tmp32 == tmp32);
      const auto tmp33 = tmp32 / tmp31;
        assert(tmp33 == tmp33);
      const auto tmp34 = -1 * tmp30;
        assert(tmp34 == tmp34);
      const auto tmp35 = 1.0 + tmp34;
        assert(tmp35 == tmp35);
      const auto tmp36 = tmp35 / tmp31;
        assert(tmp36 == tmp36);
      const auto tmp37 = 0.02 * tmp36;
        assert(tmp37 == tmp37);
      const auto tmp38 = 1.0 + tmp37;
        assert(tmp38 == tmp38);
      const auto tmp39 = -1 * tmp38;
        assert(tmp39 == tmp39);
      const auto tmp40 = 0.02 * tmp39;
        assert(tmp40 == tmp40);
      const auto tmp41 = tmp40 * tmp33;
        assert(tmp41 == tmp41);
      const auto tmp42 = std::get< 1 >( u );
        assert(tmp42 == tmp42);
      const auto tmp43 = (tmp42[ 0 ])[ 1 ] * tmp41;
        assert(tmp43 == tmp43);
      const auto tmp44 = tmp38 * tmp38;
        assert(tmp44 == tmp44);
      const auto tmp45 = (tmp42[ 0 ])[ 0 ] * tmp44;
        assert(tmp45 == tmp45);
      const auto tmp46 = tmp45 + tmp43;
        assert(tmp46 == tmp46);
      const auto tmp47 = tmp46 * tmp24;
        assert(tmp47 == tmp47);
      const auto tmp48 = -1 * tmp41;
        assert(tmp48 == tmp48);
      const auto tmp49 = (tmp42[ 0 ])[ 0 ] * tmp48;
        assert(tmp49 == tmp49);
      const auto tmp50 = (tmp42[ 0 ])[ 1 ] * tmp44;
        assert(tmp50 == tmp50);
      const auto tmp51 = tmp50 + tmp49;
        assert(tmp51 == tmp51);
      const auto tmp52 = tmp51 * tmp24;
        assert(tmp52 == tmp52);
      const auto tmp53 = 2.25 * tmp2;
        assert(tmp53 == tmp53);
      const auto tmp54 = (tmp42[ 1 ])[ 0 ] * tmp53;
        assert(tmp54 == tmp54);
      const auto tmp55 = (tmp42[ 1 ])[ 1 ] * tmp53;
        assert(tmp55 == tmp55);
      return RangeValueType{ { tmp18, tmp22 }, { { tmp47, tmp52 }, { tmp54, tmp55 } } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      using std::atan;
      using std::atan2;
      using std::pow;
      using std::tan;
      double tmp0 = constant< 0 >();
        assert(tmp0 == tmp0);
      const auto tmp1 = tmp0 / 0.0003;
        assert(tmp1 == tmp1);
      const auto tmp2 = std::get< 0 >( u );
        assert(tmp2 == tmp2);
      const auto tmp3 = tmp2[ 0 ] * tmp1;
        assert(tmp3 == tmp3);
      const auto tmp4 = -1 * tmp2[ 0 ];
        assert(tmp4 == tmp4);
      const auto tmp5 = 1.0 + tmp4;
        assert(tmp5 == tmp5);
      const auto tmp6 = tmp5 * tmp3;
        assert(tmp6 == tmp6);
      const auto tmp7 = tmp5 * tmp1;
        assert(tmp7 == tmp7);
      const auto tmp8 = -1 * tmp3;
        assert(tmp8 == tmp8);
      const auto tmp9 = tmp8 + tmp7;
        assert(tmp9 == tmp9);
      const auto tmp10 = 20.0 * tmp2[ 1 ];
        assert(tmp10 == tmp10);
      const auto tmp11 = std::atan( tmp10 );
        assert(tmp11 == tmp11);
      const auto tmp12 = 0.2864788975654116 * tmp11;
        assert(tmp12 == tmp12);
      const auto tmp13 = -1 * tmp12;
        assert(tmp13 == tmp13);
      const auto tmp14 = -0.5 + tmp2[ 0 ];
        assert(tmp14 == tmp14);
      const auto tmp15 = tmp14 + tmp13;
        assert(tmp15 == tmp15);
      const auto tmp16 = tmp15 * tmp9;
        assert(tmp16 == tmp16);
      const auto tmp17 = tmp16 + tmp6;
        assert(tmp17 == tmp17);
      const auto tmp18 = -1 * tmp17;
        assert(tmp18 == tmp18);
      const auto tmp19 = 1 + tmp18;
        assert(tmp19 == tmp19);
      const auto tmp20 = std::pow( tmp10, 2 );
        assert(tmp20 == tmp20);
      const auto tmp21 = 1.0 + tmp20;
        assert(tmp21 == tmp21);
      const auto tmp22 = 20.0 / tmp21;
        assert(tmp22 == tmp22);
      const auto tmp23 = 0.2864788975654116 * tmp22;
        assert(tmp23 == tmp23);
      const auto tmp24 = -1 * tmp23;
        assert(tmp24 == tmp24);
      const auto tmp25 = tmp6 * tmp24;
        assert(tmp25 == tmp25);
      const auto tmp26 = -1 * tmp25;
        assert(tmp26 == tmp26);
      const auto tmp27 = 0.000225 * tmp0;
        assert(tmp27 == tmp27);
      const auto tmp28 = tmp27 / 0.0003;
        assert(tmp28 == tmp28);
      typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp29 = jacobianCoefficient< 0 >( x );
        assert(tmp29 == tmp29);
      const auto tmp30 = std::atan2( (tmp29[ 0 ])[ 1 ], (tmp29[ 0 ])[ 0 ] );
        assert(tmp30 == tmp30);
      const auto tmp31 = 0.39269908169872414 + tmp30;
        assert(tmp31 == tmp31);
      const auto tmp32 = 3.0 * tmp31;
        assert(tmp32 == tmp32);
      const auto tmp33 = std::tan( tmp32 );
        assert(tmp33 == tmp33);
      const auto tmp34 = tmp33 * tmp33;
        assert(tmp34 == tmp34);
      const auto tmp35 = 1.0 + tmp34;
        assert(tmp35 == tmp35);
      const auto tmp36 = -1 * tmp34;
        assert(tmp36 == tmp36);
      const auto tmp37 = 1.0 + tmp36;
        assert(tmp37 == tmp37);
      const auto tmp38 = tmp37 / tmp35;
        assert(tmp38 == tmp38);
      const auto tmp39 = 0.02 * tmp38;
        assert(tmp39 == tmp39);
      const auto tmp40 = 1.0 + tmp39;
        assert(tmp40 == tmp40);
      const auto tmp41 = tmp40 * tmp40;
        assert(tmp41 == tmp41);
      const auto tmp42 = tmp41 * tmp28;
        assert(tmp42 == tmp42);
      const auto tmp43 = -12.0 * tmp33;
        assert(tmp43 == tmp43);
      const auto tmp44 = tmp43 / tmp35;
        assert(tmp44 == tmp44);
      const auto tmp45 = -1 * tmp40;
        assert(tmp45 == tmp45);
      const auto tmp46 = 0.02 * tmp45;
        assert(tmp46 == tmp46);
      const auto tmp47 = tmp46 * tmp44;
        assert(tmp47 == tmp47);
      const auto tmp48 = tmp47 * tmp28;
        assert(tmp48 == tmp48);
      const auto tmp49 = -1 * tmp47;
        assert(tmp49 == tmp49);
      const auto tmp50 = tmp49 * tmp28;
        assert(tmp50 == tmp50);
      const auto tmp51 = 2.25 * tmp0;
        assert(tmp51 == tmp51);
      return [ tmp19, tmp26, tmp42, tmp48, tmp50, tmp51 ] ( const DomainValueType &phi ) {
          return RangeValueType{ { tmp19 * (std::get< 0 >( phi ))[ 0 ] + tmp26 * (std::get< 0 >( phi ))[ 1 ], -1 * (std::get< 0 >( phi ))[ 0 ] + (std::get< 0 >( phi ))[ 1 ] }, { { tmp42 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ] + tmp48 * ((std::get< 1 >( phi ))[ 0 ])[ 1 ], tmp50 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ] + tmp42 * ((std::get< 1 >( phi ))[ 0 ])[ 1 ] }, { tmp51 * ((std::get< 1 >( phi ))[ 1 ])[ 0 ], tmp51 * ((std::get< 1 >( phi ))[ 1 ])[ 1 ] } } };
        };
    }
    typedef Dune::FieldVector< double, 2 >  RRangeType;
    typedef Dune::FieldMatrix< double, 2, GridPartType::dimension >  RJacobianRangeType;
    typedef Dune::Fem::BoundaryIdProvider< typename GridPartType::GridType > BoundaryIdProviderType;
    typedef std::array<int,2> DirichletComponentType;

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

    const EntityType &entity () const
    {
      return entity_;
    }

    template< std::size_t i, class Point >
    typename CoefficientFunctionSpaceType< i >::RangeType evaluateCoefficient ( const Point &x ) const
    {
      typename CoefficientFunctionSpaceType< i >::RangeType result;
      std::get< i >( coefficients_ ).evaluate( x, result );
      return result;
    }

    template< std::size_t i, class Point >
    typename CoefficientFunctionSpaceType< i >::JacobianRangeType jacobianCoefficient ( const Point &x ) const
    {
      typename CoefficientFunctionSpaceType< i >::JacobianRangeType result;
      std::get< i >( coefficients_ ).jacobian( x, result );
      return result;
    }

    template< std::size_t i, class Point >
    typename CoefficientFunctionSpaceType< i >::HessianRangeType hessianCoefficient ( const Point &x ) const
    {
      typename CoefficientFunctionSpaceType< i >::HessianRangeType result;
      std::get< i >( coefficients_ ).hessian( x, result );
      return result;
    }
    EntityType entity_;
    IntersectionType intersection_;
    ConstantTupleType constants_;
    std::tuple< Dune::Fem::ConstLocalFunction< Coeffcopy_of_solution > > coefficients_;
  };

} // namespace Integrands_e9bd00db5dc36e10800a5efb50503494v1_4

PYBIND11_MODULE( integrands_e9bd00db5dc36e10800a5efb50503494v1_4_af122c1df944c95cd395ec0f91d0f970, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > > > GridPart;
  typedef Integrands_e9bd00db5dc36e10800a5efb50503494v1_4::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > >>, Dune::FieldVector< double, 2 > > > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_e9bd00db5dc36e10800a5efb50503494v1_4::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > >>, Dune::FieldVector< double, 2 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_e9bd00db5dc36e10800a5efb50503494v1_4_af122c1df944c95cd395ec0f91d0f970.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] ( const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > >>, Dune::FieldVector< double, 2 > > &coeffcopy_of_solution ) { return new Integrands( coeffcopy_of_solution ); } ), pybind11::keep_alive< 1, 2 >() );
      cls.def_property( "dt", [] ( Integrands &self ) -> Integrands::Condt { return self.condt(); }, [] ( Integrands &self, const Integrands::Condt &v ) { self.condt() = v; } );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return false;});
  }
}
#endif // GuardIntegrands_e9bd00db5dc36e10800a5efb50503494v1_4
