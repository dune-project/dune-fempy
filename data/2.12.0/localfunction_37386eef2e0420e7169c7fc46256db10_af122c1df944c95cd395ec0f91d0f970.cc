#ifndef GUARD_37386eef2e0420e7169c7fc46256db10
#define GUARD_37386eef2e0420e7169c7fc46256db10

#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fempy/py/gridview.hh>
#include <dune/python/alugrid/hierarchical.hh>
#include <dune/python/grid/gridview.hh>
#include <dune/python/grid/hierarchical.hh>
#include <dune/fem/function/localfunction/bindable.hh>
#include <dune/fem/common/intersectionside.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/extensions.h>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/common/exceptions.hh>
#include <dune/fempy/py/ufllocalfunction.hh>

namespace UFLLocalFunctions_37386eef2e0420e7169c7fc46256db10
{

  // UFLLocalFunction
// ----------------

template< class GridPart >
struct UFLLocalFunction
  : public Dune::Fem::BindableGridFunctionWithSpace<GridPart,Dune::Dim<1>>
{
  typedef GridPart GridPartType;
  typedef typename GridPartType::GridViewType GridView;
  typedef typename GridView::ctype ctype;
  typedef Dune::Fem::BindableGridFunctionWithSpace<GridPart,Dune::Dim<1>> BaseType;
  typedef Dune::Fem::GridFunctionSpace<GridPartType,Dune::Dim<1>> FunctionSpaceType;
  typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
  typedef typename GridPartType::IntersectionType IntersectionType;
  typedef typename EntityType::Geometry Geometry;
  typedef typename Geometry::GlobalCoordinate GlobalCoordinateType;
  typedef Dune::Fem::IntersectionSide Side;
  typedef double Contime;
  typedef std::tuple< std::shared_ptr< Contime > > ConstantTupleType;
  template< std::size_t i >
  using ConstantsRangeType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;
  typedef std::tuple<> CoefficientTupleType;
  static constexpr bool gridPartValid = true;
  template< std::size_t i >
  using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
  template< std::size_t i >
  using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;
  using BaseType::entity;
  using BaseType::geometry;

  UFLLocalFunction ( const GridPartType &gridPart, const std::string &name, int order, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    : BaseType(gridPart,name,order)
  {
    std::get< 0 >( constants_ ) = std::make_shared< Contime >( (Contime(0)) );
  }

  void decouple ()
  {
    std::get< 0 >( constants_ ) = std::make_shared< Contime >( (constant<0>()) );
  }

  void bind ( const IntersectionType &intersection, Side side )
  {
    BaseType::bind(intersection, side);
  }

  void bind ( const EntityType &entity )
  {
    BaseType::bind(entity);
  }

  void unbind ()
  {
    BaseType::unbind();
  }

  template< class Point >
  void evaluate ( const Point &x, typename FunctionSpaceType::RangeType &result ) const
  {
    using std::tanh;
    double tmp0 = constant< 0 >();
      assert(tmp0 == tmp0);
    GlobalCoordinateType tmp1 = geometry().global( Dune::Fem::coordinate( x ) );
      assert(tmp1 == tmp1);
    const auto tmp2 = tmp1[ 0 ] * tmp0;
      assert(tmp2 == tmp2);
    const auto tmp3 = -1 * tmp1[ 1 ];
      assert(tmp3 == tmp3);
    const auto tmp4 = tmp3 + tmp2;
      assert(tmp4 == tmp4);
    const auto tmp5 = 2 * tmp4;
      assert(tmp5 == tmp5);
    const auto tmp6 = std::tanh( tmp5 );
      assert(tmp6 == tmp6);
    const auto tmp7 = tmp1[ 1 ] * tmp1[ 1 ];
      assert(tmp7 == tmp7);
    const auto tmp8 = tmp1[ 0 ] * tmp1[ 0 ];
      assert(tmp8 == tmp8);
    const auto tmp9 = tmp8 + tmp7;
      assert(tmp9 == tmp9);
    const auto tmp10 = tmp9 * tmp6;
      assert(tmp10 == tmp10);
    result[ 0 ] = tmp10;
  }

  template< class Point >
  void jacobian ( const Point &x, typename FunctionSpaceType::JacobianRangeType &result ) const
  {
    using std::cosh;
    using std::pow;
    using std::tanh;
    double tmp0 = constant< 0 >();
      assert(tmp0 == tmp0);
    GlobalCoordinateType tmp1 = geometry().global( Dune::Fem::coordinate( x ) );
      assert(tmp1 == tmp1);
    const auto tmp2 = tmp1[ 0 ] * tmp0;
      assert(tmp2 == tmp2);
    const auto tmp3 = -1 * tmp1[ 1 ];
      assert(tmp3 == tmp3);
    const auto tmp4 = tmp3 + tmp2;
      assert(tmp4 == tmp4);
    const auto tmp5 = 2 * tmp4;
      assert(tmp5 == tmp5);
    const auto tmp6 = std::tanh( tmp5 );
      assert(tmp6 == tmp6);
    const auto tmp7 = tmp1[ 0 ] + tmp1[ 0 ];
      assert(tmp7 == tmp7);
    const auto tmp8 = tmp7 * tmp6;
      assert(tmp8 == tmp8);
    const auto tmp9 = 2.0 * tmp5;
      assert(tmp9 == tmp9);
    const auto tmp10 = std::cosh( tmp9 );
      assert(tmp10 == tmp10);
    const auto tmp11 = 1.0 + tmp10;
      assert(tmp11 == tmp11);
    const auto tmp12 = std::cosh( tmp5 );
      assert(tmp12 == tmp12);
    const auto tmp13 = 2.0 * tmp12;
      assert(tmp13 == tmp13);
    const auto tmp14 = tmp13 / tmp11;
      assert(tmp14 == tmp14);
    const auto tmp15 = std::pow( tmp14, 2 );
      assert(tmp15 == tmp15);
    const auto tmp16 = 2 * tmp0;
      assert(tmp16 == tmp16);
    const auto tmp17 = tmp16 * tmp15;
      assert(tmp17 == tmp17);
    const auto tmp18 = tmp1[ 1 ] * tmp1[ 1 ];
      assert(tmp18 == tmp18);
    const auto tmp19 = tmp1[ 0 ] * tmp1[ 0 ];
      assert(tmp19 == tmp19);
    const auto tmp20 = tmp19 + tmp18;
      assert(tmp20 == tmp20);
    const auto tmp21 = tmp20 * tmp17;
      assert(tmp21 == tmp21);
    const auto tmp22 = tmp21 + tmp8;
      assert(tmp22 == tmp22);
    const auto tmp23 = tmp1[ 1 ] + tmp1[ 1 ];
      assert(tmp23 == tmp23);
    const auto tmp24 = tmp23 * tmp6;
      assert(tmp24 == tmp24);
    const auto tmp25 = -2 * tmp15;
      assert(tmp25 == tmp25);
    const auto tmp26 = tmp20 * tmp25;
      assert(tmp26 == tmp26);
    const auto tmp27 = tmp26 + tmp24;
      assert(tmp27 == tmp27);
    (result[ 0 ])[ 0 ] = tmp22;
    (result[ 0 ])[ 1 ] = tmp27;
  }

  template< class Point >
  void hessian ( const Point &x, typename FunctionSpaceType::HessianRangeType &result ) const
  {
    using std::cosh;
    using std::pow;
    using std::sinh;
    using std::tanh;
    double tmp0 = constant< 0 >();
      assert(tmp0 == tmp0);
    GlobalCoordinateType tmp1 = geometry().global( Dune::Fem::coordinate( x ) );
      assert(tmp1 == tmp1);
    const auto tmp2 = tmp1[ 0 ] * tmp0;
      assert(tmp2 == tmp2);
    const auto tmp3 = -1 * tmp1[ 1 ];
      assert(tmp3 == tmp3);
    const auto tmp4 = tmp3 + tmp2;
      assert(tmp4 == tmp4);
    const auto tmp5 = 2 * tmp4;
      assert(tmp5 == tmp5);
    const auto tmp6 = std::tanh( tmp5 );
      assert(tmp6 == tmp6);
    const auto tmp7 = 2 * tmp6;
      assert(tmp7 == tmp7);
    const auto tmp8 = 2.0 * tmp5;
      assert(tmp8 == tmp8);
    const auto tmp9 = std::cosh( tmp8 );
      assert(tmp9 == tmp9);
    const auto tmp10 = 1.0 + tmp9;
      assert(tmp10 == tmp10);
    const auto tmp11 = std::cosh( tmp5 );
      assert(tmp11 == tmp11);
    const auto tmp12 = 2.0 * tmp11;
      assert(tmp12 == tmp12);
    const auto tmp13 = tmp12 / tmp10;
      assert(tmp13 == tmp13);
    const auto tmp14 = std::pow( tmp13, 2 );
      assert(tmp14 == tmp14);
    const auto tmp15 = 2 * tmp0;
      assert(tmp15 == tmp15);
    const auto tmp16 = tmp15 * tmp14;
      assert(tmp16 == tmp16);
    const auto tmp17 = tmp1[ 0 ] + tmp1[ 0 ];
      assert(tmp17 == tmp17);
    const auto tmp18 = tmp17 * tmp16;
      assert(tmp18 == tmp18);
    const auto tmp19 = tmp18 + tmp7;
      assert(tmp19 == tmp19);
    const auto tmp20 = std::sinh( tmp5 );
      assert(tmp20 == tmp20);
    const auto tmp21 = tmp15 * tmp20;
      assert(tmp21 == tmp21);
    const auto tmp22 = 2.0 * tmp21;
      assert(tmp22 == tmp22);
    const auto tmp23 = std::sinh( tmp8 );
      assert(tmp23 == tmp23);
    const auto tmp24 = 2.0 * tmp15;
      assert(tmp24 == tmp24);
    const auto tmp25 = tmp24 * tmp23;
      assert(tmp25 == tmp25);
    const auto tmp26 = tmp25 * tmp13;
      assert(tmp26 == tmp26);
    const auto tmp27 = -1 * tmp26;
      assert(tmp27 == tmp27);
    const auto tmp28 = tmp27 + tmp22;
      assert(tmp28 == tmp28);
    const auto tmp29 = tmp28 / tmp10;
      assert(tmp29 == tmp29);
    const auto tmp30 = 2 * tmp29;
      assert(tmp30 == tmp30);
    const auto tmp31 = tmp30 * tmp13;
      assert(tmp31 == tmp31);
    const auto tmp32 = tmp31 * tmp15;
      assert(tmp32 == tmp32);
    const auto tmp33 = tmp1[ 1 ] * tmp1[ 1 ];
      assert(tmp33 == tmp33);
    const auto tmp34 = tmp1[ 0 ] * tmp1[ 0 ];
      assert(tmp34 == tmp34);
    const auto tmp35 = tmp34 + tmp33;
      assert(tmp35 == tmp35);
    const auto tmp36 = tmp35 * tmp32;
      assert(tmp36 == tmp36);
    const auto tmp37 = tmp36 + tmp18;
      assert(tmp37 == tmp37);
    const auto tmp38 = tmp37 + tmp19;
      assert(tmp38 == tmp38);
    const auto tmp39 = -2 * tmp14;
      assert(tmp39 == tmp39);
    const auto tmp40 = tmp17 * tmp39;
      assert(tmp40 == tmp40);
    const auto tmp41 = tmp1[ 1 ] + tmp1[ 1 ];
      assert(tmp41 == tmp41);
    const auto tmp42 = tmp41 * tmp16;
      assert(tmp42 == tmp42);
    const auto tmp43 = -2 * tmp20;
      assert(tmp43 == tmp43);
    const auto tmp44 = 2.0 * tmp43;
      assert(tmp44 == tmp44);
    const auto tmp45 = -4.0 * tmp23;
      assert(tmp45 == tmp45);
    const auto tmp46 = tmp45 * tmp13;
      assert(tmp46 == tmp46);
    const auto tmp47 = -1 * tmp46;
      assert(tmp47 == tmp47);
    const auto tmp48 = tmp47 + tmp44;
      assert(tmp48 == tmp48);
    const auto tmp49 = tmp48 / tmp10;
      assert(tmp49 == tmp49);
    const auto tmp50 = 2 * tmp49;
      assert(tmp50 == tmp50);
    const auto tmp51 = tmp50 * tmp13;
      assert(tmp51 == tmp51);
    const auto tmp52 = tmp51 * tmp15;
      assert(tmp52 == tmp52);
    const auto tmp53 = tmp35 * tmp52;
      assert(tmp53 == tmp53);
    const auto tmp54 = tmp53 + tmp42;
      assert(tmp54 == tmp54);
    const auto tmp55 = tmp54 + tmp40;
      assert(tmp55 == tmp55);
    const auto tmp56 = -2 * tmp31;
      assert(tmp56 == tmp56);
    const auto tmp57 = tmp35 * tmp56;
      assert(tmp57 == tmp57);
    const auto tmp58 = tmp57 + tmp40;
      assert(tmp58 == tmp58);
    const auto tmp59 = tmp58 + tmp42;
      assert(tmp59 == tmp59);
    const auto tmp60 = tmp41 * tmp39;
      assert(tmp60 == tmp60);
    const auto tmp61 = tmp60 + tmp7;
      assert(tmp61 == tmp61);
    const auto tmp62 = -2 * tmp51;
      assert(tmp62 == tmp62);
    const auto tmp63 = tmp35 * tmp62;
      assert(tmp63 == tmp63);
    const auto tmp64 = tmp63 + tmp60;
      assert(tmp64 == tmp64);
    const auto tmp65 = tmp64 + tmp61;
      assert(tmp65 == tmp65);
    ((result[ 0 ])[ 0 ])[ 0 ] = tmp38;
    ((result[ 0 ])[ 0 ])[ 1 ] = tmp55;
    ((result[ 0 ])[ 1 ])[ 0 ] = tmp59;
    ((result[ 0 ])[ 1 ])[ 1 ] = tmp65;
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

  const Contime &contime () const
  {
    return *std::get< 0 >( constants_ );
  }

  Contime &contime ()
  {
    return *std::get< 0 >( constants_ );
  }
  ConstantTupleType constants_;
  std::tuple<  > coefficients_;
};

} // namespace UFLLocalFunctions_37386eef2e0420e7169c7fc46256db10

PYBIND11_MODULE( localfunction_37386eef2e0420e7169c7fc46256db10_af122c1df944c95cd395ec0f91d0f970, module )
{
  typedef UFLLocalFunctions_37386eef2e0420e7169c7fc46256db10::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > > > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_37386eef2e0420e7169c7fc46256db10::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > > > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_37386eef2e0420e7169c7fc46256db10_af122c1df944c95cd395ec0f91d0f970.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > >>(gridView),name,order); } ), pybind11::keep_alive< 1, 2 >() );
      cls.def_property( "time", [] ( LocalFunctionType &self ) -> UFLLocalFunctions_37386eef2e0420e7169c7fc46256db10::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > > > >::Contime { return self.contime(); }, [] ( LocalFunctionType &self, const UFLLocalFunctions_37386eef2e0420e7169c7fc46256db10::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::Fem::AdaptiveLeafGridPart< Dune::ALUGrid< 2, 2, Dune::simplex > > > >::Contime &v ) { self.contime() = v; } );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
