#ifndef GUARD_7f2b0a245e7a0769f9f89b77b394214a
#define GUARD_7f2b0a245e7a0769f9f89b77b394214a

#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/python/grid/hierarchical.hh>
#include <dune/fem/function/localfunction/bindable.hh>
#include <dune/fem/common/intersectionside.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/extensions.h>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/common/exceptions.hh>
#include <dune/fempy/py/ufllocalfunction.hh>

namespace UFLLocalFunctions_7f2b0a245e7a0769f9f89b77b394214a
{

  // UFLLocalFunction
// ----------------

template< class GridPart >
struct UFLLocalFunction
  : public Dune::Fem::BindableGridFunctionWithSpace<GridPart,Dune::Dim<2>>
{
  typedef GridPart GridPartType;
  typedef typename GridPartType::GridViewType GridView;
  typedef typename GridView::ctype ctype;
  typedef Dune::Fem::BindableGridFunctionWithSpace<GridPart,Dune::Dim<2>> BaseType;
  typedef Dune::Fem::GridFunctionSpace<GridPartType,Dune::Dim<2>> FunctionSpaceType;
  typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
  typedef typename GridPartType::IntersectionType IntersectionType;
  typedef typename EntityType::Geometry Geometry;
  typedef typename Geometry::GlobalCoordinate GlobalCoordinateType;
  typedef Dune::Fem::IntersectionSide Side;
  typedef std::tuple<> ConstantTupleType;
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
  {}

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
    using std::pow;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
    const auto tmp1 = -1 + tmp0[ 0 ];
    const auto tmp2 = std::pow( tmp1, 2 );
    const auto tmp3 = std::pow( tmp0[ 0 ], 2 );
    const auto tmp4 = tmp3 * tmp2;
    const auto tmp5 = 12 * tmp0[ 1 ];
    const auto tmp6 = -6 + tmp5;
    const auto tmp7 = tmp6 * tmp4;
    const auto tmp8 = -1 + tmp0[ 1 ];
    const auto tmp9 = tmp0[ 1 ] * tmp8;
    const auto tmp10 = 2 * tmp0[ 1 ];
    const auto tmp11 = -1 + tmp10;
    const auto tmp12 = tmp11 * tmp9;
    const auto tmp13 = 12 * tmp3;
    const auto tmp14 = 12 * tmp0[ 0 ];
    const auto tmp15 = -1 * tmp14;
    const auto tmp16 = tmp15 + tmp13;
    const auto tmp17 = 2 + tmp16;
    const auto tmp18 = tmp17 * tmp12;
    const auto tmp19 = tmp18 + tmp7;
    const auto tmp20 = 256 * tmp19;
    const auto tmp21 = -0.5 + tmp0[ 1 ];
    const auto tmp22 = tmp21 + tmp20;
    const auto tmp23 = std::pow( tmp8, 2 );
    const auto tmp24 = std::pow( tmp0[ 1 ], 2 );
    const auto tmp25 = tmp24 * tmp23;
    const auto tmp26 = -6 + tmp14;
    const auto tmp27 = tmp26 * tmp25;
    const auto tmp28 = tmp0[ 0 ] * tmp1;
    const auto tmp29 = 2 * tmp0[ 0 ];
    const auto tmp30 = -1 + tmp29;
    const auto tmp31 = tmp30 * tmp28;
    const auto tmp32 = 12 * tmp24;
    const auto tmp33 = -1 * tmp5;
    const auto tmp34 = tmp33 + tmp32;
    const auto tmp35 = 2 + tmp34;
    const auto tmp36 = tmp35 * tmp31;
    const auto tmp37 = tmp36 + tmp27;
    const auto tmp38 = 256 * tmp37;
    const auto tmp39 = -1 * tmp38;
    const auto tmp40 = -0.5 + tmp0[ 0 ];
    const auto tmp41 = tmp40 + tmp39;
    result[ 0 ] = tmp22;
    result[ 1 ] = tmp41;
  }

  template< class Point >
  void jacobian ( const Point &x, typename FunctionSpaceType::JacobianRangeType &result ) const
  {
    using std::pow;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
    const auto tmp1 = -1 + tmp0[ 1 ];
    const auto tmp2 = tmp0[ 1 ] * tmp1;
    const auto tmp3 = 2 * tmp0[ 1 ];
    const auto tmp4 = -1 + tmp3;
    const auto tmp5 = tmp4 * tmp2;
    const auto tmp6 = 2 * tmp0[ 0 ];
    const auto tmp7 = 12 * tmp6;
    const auto tmp8 = -12 + tmp7;
    const auto tmp9 = tmp8 * tmp5;
    const auto tmp10 = -1 + tmp0[ 0 ];
    const auto tmp11 = std::pow( tmp10, 2 );
    const auto tmp12 = tmp6 * tmp11;
    const auto tmp13 = std::pow( tmp0[ 0 ], 2 );
    const auto tmp14 = 2 * tmp10;
    const auto tmp15 = tmp14 * tmp13;
    const auto tmp16 = tmp15 + tmp12;
    const auto tmp17 = 12 * tmp0[ 1 ];
    const auto tmp18 = -6 + tmp17;
    const auto tmp19 = tmp18 * tmp16;
    const auto tmp20 = tmp19 + tmp9;
    const auto tmp21 = 256 * tmp20;
    const auto tmp22 = tmp13 * tmp11;
    const auto tmp23 = 12 * tmp22;
    const auto tmp24 = 2 * tmp2;
    const auto tmp25 = tmp0[ 1 ] + tmp1;
    const auto tmp26 = tmp25 * tmp4;
    const auto tmp27 = tmp26 + tmp24;
    const auto tmp28 = 12 * tmp13;
    const auto tmp29 = 12 * tmp0[ 0 ];
    const auto tmp30 = -1 * tmp29;
    const auto tmp31 = tmp30 + tmp28;
    const auto tmp32 = 2 + tmp31;
    const auto tmp33 = tmp32 * tmp27;
    const auto tmp34 = tmp33 + tmp23;
    const auto tmp35 = 256 * tmp34;
    const auto tmp36 = 1 + tmp35;
    const auto tmp37 = std::pow( tmp1, 2 );
    const auto tmp38 = std::pow( tmp0[ 1 ], 2 );
    const auto tmp39 = tmp38 * tmp37;
    const auto tmp40 = 12 * tmp39;
    const auto tmp41 = tmp0[ 0 ] * tmp10;
    const auto tmp42 = 2 * tmp41;
    const auto tmp43 = -1 + tmp6;
    const auto tmp44 = tmp0[ 0 ] + tmp10;
    const auto tmp45 = tmp44 * tmp43;
    const auto tmp46 = tmp45 + tmp42;
    const auto tmp47 = 12 * tmp38;
    const auto tmp48 = -1 * tmp17;
    const auto tmp49 = tmp48 + tmp47;
    const auto tmp50 = 2 + tmp49;
    const auto tmp51 = tmp50 * tmp46;
    const auto tmp52 = tmp51 + tmp40;
    const auto tmp53 = 256 * tmp52;
    const auto tmp54 = -1 * tmp53;
    const auto tmp55 = 1 + tmp54;
    const auto tmp56 = tmp43 * tmp41;
    const auto tmp57 = 12 * tmp3;
    const auto tmp58 = -12 + tmp57;
    const auto tmp59 = tmp58 * tmp56;
    const auto tmp60 = tmp3 * tmp37;
    const auto tmp61 = 2 * tmp1;
    const auto tmp62 = tmp61 * tmp38;
    const auto tmp63 = tmp62 + tmp60;
    const auto tmp64 = -6 + tmp29;
    const auto tmp65 = tmp64 * tmp63;
    const auto tmp66 = tmp65 + tmp59;
    const auto tmp67 = 256 * tmp66;
    const auto tmp68 = -1 * tmp67;
    (result[ 0 ])[ 0 ] = tmp21;
    (result[ 0 ])[ 1 ] = tmp36;
    (result[ 1 ])[ 0 ] = tmp55;
    (result[ 1 ])[ 1 ] = tmp68;
  }

  template< class Point >
  void hessian ( const Point &x, typename FunctionSpaceType::HessianRangeType &result ) const
  {
    using std::pow;
    GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
    const auto tmp1 = -1 + tmp0[ 1 ];
    const auto tmp2 = tmp0[ 1 ] * tmp1;
    const auto tmp3 = 2 * tmp0[ 1 ];
    const auto tmp4 = -1 + tmp3;
    const auto tmp5 = tmp4 * tmp2;
    const auto tmp6 = 24 * tmp5;
    const auto tmp7 = 12 * tmp0[ 1 ];
    const auto tmp8 = -6 + tmp7;
    const auto tmp9 = -1 + tmp0[ 0 ];
    const auto tmp10 = std::pow( tmp9, 2 );
    const auto tmp11 = 2 * tmp10;
    const auto tmp12 = 2 * tmp9;
    const auto tmp13 = 2 * tmp0[ 0 ];
    const auto tmp14 = tmp13 * tmp12;
    const auto tmp15 = tmp14 + tmp11;
    const auto tmp16 = std::pow( tmp0[ 0 ], 2 );
    const auto tmp17 = 2 * tmp16;
    const auto tmp18 = tmp14 + tmp17;
    const auto tmp19 = tmp18 + tmp15;
    const auto tmp20 = tmp19 * tmp8;
    const auto tmp21 = tmp20 + tmp6;
    const auto tmp22 = 256 * tmp21;
    const auto tmp23 = tmp13 * tmp10;
    const auto tmp24 = tmp12 * tmp16;
    const auto tmp25 = tmp24 + tmp23;
    const auto tmp26 = 12 * tmp25;
    const auto tmp27 = 2 * tmp2;
    const auto tmp28 = tmp0[ 1 ] + tmp1;
    const auto tmp29 = tmp28 * tmp4;
    const auto tmp30 = tmp29 + tmp27;
    const auto tmp31 = 12 * tmp13;
    const auto tmp32 = -12 + tmp31;
    const auto tmp33 = tmp32 * tmp30;
    const auto tmp34 = tmp33 + tmp26;
    const auto tmp35 = 256 * tmp34;
    const auto tmp36 = 2 * tmp28;
    const auto tmp37 = 2 * tmp4;
    const auto tmp38 = tmp36 + tmp37;
    const auto tmp39 = tmp38 + tmp36;
    const auto tmp40 = 12 * tmp16;
    const auto tmp41 = 12 * tmp0[ 0 ];
    const auto tmp42 = -1 * tmp41;
    const auto tmp43 = tmp42 + tmp40;
    const auto tmp44 = 2 + tmp43;
    const auto tmp45 = tmp44 * tmp39;
    const auto tmp46 = 256 * tmp45;
    const auto tmp47 = tmp0[ 0 ] + tmp9;
    const auto tmp48 = 2 * tmp47;
    const auto tmp49 = -1 + tmp13;
    const auto tmp50 = 2 * tmp49;
    const auto tmp51 = tmp48 + tmp50;
    const auto tmp52 = tmp51 + tmp48;
    const auto tmp53 = std::pow( tmp0[ 1 ], 2 );
    const auto tmp54 = 12 * tmp53;
    const auto tmp55 = -1 * tmp7;
    const auto tmp56 = tmp55 + tmp54;
    const auto tmp57 = 2 + tmp56;
    const auto tmp58 = tmp57 * tmp52;
    const auto tmp59 = 256 * tmp58;
    const auto tmp60 = -1 * tmp59;
    const auto tmp61 = std::pow( tmp1, 2 );
    const auto tmp62 = tmp3 * tmp61;
    const auto tmp63 = 2 * tmp1;
    const auto tmp64 = tmp63 * tmp53;
    const auto tmp65 = tmp64 + tmp62;
    const auto tmp66 = 12 * tmp65;
    const auto tmp67 = tmp0[ 0 ] * tmp9;
    const auto tmp68 = 2 * tmp67;
    const auto tmp69 = tmp47 * tmp49;
    const auto tmp70 = tmp69 + tmp68;
    const auto tmp71 = 12 * tmp3;
    const auto tmp72 = -12 + tmp71;
    const auto tmp73 = tmp72 * tmp70;
    const auto tmp74 = tmp73 + tmp66;
    const auto tmp75 = 256 * tmp74;
    const auto tmp76 = -1 * tmp75;
    const auto tmp77 = tmp49 * tmp67;
    const auto tmp78 = 24 * tmp77;
    const auto tmp79 = -6 + tmp41;
    const auto tmp80 = 2 * tmp61;
    const auto tmp81 = tmp3 * tmp63;
    const auto tmp82 = tmp81 + tmp80;
    const auto tmp83 = 2 * tmp53;
    const auto tmp84 = tmp81 + tmp83;
    const auto tmp85 = tmp84 + tmp82;
    const auto tmp86 = tmp85 * tmp79;
    const auto tmp87 = tmp86 + tmp78;
    const auto tmp88 = 256 * tmp87;
    const auto tmp89 = -1 * tmp88;
    ((result[ 0 ])[ 0 ])[ 0 ] = tmp22;
    ((result[ 0 ])[ 0 ])[ 1 ] = tmp35;
    ((result[ 0 ])[ 1 ])[ 0 ] = tmp35;
    ((result[ 0 ])[ 1 ])[ 1 ] = tmp46;
    ((result[ 1 ])[ 0 ])[ 0 ] = tmp60;
    ((result[ 1 ])[ 0 ])[ 1 ] = tmp76;
    ((result[ 1 ])[ 1 ])[ 0 ] = tmp76;
    ((result[ 1 ])[ 1 ])[ 1 ] = tmp89;
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
  ConstantTupleType constants_;
  std::tuple<  > coefficients_;
};

} // namespace UFLLocalFunctions_7f2b0a245e7a0769f9f89b77b394214a

PYBIND11_MODULE( localfunction_7f2b0a245e7a0769f9f89b77b394214a_19659fe2f2f74bcfb5d532f1a0d38be9, module )
{
  typedef UFLLocalFunctions_7f2b0a245e7a0769f9f89b77b394214a::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView > > LocalFunctionType;
  if constexpr( LocalFunctionType::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("UFLLocalFunctions_7f2b0a245e7a0769f9f89b77b394214a::UFLLocalFunction< typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView > >"), Dune::Python::IncludeFiles({"python/dune/generated/localfunction_7f2b0a245e7a0769f9f89b77b394214a_19659fe2f2f74bcfb5d532f1a0d38be9.cc"})).first;
      Dune::FemPy::registerUFLLocalFunction( module, cls );
      cls.def( pybind11::init( [] ( pybind11::object gridView, const std::string &name, int order ) {return new LocalFunctionType( Dune::FemPy::gridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>(gridView),name,order); } ), pybind11::keep_alive< 1, 2 >() );
      cls.def_property_readonly( "virtualized", [] ( LocalFunctionType& ) -> bool { return true;});
  }
}

#endif
