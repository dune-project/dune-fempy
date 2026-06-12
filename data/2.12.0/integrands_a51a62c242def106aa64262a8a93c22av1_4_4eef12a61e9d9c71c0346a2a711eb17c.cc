#ifndef GuardIntegrands_a51a62c242def106aa64262a8a93c22av1_4
#define GuardIntegrands_a51a62c242def106aa64262a8a93c22av1_4
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
#include <dune/fempy/geometry/celldiameter.hh>
#include <dune/fempy/py/integrands.hh>

namespace Integrands_a51a62c242def106aa64262a8a93c22av1_4
{




  // Integrands
  // ----------

  template< class GridPart >
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
    typedef std::tuple<> ConstantTupleType;
    typedef std::tuple<> CoefficientTupleType;
    static constexpr bool gridPartValid = true;
    template< std::size_t i >
    using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
    template< std::size_t i >
    using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;

    Integrands ( const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    {}

    void decouple ()
    {}

    bool init ( const EntityType &entity )
    {
      entity_ = entity;
      geometry_.emplace( this->entity().geometry() );
      return true;
    }

    void unbind ()
    {}

    bool init ( const IntersectionType &intersection )
    {
      intersection_ = intersection;
      return (intersection.boundary() && init( intersection.inside() ));
    }
    typedef std::tuple< Dune::FieldVector< double, 1 > > DomainValueType;
    typedef std::tuple< Dune::FieldVector< double, 1 > > RangeValueType;
    static constexpr bool _nonlinear = false;

    bool nonlinear () const
    {
      return _nonlinear;
    }

    template< class Point >
    RangeValueType interior ( const Point &x, const DomainValueType &u ) const
    {
      using std::sin;
      GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
        assert(tmp0 == tmp0);
      const auto tmp1 = 3.141592653589793 * tmp0[ 1 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = std::sin( tmp1 );
        assert(tmp2 == tmp2);
      const auto tmp3 = 3.141592653589793 * tmp0[ 0 ];
        assert(tmp3 == tmp3);
      const auto tmp4 = std::sin( tmp3 );
        assert(tmp4 == tmp4);
      const auto tmp5 = tmp4 * tmp2;
        assert(tmp5 == tmp5);
      const auto tmp6 = 3.141592653589793 * tmp0[ 2 ];
        assert(tmp6 == tmp6);
      const auto tmp7 = std::sin( tmp6 );
        assert(tmp7 == tmp7);
      const auto tmp8 = -1 * tmp7;
        assert(tmp8 == tmp8);
      const auto tmp9 = 3.141592653589793 * tmp8;
        assert(tmp9 == tmp9);
      const auto tmp10 = 3.141592653589793 * tmp9;
        assert(tmp10 == tmp10);
      const auto tmp11 = tmp10 * tmp5;
        assert(tmp11 == tmp11);
      const auto tmp12 = -1 * tmp4;
        assert(tmp12 == tmp12);
      const auto tmp13 = 3.141592653589793 * tmp12;
        assert(tmp13 == tmp13);
      const auto tmp14 = 3.141592653589793 * tmp13;
        assert(tmp14 == tmp14);
      const auto tmp15 = tmp14 * tmp2;
        assert(tmp15 == tmp15);
      const auto tmp16 = tmp15 * tmp7;
        assert(tmp16 == tmp16);
      const auto tmp17 = -1 * tmp2;
        assert(tmp17 == tmp17);
      const auto tmp18 = 3.141592653589793 * tmp17;
        assert(tmp18 == tmp18);
      const auto tmp19 = 3.141592653589793 * tmp18;
        assert(tmp19 == tmp19);
      const auto tmp20 = tmp19 * tmp4;
        assert(tmp20 == tmp20);
      const auto tmp21 = tmp20 * tmp7;
        assert(tmp21 == tmp21);
      const auto tmp22 = tmp21 + tmp16;
        assert(tmp22 == tmp22);
      const auto tmp23 = tmp22 + tmp11;
        assert(tmp23 == tmp23);
      const auto tmp24 = -1 * tmp23;
        assert(tmp24 == tmp24);
      const auto tmp25 = -1 * tmp24;
        assert(tmp25 == tmp25);
      const auto tmp26 = std::get< 0 >( u );
        assert(tmp26 == tmp26);
      const auto tmp27 = tmp26[ 0 ] + tmp25;
        assert(tmp27 == tmp27);
      return RangeValueType{ { tmp27 } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      return [] ( const DomainValueType &phi ) {
          return RangeValueType{ { (std::get< 0 >( phi ))[ 0 ] } };
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

    const EntityType &entity () const
    {
      return entity_;
    }

    const Geometry &geometry () const
    {
      return *geometry_;
    }
    EntityType entity_;
    IntersectionType intersection_;
    std::optional< Geometry > geometry_;
    ConstantTupleType constants_;
    std::tuple<  > coefficients_;
  };

} // namespace Integrands_a51a62c242def106aa64262a8a93c22av1_4

PYBIND11_MODULE( integrands_a51a62c242def106aa64262a8a93c22av1_4_4eef12a61e9d9c71c0346a2a711eb17c, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::YaspGrid< 3, Dune::EquidistantOffsetCoordinates< double, 3 > >::LeafGridView > GridPart;
  typedef Integrands_a51a62c242def106aa64262a8a93c22av1_4::Integrands< GridPart > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_a51a62c242def106aa64262a8a93c22av1_4::Integrands< GridPart >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_a51a62c242def106aa64262a8a93c22av1_4_4eef12a61e9d9c71c0346a2a711eb17c.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] () { return new Integrands(); } ) );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return false;});
  }
}
#endif // GuardIntegrands_a51a62c242def106aa64262a8a93c22av1_4
