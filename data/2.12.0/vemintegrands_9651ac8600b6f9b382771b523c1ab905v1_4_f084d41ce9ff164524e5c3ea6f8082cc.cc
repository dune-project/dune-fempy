#ifndef GuardIntegrands_9651ac8600b6f9b382771b523c1ab905v1_4
#define GuardIntegrands_9651ac8600b6f9b382771b523c1ab905v1_4
#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/python/grid/hierarchical.hh>
#include <dune/vem/misc/grid.hh>
#include <cmath>
#include <tuple>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/extensions.h>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/geometry/edgelength.hh>
#include <dune/fempy/geometry/celldiameter.hh>
#include <dune/vem/py/integrands.hh>

namespace Integrands_9651ac8600b6f9b382771b523c1ab905v1_4
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
    typedef std::tuple< Dune::Fem::ExplicitFieldVector< Dune::FieldMatrix< double, 2, 2 >, 1 > > DomainValueType;
    typedef std::tuple< Dune::FieldVector< double, 1 >, Dune::Fem::ExplicitFieldVector< Dune::FieldMatrix< double, 2, 2 >, 1 > > RangeValueType;
    static constexpr bool _nonlinear = false;

    bool nonlinear () const
    {
      return _nonlinear;
    }

    template< class Point >
    RangeValueType interior ( const Point &x, const DomainValueType &u ) const
    {
      using std::cos;
      using std::pow;
      using std::sin;
      GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
        assert(tmp0 == tmp0);
      const auto tmp1 = 6.283185307179586 * tmp0[ 1 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = std::sin( tmp1 );
        assert(tmp2 == tmp2);
      const auto tmp3 = std::pow( tmp2, 2 );
        assert(tmp3 == tmp3);
      const auto tmp4 = 6.283185307179586 * tmp0[ 0 ];
        assert(tmp4 == tmp4);
      const auto tmp5 = std::sin( tmp4 );
        assert(tmp5 == tmp5);
      const auto tmp6 = -1 * tmp5;
        assert(tmp6 == tmp6);
      const auto tmp7 = 6.283185307179586 * tmp6;
        assert(tmp7 == tmp7);
      const auto tmp8 = 6.283185307179586 * tmp7;
        assert(tmp8 == tmp8);
      const auto tmp9 = -1 * tmp8;
        assert(tmp9 == tmp9);
      const auto tmp10 = 6.283185307179586 * tmp9;
        assert(tmp10 == tmp10);
      const auto tmp11 = 6.283185307179586 * tmp10;
        assert(tmp11 == tmp11);
      const auto tmp12 = 2 * tmp11;
        assert(tmp12 == tmp12);
      const auto tmp13 = tmp12 * tmp5;
        assert(tmp13 == tmp13);
      const auto tmp14 = std::cos( tmp4 );
        assert(tmp14 == tmp14);
      const auto tmp15 = 6.283185307179586 * tmp14;
        assert(tmp15 == tmp15);
      const auto tmp16 = -1 * tmp15;
        assert(tmp16 == tmp16);
      const auto tmp17 = 6.283185307179586 * tmp16;
        assert(tmp17 == tmp17);
      const auto tmp18 = 6.283185307179586 * tmp17;
        assert(tmp18 == tmp18);
      const auto tmp19 = 2 * tmp18;
        assert(tmp19 == tmp19);
      const auto tmp20 = tmp19 * tmp15;
        assert(tmp20 == tmp20);
      const auto tmp21 = tmp20 + tmp13;
        assert(tmp21 == tmp21);
      const auto tmp22 = 2 * tmp8;
        assert(tmp22 == tmp22);
      const auto tmp23 = tmp22 * tmp8;
        assert(tmp23 == tmp23);
      const auto tmp24 = tmp23 + tmp20;
        assert(tmp24 == tmp24);
      const auto tmp25 = tmp24 + tmp21;
        assert(tmp25 == tmp25);
      const auto tmp26 = 2 * tmp15;
        assert(tmp26 == tmp26);
      const auto tmp27 = tmp18 * tmp26;
        assert(tmp27 == tmp27);
      const auto tmp28 = tmp23 + tmp27;
        assert(tmp28 == tmp28);
      const auto tmp29 = tmp28 + tmp24;
        assert(tmp29 == tmp29);
      const auto tmp30 = tmp29 + tmp25;
        assert(tmp30 == tmp30);
      const auto tmp31 = tmp30 * tmp3;
        assert(tmp31 == tmp31);
      const auto tmp32 = -1 * tmp2;
        assert(tmp32 == tmp32);
      const auto tmp33 = 6.283185307179586 * tmp32;
        assert(tmp33 == tmp33);
      const auto tmp34 = 6.283185307179586 * tmp33;
        assert(tmp34 == tmp34);
      const auto tmp35 = 2 * tmp34;
        assert(tmp35 == tmp35);
      const auto tmp36 = tmp35 * tmp2;
        assert(tmp36 == tmp36);
      const auto tmp37 = std::cos( tmp1 );
        assert(tmp37 == tmp37);
      const auto tmp38 = 6.283185307179586 * tmp37;
        assert(tmp38 == tmp38);
      const auto tmp39 = 2 * tmp38;
        assert(tmp39 == tmp39);
      const auto tmp40 = tmp39 * tmp38;
        assert(tmp40 == tmp40);
      const auto tmp41 = tmp40 + tmp36;
        assert(tmp41 == tmp41);
      const auto tmp42 = tmp22 * tmp5;
        assert(tmp42 == tmp42);
      const auto tmp43 = tmp26 * tmp15;
        assert(tmp43 == tmp43);
      const auto tmp44 = tmp43 + tmp42;
        assert(tmp44 == tmp44);
      const auto tmp45 = tmp44 * tmp41;
        assert(tmp45 == tmp45);
      const auto tmp46 = tmp45 + tmp31;
        assert(tmp46 == tmp46);
      const auto tmp47 = std::pow( tmp5, 2 );
        assert(tmp47 == tmp47);
      const auto tmp48 = -1 * tmp34;
        assert(tmp48 == tmp48);
      const auto tmp49 = 6.283185307179586 * tmp48;
        assert(tmp49 == tmp49);
      const auto tmp50 = 6.283185307179586 * tmp49;
        assert(tmp50 == tmp50);
      const auto tmp51 = 2 * tmp50;
        assert(tmp51 == tmp51);
      const auto tmp52 = tmp51 * tmp2;
        assert(tmp52 == tmp52);
      const auto tmp53 = -1 * tmp38;
        assert(tmp53 == tmp53);
      const auto tmp54 = 6.283185307179586 * tmp53;
        assert(tmp54 == tmp54);
      const auto tmp55 = 6.283185307179586 * tmp54;
        assert(tmp55 == tmp55);
      const auto tmp56 = 2 * tmp55;
        assert(tmp56 == tmp56);
      const auto tmp57 = tmp56 * tmp38;
        assert(tmp57 == tmp57);
      const auto tmp58 = tmp57 + tmp52;
        assert(tmp58 == tmp58);
      const auto tmp59 = tmp35 * tmp34;
        assert(tmp59 == tmp59);
      const auto tmp60 = tmp59 + tmp57;
        assert(tmp60 == tmp60);
      const auto tmp61 = tmp60 + tmp58;
        assert(tmp61 == tmp61);
      const auto tmp62 = tmp55 * tmp39;
        assert(tmp62 == tmp62);
      const auto tmp63 = tmp59 + tmp62;
        assert(tmp63 == tmp63);
      const auto tmp64 = tmp63 + tmp60;
        assert(tmp64 == tmp64);
      const auto tmp65 = tmp64 + tmp61;
        assert(tmp65 == tmp65);
      const auto tmp66 = tmp65 * tmp47;
        assert(tmp66 == tmp66);
      const auto tmp67 = tmp45 + tmp66;
        assert(tmp67 == tmp67);
      const auto tmp68 = tmp67 + tmp46;
        assert(tmp68 == tmp68);
      const auto tmp69 = -1 * tmp68;
        assert(tmp69 == tmp69);
      const auto tmp70 = std::get< 0 >( u );
        assert(tmp70 == tmp70);
      return RangeValueType{ { tmp69 }, { { { ((tmp70[ 0 ])[ 0 ])[ 0 ], ((tmp70[ 0 ])[ 0 ])[ 1 ] }, { ((tmp70[ 0 ])[ 1 ])[ 0 ], ((tmp70[ 0 ])[ 1 ])[ 1 ] } } } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      return [] ( const DomainValueType &phi ) {
          return RangeValueType{ { 0 }, { { { (((std::get< 0 >( phi ))[ 0 ])[ 0 ])[ 0 ], (((std::get< 0 >( phi ))[ 0 ])[ 0 ])[ 1 ] }, { (((std::get< 0 >( phi ))[ 0 ])[ 1 ])[ 0 ], (((std::get< 0 >( phi ))[ 0 ])[ 1 ])[ 1 ] } } } };
        };
    }
    typedef Dune::FieldVector< double, 1 >  RRangeType;
    typedef Dune::FieldMatrix< double, 1, GridPartType::dimension >  RJacobianRangeType;
    typedef Dune::Fem::BoundaryIdProvider< typename GridPartType::GridType > BoundaryIdProviderType;
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
          result[ 0 ] = 0.0;
        }
        break;
      case 2:
        {
          result[ 0 ] = 0.0;
        }
        break;
      case 3:
        {
          result[ 0 ] = 0.0;
        }
        break;
      case 4:
        {
          result[ 0 ] = 0.0;
        }
        break;
      default:
        {
          result = RRangeType( 0 );
        }
      }
    }

    template< class Point >
    void dDirichlet ( int bndId, const Point &x, RJacobianRangeType &result ) const
    {
      switch( bndId )
      {
      case 1:
        {
          (result[ 0 ])[ 0 ] = 0.0;
          (result[ 0 ])[ 1 ] = 0.0;
        }
        break;
      case 2:
        {
          (result[ 0 ])[ 0 ] = 0.0;
          (result[ 0 ])[ 1 ] = 0.0;
        }
        break;
      case 3:
        {
          (result[ 0 ])[ 0 ] = 0.0;
          (result[ 0 ])[ 1 ] = 0.0;
        }
        break;
      case 4:
        {
          (result[ 0 ])[ 0 ] = 0.0;
          (result[ 0 ])[ 1 ] = 0.0;
        }
        break;
      default:
        {
          result = RJacobianRangeType( 0 );
        }
      }
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

  public:
    template< class Point, class DRangeType >
    RRangeType hessStabilization ( const Point &x, const DRangeType &u ) const
    {
      RRangeType result;
      result[ 0 ] = 1;
      return result;
    }

    template< class Point, class DRangeType >
    RRangeType linHessStabilization ( const Point &x, const DRangeType &u ) const
    {
      RRangeType result;
      result[ 0 ] = 0.0;
      return result;
    }

    template< class Point, class DRangeType >
    RRangeType gradStabilization ( const Point &x, const DRangeType &u ) const
    {
      RRangeType result;
      result = RRangeType( (0) );
      return result;
    }

    template< class Point, class DRangeType >
    RRangeType linGradStabilization ( const Point &x, const DRangeType &u ) const
    {
      RRangeType result;
      result = RRangeType( (0) );
      return result;
    }

    template< class Point, class DRangeType >
    RRangeType massStabilization ( const Point &x, const DRangeType &u ) const
    {
      RRangeType result;
      result = RRangeType( (0) );
      return result;
    }

    template< class Point, class DRangeType >
    RRangeType linMassStabilization ( const Point &x, const DRangeType &u ) const
    {
      RRangeType result;
      result = RRangeType( (0) );
      return result;
    }
  };

} // namespace Integrands_9651ac8600b6f9b382771b523c1ab905v1_4

PYBIND11_MODULE( vemintegrands_9651ac8600b6f9b382771b523c1ab905v1_4_f084d41ce9ff164524e5c3ea6f8082cc, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView > GridPart;
  typedef Integrands_9651ac8600b6f9b382771b523c1ab905v1_4::Integrands< GridPart > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_9651ac8600b6f9b382771b523c1ab905v1_4::Integrands< GridPart >"), Dune::Python::IncludeFiles({"python/dune/generated/vemintegrands_9651ac8600b6f9b382771b523c1ab905v1_4_f084d41ce9ff164524e5c3ea6f8082cc.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] () { return new Integrands(); } ) );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return true;});
  }
}
#endif // GuardIntegrands_9651ac8600b6f9b382771b523c1ab905v1_4
