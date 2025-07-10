#ifndef GuardIntegrands_5dbc63b9691d63201b9a98e2aa8071bdv1_3
#define GuardIntegrands_5dbc63b9691d63201b9a98e2aa8071bdv1_3
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
#include <dune/fempy/py/integrands.hh>

namespace Integrands_5dbc63b9691d63201b9a98e2aa8071bdv1_3
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
    typedef std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 > > DomainValueType;
    typedef std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 > > RangeValueType;
    static constexpr bool _nonlinear = false;

    bool nonlinear () const
    {
      return _nonlinear;
    }

    template< class Point >
    RangeValueType interior ( const Point &x, const DomainValueType &u ) const
    {
      using std::cos;
      using std::sin;
      GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
      const auto tmp1 = tmp0[ 1 ] * tmp0[ 1 ];
      const auto tmp2 = tmp0[ 0 ] * tmp0[ 0 ];
      const auto tmp3 = tmp2 + tmp1;
      const auto tmp4 = std::sin( tmp3 );
      const auto tmp5 = 1 + tmp4;
      const auto tmp6 = std::get< 0 >( u );
      const auto tmp7 = tmp6[ 0 ] * tmp5;
      const auto tmp8 = 3.141592653589793 * tmp0[ 0 ];
      const auto tmp9 = tmp0[ 1 ] * tmp8;
      const auto tmp10 = std::cos( tmp9 );
      const auto tmp11 = tmp0[ 0 ] * tmp0[ 1 ];
      const auto tmp12 = tmp11 * tmp10;
      const auto tmp13 = tmp5 * tmp12;
      const auto tmp14 = -1 * tmp4;
      const auto tmp15 = tmp0[ 1 ] + tmp0[ 1 ];
      const auto tmp16 = tmp15 * tmp14;
      const auto tmp17 = 0.9 * tmp16;
      const auto tmp18 = -1 * tmp17;
      const auto tmp19 = tmp0[ 0 ] * tmp10;
      const auto tmp20 = std::sin( tmp9 );
      const auto tmp21 = -1 * tmp20;
      const auto tmp22 = tmp8 * tmp21;
      const auto tmp23 = tmp11 * tmp22;
      const auto tmp24 = tmp23 + tmp19;
      const auto tmp25 = tmp24 * tmp18;
      const auto tmp26 = std::cos( tmp3 );
      const auto tmp27 = 0.9 * tmp26;
      const auto tmp28 = -1 * tmp27;
      const auto tmp29 = 1 + tmp28;
      const auto tmp30 = tmp0[ 0 ] * tmp22;
      const auto tmp31 = tmp8 * tmp10;
      const auto tmp32 = -1 * tmp31;
      const auto tmp33 = tmp8 * tmp32;
      const auto tmp34 = tmp11 * tmp33;
      const auto tmp35 = tmp34 + tmp30;
      const auto tmp36 = tmp35 + tmp30;
      const auto tmp37 = tmp36 * tmp29;
      const auto tmp38 = tmp37 + tmp25;
      const auto tmp39 = tmp0[ 0 ] + tmp0[ 0 ];
      const auto tmp40 = tmp39 * tmp14;
      const auto tmp41 = 0.9 * tmp40;
      const auto tmp42 = -1 * tmp41;
      const auto tmp43 = tmp0[ 1 ] * tmp10;
      const auto tmp44 = 3.141592653589793 * tmp0[ 1 ];
      const auto tmp45 = tmp44 * tmp21;
      const auto tmp46 = tmp11 * tmp45;
      const auto tmp47 = tmp46 + tmp43;
      const auto tmp48 = tmp47 * tmp42;
      const auto tmp49 = tmp0[ 1 ] * tmp45;
      const auto tmp50 = tmp44 * tmp10;
      const auto tmp51 = -1 * tmp50;
      const auto tmp52 = tmp44 * tmp51;
      const auto tmp53 = tmp11 * tmp52;
      const auto tmp54 = tmp53 + tmp49;
      const auto tmp55 = tmp54 + tmp49;
      const auto tmp56 = tmp55 * tmp29;
      const auto tmp57 = tmp56 + tmp48;
      const auto tmp58 = tmp57 + tmp38;
      const auto tmp59 = -1 * tmp58;
      const auto tmp60 = tmp59 + tmp13;
      const auto tmp61 = -1 * tmp60;
      const auto tmp62 = tmp61 + tmp7;
      const auto tmp63 = std::get< 1 >( u );
      const auto tmp64 = (tmp63[ 0 ])[ 0 ] * tmp29;
      const auto tmp65 = (tmp63[ 0 ])[ 1 ] * tmp29;
      return RangeValueType{ { tmp62 }, { { tmp64, tmp65 } } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      using std::cos;
      using std::sin;
      GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
      const auto tmp1 = tmp0[ 1 ] * tmp0[ 1 ];
      const auto tmp2 = tmp0[ 0 ] * tmp0[ 0 ];
      const auto tmp3 = tmp2 + tmp1;
      const auto tmp4 = std::sin( tmp3 );
      const auto tmp5 = 1 + tmp4;
      const auto tmp6 = std::cos( tmp3 );
      const auto tmp7 = 0.9 * tmp6;
      const auto tmp8 = -1 * tmp7;
      const auto tmp9 = 1 + tmp8;
      return [ tmp5, tmp9 ] ( const DomainValueType &phi ) {
          return RangeValueType{ { tmp5 * (std::get< 0 >( phi ))[ 0 ] }, { { tmp9 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ], tmp9 * ((std::get< 1 >( phi ))[ 0 ])[ 1 ] } } };
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
          using std::cos;
          auto tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
          const auto tmp1 = 3.141592653589793 * tmp0[ 0 ];
          const auto tmp2 = tmp0[ 1 ] * tmp1;
          const auto tmp3 = std::cos( tmp2 );
          const auto tmp4 = tmp0[ 0 ] * tmp0[ 1 ];
          const auto tmp5 = tmp4 * tmp3;
          result[ 0 ] = tmp5;
        }
        break;
      case 2:
        {
          using std::cos;
          auto tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
          const auto tmp1 = 3.141592653589793 * tmp0[ 0 ];
          const auto tmp2 = tmp0[ 1 ] * tmp1;
          const auto tmp3 = std::cos( tmp2 );
          const auto tmp4 = tmp0[ 0 ] * tmp0[ 1 ];
          const auto tmp5 = tmp4 * tmp3;
          result[ 0 ] = tmp5;
        }
        break;
      case 3:
        {
          using std::cos;
          auto tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
          const auto tmp1 = 3.141592653589793 * tmp0[ 0 ];
          const auto tmp2 = tmp0[ 1 ] * tmp1;
          const auto tmp3 = std::cos( tmp2 );
          const auto tmp4 = tmp0[ 0 ] * tmp0[ 1 ];
          const auto tmp5 = tmp4 * tmp3;
          result[ 0 ] = tmp5;
        }
        break;
      case 4:
        {
          using std::cos;
          auto tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
          const auto tmp1 = 3.141592653589793 * tmp0[ 0 ];
          const auto tmp2 = tmp0[ 1 ] * tmp1;
          const auto tmp3 = std::cos( tmp2 );
          const auto tmp4 = tmp0[ 0 ] * tmp0[ 1 ];
          const auto tmp5 = tmp4 * tmp3;
          result[ 0 ] = tmp5;
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
          using std::cos;
          using std::sin;
          auto tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
          const auto tmp1 = 3.141592653589793 * tmp0[ 0 ];
          const auto tmp2 = tmp0[ 1 ] * tmp1;
          const auto tmp3 = std::cos( tmp2 );
          const auto tmp4 = tmp0[ 1 ] * tmp3;
          const auto tmp5 = std::sin( tmp2 );
          const auto tmp6 = -1 * tmp5;
          const auto tmp7 = 3.141592653589793 * tmp0[ 1 ];
          const auto tmp8 = tmp7 * tmp6;
          const auto tmp9 = tmp0[ 0 ] * tmp0[ 1 ];
          const auto tmp10 = tmp9 * tmp8;
          const auto tmp11 = tmp10 + tmp4;
          const auto tmp12 = tmp0[ 0 ] * tmp3;
          const auto tmp13 = tmp1 * tmp6;
          const auto tmp14 = tmp9 * tmp13;
          const auto tmp15 = tmp14 + tmp12;
          (result[ 0 ])[ 0 ] = tmp11;
          (result[ 0 ])[ 1 ] = tmp15;
        }
        break;
      case 2:
        {
          using std::cos;
          using std::sin;
          auto tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
          const auto tmp1 = 3.141592653589793 * tmp0[ 0 ];
          const auto tmp2 = tmp0[ 1 ] * tmp1;
          const auto tmp3 = std::cos( tmp2 );
          const auto tmp4 = tmp0[ 1 ] * tmp3;
          const auto tmp5 = std::sin( tmp2 );
          const auto tmp6 = -1 * tmp5;
          const auto tmp7 = 3.141592653589793 * tmp0[ 1 ];
          const auto tmp8 = tmp7 * tmp6;
          const auto tmp9 = tmp0[ 0 ] * tmp0[ 1 ];
          const auto tmp10 = tmp9 * tmp8;
          const auto tmp11 = tmp10 + tmp4;
          const auto tmp12 = tmp0[ 0 ] * tmp3;
          const auto tmp13 = tmp1 * tmp6;
          const auto tmp14 = tmp9 * tmp13;
          const auto tmp15 = tmp14 + tmp12;
          (result[ 0 ])[ 0 ] = tmp11;
          (result[ 0 ])[ 1 ] = tmp15;
        }
        break;
      case 3:
        {
          using std::cos;
          using std::sin;
          auto tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
          const auto tmp1 = 3.141592653589793 * tmp0[ 0 ];
          const auto tmp2 = tmp0[ 1 ] * tmp1;
          const auto tmp3 = std::cos( tmp2 );
          const auto tmp4 = tmp0[ 1 ] * tmp3;
          const auto tmp5 = std::sin( tmp2 );
          const auto tmp6 = -1 * tmp5;
          const auto tmp7 = 3.141592653589793 * tmp0[ 1 ];
          const auto tmp8 = tmp7 * tmp6;
          const auto tmp9 = tmp0[ 0 ] * tmp0[ 1 ];
          const auto tmp10 = tmp9 * tmp8;
          const auto tmp11 = tmp10 + tmp4;
          const auto tmp12 = tmp0[ 0 ] * tmp3;
          const auto tmp13 = tmp1 * tmp6;
          const auto tmp14 = tmp9 * tmp13;
          const auto tmp15 = tmp14 + tmp12;
          (result[ 0 ])[ 0 ] = tmp11;
          (result[ 0 ])[ 1 ] = tmp15;
        }
        break;
      case 4:
        {
          using std::cos;
          using std::sin;
          auto tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
          const auto tmp1 = 3.141592653589793 * tmp0[ 0 ];
          const auto tmp2 = tmp0[ 1 ] * tmp1;
          const auto tmp3 = std::cos( tmp2 );
          const auto tmp4 = tmp0[ 1 ] * tmp3;
          const auto tmp5 = std::sin( tmp2 );
          const auto tmp6 = -1 * tmp5;
          const auto tmp7 = 3.141592653589793 * tmp0[ 1 ];
          const auto tmp8 = tmp7 * tmp6;
          const auto tmp9 = tmp0[ 0 ] * tmp0[ 1 ];
          const auto tmp10 = tmp9 * tmp8;
          const auto tmp11 = tmp10 + tmp4;
          const auto tmp12 = tmp0[ 0 ] * tmp3;
          const auto tmp13 = tmp1 * tmp6;
          const auto tmp14 = tmp9 * tmp13;
          const auto tmp15 = tmp14 + tmp12;
          (result[ 0 ])[ 0 ] = tmp11;
          (result[ 0 ])[ 1 ] = tmp15;
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
  };

} // namespace Integrands_5dbc63b9691d63201b9a98e2aa8071bdv1_3

PYBIND11_MODULE( integrands_5dbc63b9691d63201b9a98e2aa8071bdv1_3_f084d41ce9ff164524e5c3ea6f8082cc, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView > GridPart;
  typedef Integrands_5dbc63b9691d63201b9a98e2aa8071bdv1_3::Integrands< GridPart > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_5dbc63b9691d63201b9a98e2aa8071bdv1_3::Integrands< GridPart >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_5dbc63b9691d63201b9a98e2aa8071bdv1_3_f084d41ce9ff164524e5c3ea6f8082cc.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] () { return new Integrands(); } ) );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return true;});
  }
}
#endif // GuardIntegrands_5dbc63b9691d63201b9a98e2aa8071bdv1_3
