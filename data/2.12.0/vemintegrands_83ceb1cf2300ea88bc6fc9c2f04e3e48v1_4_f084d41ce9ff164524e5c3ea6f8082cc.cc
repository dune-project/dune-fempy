#ifndef GuardIntegrands_83ceb1cf2300ea88bc6fc9c2f04e3e48v1_4
#define GuardIntegrands_83ceb1cf2300ea88bc6fc9c2f04e3e48v1_4
#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/python/grid/hierarchical.hh>
#include <dune/vem/misc/grid.hh>
#include <tuple>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/extensions.h>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/geometry/edgelength.hh>
#include <dune/fempy/geometry/celldiameter.hh>
#include <dune/vem/py/integrands.hh>

namespace Integrands_83ceb1cf2300ea88bc6fc9c2f04e3e48v1_4
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
      return true;
    }

    void unbind ()
    {}

    bool init ( const IntersectionType &intersection )
    {
      intersection_ = intersection;
      return (intersection.boundary() && init( intersection.inside() ));
    }
    typedef std::tuple< Dune::FieldMatrix< double, 2, 2 > > DomainValueType;
    typedef std::tuple< Dune::FieldVector< double, 2 >, Dune::FieldMatrix< double, 2, 2 > > RangeValueType;
    static constexpr bool _nonlinear = false;

    bool nonlinear () const
    {
      return _nonlinear;
    }

    template< class Point >
    RangeValueType interior ( const Point &x, const DomainValueType &u ) const
    {
      const auto tmp0 = std::get< 0 >( u );
        assert(tmp0 == tmp0);
      const auto tmp1 = (tmp0[ 0 ])[ 0 ] + (tmp0[ 0 ])[ 0 ];
        assert(tmp1 == tmp1);
      const auto tmp2 = 0.5 * tmp1;
        assert(tmp2 == tmp2);
      const auto tmp3 = 2 * tmp2;
        assert(tmp3 == tmp3);
      const auto tmp4 = (tmp0[ 0 ])[ 0 ] + (tmp0[ 1 ])[ 1 ];
        assert(tmp4 == tmp4);
      const auto tmp5 = 1.25 * tmp4;
        assert(tmp5 == tmp5);
      const auto tmp6 = tmp5 + tmp3;
        assert(tmp6 == tmp6);
      const auto tmp7 = (tmp0[ 0 ])[ 1 ] + (tmp0[ 1 ])[ 0 ];
        assert(tmp7 == tmp7);
      const auto tmp8 = 0.5 * tmp7;
        assert(tmp8 == tmp8);
      const auto tmp9 = 2 * tmp8;
        assert(tmp9 == tmp9);
      const auto tmp10 = 0.5 * tmp9;
        assert(tmp10 == tmp10);
      const auto tmp11 = tmp10 + tmp10;
        assert(tmp11 == tmp11);
      const auto tmp12 = (tmp0[ 1 ])[ 1 ] + (tmp0[ 1 ])[ 1 ];
        assert(tmp12 == tmp12);
      const auto tmp13 = 0.5 * tmp12;
        assert(tmp13 == tmp13);
      const auto tmp14 = 2 * tmp13;
        assert(tmp14 == tmp14);
      const auto tmp15 = tmp5 + tmp14;
        assert(tmp15 == tmp15);
      return RangeValueType{ { 0.0, 0.016000000000000004 }, { { tmp6, tmp11 }, { tmp11, tmp15 } } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      return [] ( const DomainValueType &phi ) {
          return RangeValueType{ { 0, 0 }, { { 3.25 * ((std::get< 0 >( phi ))[ 0 ])[ 0 ] + 1.25 * ((std::get< 0 >( phi ))[ 1 ])[ 1 ], 1.0 * ((std::get< 0 >( phi ))[ 0 ])[ 1 ] + 1.0 * ((std::get< 0 >( phi ))[ 1 ])[ 0 ] }, { 1.0 * ((std::get< 0 >( phi ))[ 0 ])[ 1 ] + 1.0 * ((std::get< 0 >( phi ))[ 1 ])[ 0 ], 1.25 * ((std::get< 0 >( phi ))[ 0 ])[ 0 ] + 3.25 * ((std::get< 0 >( phi ))[ 1 ])[ 1 ] } } };
        };
    }
    typedef Dune::FieldVector< double, 2 >  RRangeType;
    typedef Dune::FieldMatrix< double, 2, GridPartType::dimension >  RJacobianRangeType;
    typedef Dune::Fem::BoundaryIdProvider< typename GridPartType::GridType > BoundaryIdProviderType;
    typedef std::array<int,2> DirichletComponentType;

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
      default:
        {
          int domainId;
          {
              auto tmp0 = intersection.geometry().center( );
                assert(tmp0 == tmp0);
              domainId = tmp0[ 0 ] < 1e-10;
              if (domainId)
              {
                std::fill( dirichletComponent.begin(), dirichletComponent.end(), 2 );
                return true;
              }
          }
          return false;
        }
      }
    }

    template< class Point >
    void dirichlet ( int bndId, const Point &x, RRangeType &result ) const
    {
      switch( bndId )
      {
      case 2:
        {
          result[ 0 ] = 0.0;
          result[ 1 ] = 0.0;
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
      case 2:
        {
          (result[ 0 ])[ 0 ] = 0.0;
          (result[ 0 ])[ 1 ] = 0.0;
          (result[ 1 ])[ 0 ] = 0.0;
          (result[ 1 ])[ 1 ] = 0.0;
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
    EntityType entity_;
    IntersectionType intersection_;
    ConstantTupleType constants_;
    std::tuple<  > coefficients_;

  public:
    template< class Point, class DRangeType >
    RRangeType hessStabilization ( const Point &x, const DRangeType &u ) const
    {
      RRangeType result;
      result = RRangeType( (0) );
      return result;
    }

    template< class Point, class DRangeType >
    RRangeType linHessStabilization ( const Point &x, const DRangeType &u ) const
    {
      RRangeType result;
      result = RRangeType( (0) );
      return result;
    }

    template< class Point, class DRangeType >
    RRangeType gradStabilization ( const Point &x, const DRangeType &u ) const
    {
      RRangeType result;
      result[ 0 ] = 3.25;
      result[ 1 ] = 3.25;
      return result;
    }

    template< class Point, class DRangeType >
    RRangeType linGradStabilization ( const Point &x, const DRangeType &u ) const
    {
      RRangeType result;
      result[ 0 ] = 0.0;
      result[ 1 ] = 0.0;
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

} // namespace Integrands_83ceb1cf2300ea88bc6fc9c2f04e3e48v1_4

PYBIND11_MODULE( vemintegrands_83ceb1cf2300ea88bc6fc9c2f04e3e48v1_4_f084d41ce9ff164524e5c3ea6f8082cc, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::Vem::Grid<2,2>::LeafGridView > GridPart;
  typedef Integrands_83ceb1cf2300ea88bc6fc9c2f04e3e48v1_4::Integrands< GridPart > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_83ceb1cf2300ea88bc6fc9c2f04e3e48v1_4::Integrands< GridPart >"), Dune::Python::IncludeFiles({"python/dune/generated/vemintegrands_83ceb1cf2300ea88bc6fc9c2f04e3e48v1_4_f084d41ce9ff164524e5c3ea6f8082cc.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] () { return new Integrands(); } ) );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return true;});
  }
}
#endif // GuardIntegrands_83ceb1cf2300ea88bc6fc9c2f04e3e48v1_4
