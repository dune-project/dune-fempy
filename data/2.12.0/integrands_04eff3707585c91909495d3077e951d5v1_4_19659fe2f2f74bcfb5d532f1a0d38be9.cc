#ifndef GuardIntegrands_04eff3707585c91909495d3077e951d5v1_4
#define GuardIntegrands_04eff3707585c91909495d3077e951d5v1_4
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

namespace Integrands_04eff3707585c91909495d3077e951d5v1_4
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
    typedef double Conmass;
    typedef std::tuple< std::shared_ptr< Conmass > > ConstantTupleType;
    template< std::size_t i >
    using ConstantsRangeType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;
    typedef std::tuple<> CoefficientTupleType;
    static constexpr bool gridPartValid = true;
    template< std::size_t i >
    using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
    template< std::size_t i >
    using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;

    Integrands ( const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
    {
      std::get< 0 >( constants_ ) = std::make_shared< Conmass >( (Conmass(0)) );
    }

    void decouple ()
    {
      std::get< 0 >( constants_ ) = std::make_shared< Conmass >( (constant<0>()) );
    }

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
      double tmp0 = constant< 0 >();
        assert(tmp0 == tmp0);
      const auto tmp1 = std::get< 0 >( u );
        assert(tmp1 == tmp1);
      const auto tmp2 = tmp1[ 0 ] * tmp0;
        assert(tmp2 == tmp2);
      GlobalCoordinateType tmp3 = geometry().global( Dune::Fem::coordinate( x ) );
        assert(tmp3 == tmp3);
      const auto tmp4 = 2 * tmp3[ 1 ];
        assert(tmp4 == tmp4);
      const auto tmp5 = 3 * tmp4;
        assert(tmp5 == tmp5);
      const auto tmp6 = -1 * tmp5;
        assert(tmp6 == tmp6);
      const auto tmp7 = 0.3333333333333333 * tmp6;
        assert(tmp7 == tmp7);
      const auto tmp8 = -1 * tmp7;
        assert(tmp8 == tmp8);
      const auto tmp9 = 1.0 + tmp8;
        assert(tmp9 == tmp9);
      const auto tmp10 = 2 * tmp3[ 0 ];
        assert(tmp10 == tmp10);
      const auto tmp11 = 3 * tmp10;
        assert(tmp11 == tmp11);
      const auto tmp12 = 0.3333333333333333 * tmp11;
        assert(tmp12 == tmp12);
      const auto tmp13 = -1 * tmp12;
        assert(tmp13 == tmp13);
      const auto tmp14 = 1.0 + tmp13;
        assert(tmp14 == tmp14);
      const auto tmp15 = tmp14 + tmp9;
        assert(tmp15 == tmp15);
      const auto tmp16 = -1 * tmp15;
        assert(tmp16 == tmp16);
      const auto tmp17 = -1 * tmp16;
        assert(tmp17 == tmp17);
      const auto tmp18 = tmp17 + tmp2;
        assert(tmp18 == tmp18);
      const auto tmp19 = std::get< 1 >( u );
        assert(tmp19 == tmp19);
      return RangeValueType{ { tmp18 }, { { (tmp19[ 0 ])[ 0 ], (tmp19[ 0 ])[ 1 ] } } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      double tmp0 = constant< 0 >();
        assert(tmp0 == tmp0);
      return [ tmp0 ] ( const DomainValueType &phi ) {
          return RangeValueType{ { tmp0 * (std::get< 0 >( phi ))[ 0 ] }, { { ((std::get< 1 >( phi ))[ 0 ])[ 0 ], ((std::get< 1 >( phi ))[ 0 ])[ 1 ] } } };
        };
    }

    template< class Point >
    RangeValueType boundary ( const Point &x, const DomainValueType &u ) const
    {
      using std::pow;
      GlobalCoordinateType tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
        assert(tmp0 == tmp0);
      const auto tmp1 = std::pow( tmp0[ 0 ], 2 );
        assert(tmp1 == tmp1);
      const auto tmp2 = 3 * tmp1;
        assert(tmp2 == tmp2);
      const auto tmp3 = 0.3333333333333333 * tmp2;
        assert(tmp3 == tmp3);
      const auto tmp4 = -1 * tmp3;
        assert(tmp4 == tmp4);
      const auto tmp5 = 2 * tmp0[ 0 ];
        assert(tmp5 == tmp5);
      const auto tmp6 = 0.5 * tmp5;
        assert(tmp6 == tmp6);
      const auto tmp7 = tmp6 + tmp4;
        assert(tmp7 == tmp7);
      GlobalCoordinateType tmp8 = intersection_.unitOuterNormal( x.localPosition() );
        assert(tmp8 == tmp8);
      const auto tmp9 = tmp8[ 0 ] * tmp7;
        assert(tmp9 == tmp9);
      const auto tmp10 = std::pow( tmp0[ 1 ], 2 );
        assert(tmp10 == tmp10);
      const auto tmp11 = 3 * tmp10;
        assert(tmp11 == tmp11);
      const auto tmp12 = -1 * tmp11;
        assert(tmp12 == tmp12);
      const auto tmp13 = 0.3333333333333333 * tmp12;
        assert(tmp13 == tmp13);
      const auto tmp14 = -1 * tmp13;
        assert(tmp14 == tmp14);
      const auto tmp15 = 2 * tmp0[ 1 ];
        assert(tmp15 == tmp15);
      const auto tmp16 = 0.5 * tmp15;
        assert(tmp16 == tmp16);
      const auto tmp17 = tmp16 + tmp14;
        assert(tmp17 == tmp17);
      const auto tmp18 = tmp8[ 1 ] * tmp17;
        assert(tmp18 == tmp18);
      const auto tmp19 = tmp18 + tmp9;
        assert(tmp19 == tmp19);
      const auto tmp20 = tmp19 * (tmp0[ 0 ] >= 1e-08 ? 1 : 0.0);
        assert(tmp20 == tmp20);
      const auto tmp21 = -1 * tmp20;
        assert(tmp21 == tmp21);
      return RangeValueType{ { tmp21 }, { { 0, 0 } } };
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
              domainId = tmp0[ 0 ] <= 1e-08;
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
          using std::pow;
          auto tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
            assert(tmp0 == tmp0);
          const auto tmp1 = std::pow( tmp0[ 0 ], 3 );
            assert(tmp1 == tmp1);
          const auto tmp2 = std::pow( tmp0[ 1 ], 3 );
            assert(tmp2 == tmp2);
          const auto tmp3 = -1 * tmp2;
            assert(tmp3 == tmp3);
          const auto tmp4 = tmp3 + tmp1;
            assert(tmp4 == tmp4);
          const auto tmp5 = 0.3333333333333333 * tmp4;
            assert(tmp5 == tmp5);
          const auto tmp6 = -1 * tmp5;
            assert(tmp6 == tmp6);
          const auto tmp7 = std::pow( tmp0[ 1 ], 2 );
            assert(tmp7 == tmp7);
          const auto tmp8 = std::pow( tmp0[ 0 ], 2 );
            assert(tmp8 == tmp8);
          const auto tmp9 = tmp8 + tmp7;
            assert(tmp9 == tmp9);
          const auto tmp10 = 0.5 * tmp9;
            assert(tmp10 == tmp10);
          const auto tmp11 = tmp10 + tmp6;
            assert(tmp11 == tmp11);
          const auto tmp12 = 1 + tmp11;
            assert(tmp12 == tmp12);
          result[ 0 ] = tmp12;
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
          using std::pow;
          auto tmp0 = geometry().global( Dune::Fem::coordinate( x ) );
            assert(tmp0 == tmp0);
          const auto tmp1 = std::pow( tmp0[ 0 ], 2 );
            assert(tmp1 == tmp1);
          const auto tmp2 = 3 * tmp1;
            assert(tmp2 == tmp2);
          const auto tmp3 = 0.3333333333333333 * tmp2;
            assert(tmp3 == tmp3);
          const auto tmp4 = -1 * tmp3;
            assert(tmp4 == tmp4);
          const auto tmp5 = 2 * tmp0[ 0 ];
            assert(tmp5 == tmp5);
          const auto tmp6 = 0.5 * tmp5;
            assert(tmp6 == tmp6);
          const auto tmp7 = tmp6 + tmp4;
            assert(tmp7 == tmp7);
          const auto tmp8 = std::pow( tmp0[ 1 ], 2 );
            assert(tmp8 == tmp8);
          const auto tmp9 = 3 * tmp8;
            assert(tmp9 == tmp9);
          const auto tmp10 = -1 * tmp9;
            assert(tmp10 == tmp10);
          const auto tmp11 = 0.3333333333333333 * tmp10;
            assert(tmp11 == tmp11);
          const auto tmp12 = -1 * tmp11;
            assert(tmp12 == tmp12);
          const auto tmp13 = 2 * tmp0[ 1 ];
            assert(tmp13 == tmp13);
          const auto tmp14 = 0.5 * tmp13;
            assert(tmp14 == tmp14);
          const auto tmp15 = tmp14 + tmp12;
            assert(tmp15 == tmp15);
          (result[ 0 ])[ 0 ] = tmp7;
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

    const Conmass &conmass () const
    {
      return *std::get< 0 >( constants_ );
    }

    Conmass &conmass ()
    {
      return *std::get< 0 >( constants_ );
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

} // namespace Integrands_04eff3707585c91909495d3077e951d5v1_4

PYBIND11_MODULE( integrands_04eff3707585c91909495d3077e951d5v1_4_19659fe2f2f74bcfb5d532f1a0d38be9, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView > GridPart;
  typedef Integrands_04eff3707585c91909495d3077e951d5v1_4::Integrands< GridPart > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_04eff3707585c91909495d3077e951d5v1_4::Integrands< GridPart >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_04eff3707585c91909495d3077e951d5v1_4_19659fe2f2f74bcfb5d532f1a0d38be9.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] () { return new Integrands(); } ) );
      cls.def_property( "mass", [] ( Integrands &self ) -> Integrands::Conmass { return self.conmass(); }, [] ( Integrands &self, const Integrands::Conmass &v ) { self.conmass() = v; } );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return true;});
  }
}
#endif // GuardIntegrands_04eff3707585c91909495d3077e951d5v1_4
