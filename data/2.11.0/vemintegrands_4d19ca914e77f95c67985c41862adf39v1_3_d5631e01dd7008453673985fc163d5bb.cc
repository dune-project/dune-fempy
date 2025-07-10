#ifndef GuardIntegrands_4d19ca914e77f95c67985c41862adf39v1_3
#define GuardIntegrands_4d19ca914e77f95c67985c41862adf39v1_3
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
#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/function/simplegridfunction.hh>
#include <dune/fem/misc/gridfunctionview.hh>
#include <dune/vem/py/integrands.hh>

namespace Integrands_4d19ca914e77f95c67985c41862adf39v1_3
{




  // Integrands
  // ----------

  template< class GridPart, class CoeffoldSolution >
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
    typedef double Coneps;
    typedef std::tuple< std::shared_ptr< Condt >, std::shared_ptr< Coneps > > ConstantTupleType;
    template< std::size_t i >
    using ConstantsRangeType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;
    typedef std::tuple< Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > > > CoefficientFunctionSpaceTupleType;
    typedef std::tuple< CoeffoldSolution > CoefficientTupleType;
    template< std::size_t i >
    using CoefficientFunctionSpaceType = std::tuple_element_t< i, CoefficientFunctionSpaceTupleType >;
    template< std::size_t i >
    using CoefficientRangeType = typename CoefficientFunctionSpaceType< i >::RangeType;
    template< std::size_t i >
    using CoefficientJacobianRangeType = typename CoefficientFunctionSpaceType< i >::JacobianRangeType;
    static constexpr bool gridPartValid = Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<CoeffoldSolution>>();
    template< std::size_t i >
    using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
    template< std::size_t i >
    using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;

    Integrands ( const CoeffoldSolution &coeffoldSolution, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
      : coefficients_( Dune::Fem::ConstLocalFunction< CoeffoldSolution >( coeffoldSolution ) )
    {
      std::get< 0 >( constants_ ) = std::make_shared< Condt >( (Condt(0)) );
      std::get< 1 >( constants_ ) = std::make_shared< Coneps >( (Coneps(0)) );
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
    typedef std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 >, Dune::Fem::ExplicitFieldVector< Dune::FieldMatrix< double, 2, 2 >, 1 > > DomainValueType;
    typedef std::tuple< Dune::FieldVector< double, 1 >, Dune::FieldMatrix< double, 1, 2 >, Dune::Fem::ExplicitFieldVector< Dune::FieldMatrix< double, 2, 2 >, 1 > > RangeValueType;
    static constexpr bool _nonlinear = true;

    bool nonlinear () const
    {
      return _nonlinear;
    }

    template< class Point >
    RangeValueType interior ( const Point &x, const DomainValueType &u ) const
    {
      using std::pow;
      typename CoefficientFunctionSpaceType< 0 >::RangeType tmp0 = evaluateCoefficient< 0 >( x );
      const auto tmp1 = -1 * tmp0[ 0 ];
      const auto tmp2 = std::get< 0 >( u );
      const auto tmp3 = tmp2[ 0 ] + tmp1;
      const auto tmp4 = std::get< 1 >( u );
      const auto tmp5 = 2 * (tmp4[ 0 ])[ 0 ];
      const auto tmp6 = tmp2[ 0 ] * tmp5;
      const auto tmp7 = tmp2[ 0 ] * tmp6;
      const auto tmp8 = std::pow( tmp2[ 0 ], 2 );
      const auto tmp9 = -1 + tmp8;
      const auto tmp10 = (tmp4[ 0 ])[ 0 ] * tmp9;
      const auto tmp11 = tmp10 + tmp7;
      double tmp12 = constant< 0 >();
      const auto tmp13 = tmp12 * tmp11;
      const auto tmp14 = 2 * (tmp4[ 0 ])[ 1 ];
      const auto tmp15 = tmp2[ 0 ] * tmp14;
      const auto tmp16 = tmp2[ 0 ] * tmp15;
      const auto tmp17 = (tmp4[ 0 ])[ 1 ] * tmp9;
      const auto tmp18 = tmp17 + tmp16;
      const auto tmp19 = tmp12 * tmp18;
      double tmp20 = constant< 1 >();
      const auto tmp21 = tmp12 * tmp20;
      const auto tmp22 = tmp20 * tmp21;
      const auto tmp23 = std::get< 2 >( u );
      const auto tmp24 = ((tmp23[ 0 ])[ 0 ])[ 0 ] * tmp22;
      const auto tmp25 = ((tmp23[ 0 ])[ 0 ])[ 1 ] * tmp22;
      const auto tmp26 = ((tmp23[ 0 ])[ 1 ])[ 0 ] * tmp22;
      const auto tmp27 = ((tmp23[ 0 ])[ 1 ])[ 1 ] * tmp22;
      return RangeValueType{ { tmp3 }, { { tmp13, tmp19 } }, { { { tmp24, tmp25 }, { tmp26, tmp27 } } } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      using std::pow;
      const auto tmp0 = std::get< 0 >( u );
      const auto tmp1 = 2 * tmp0[ 0 ];
      const auto tmp2 = std::get< 1 >( u );
      const auto tmp3 = (tmp2[ 0 ])[ 0 ] * tmp1;
      const auto tmp4 = 2 * (tmp2[ 0 ])[ 0 ];
      const auto tmp5 = tmp0[ 0 ] * tmp4;
      const auto tmp6 = tmp5 + tmp5;
      const auto tmp7 = tmp6 + tmp3;
      double tmp8 = constant< 0 >();
      const auto tmp9 = tmp8 * tmp7;
      const auto tmp10 = (tmp2[ 0 ])[ 1 ] * tmp1;
      const auto tmp11 = 2 * (tmp2[ 0 ])[ 1 ];
      const auto tmp12 = tmp0[ 0 ] * tmp11;
      const auto tmp13 = tmp12 + tmp12;
      const auto tmp14 = tmp13 + tmp10;
      const auto tmp15 = tmp8 * tmp14;
      const auto tmp16 = tmp0[ 0 ] * tmp1;
      const auto tmp17 = std::pow( tmp0[ 0 ], 2 );
      const auto tmp18 = -1 + tmp17;
      const auto tmp19 = tmp18 + tmp16;
      const auto tmp20 = tmp8 * tmp19;
      double tmp21 = constant< 1 >();
      const auto tmp22 = tmp8 * tmp21;
      const auto tmp23 = tmp21 * tmp22;
      return [ tmp15, tmp20, tmp23, tmp9 ] ( const DomainValueType &phi ) {
          return RangeValueType{ { (std::get< 0 >( phi ))[ 0 ] }, { { tmp9 * (std::get< 0 >( phi ))[ 0 ] + tmp20 * ((std::get< 1 >( phi ))[ 0 ])[ 0 ], tmp15 * (std::get< 0 >( phi ))[ 0 ] + tmp20 * ((std::get< 1 >( phi ))[ 0 ])[ 1 ] } }, { { { tmp23 * (((std::get< 2 >( phi ))[ 0 ])[ 0 ])[ 0 ], tmp23 * (((std::get< 2 >( phi ))[ 0 ])[ 0 ])[ 1 ] }, { tmp23 * (((std::get< 2 >( phi ))[ 0 ])[ 1 ])[ 0 ], tmp23 * (((std::get< 2 >( phi ))[ 0 ])[ 1 ])[ 1 ] } } } };
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

    const Condt &condt () const
    {
      return *std::get< 0 >( constants_ );
    }

    Condt &condt ()
    {
      return *std::get< 0 >( constants_ );
    }

    const Coneps &coneps () const
    {
      return *std::get< 1 >( constants_ );
    }

    Coneps &coneps ()
    {
      return *std::get< 1 >( constants_ );
    }

    const EntityType &entity () const
    {
      return entity_;
    }

    template< std::size_t i, class Point >
    typename CoefficientFunctionSpaceType< i >::RangeType evaluateCoefficient ( const Point &x ) const
    {
      typename CoefficientFunctionSpaceType< i >::RangeType result;
      std::get< i >( coefficients_ ).evaluate( x, result );;
      return result;
    }

    template< std::size_t i, class Point >
    typename CoefficientFunctionSpaceType< i >::JacobianRangeType jacobianCoefficient ( const Point &x ) const
    {
      typename CoefficientFunctionSpaceType< i >::JacobianRangeType result;
      std::get< i >( coefficients_ ).jacobian( x, result );;
      return result;
    }

    template< std::size_t i, class Point >
    typename CoefficientFunctionSpaceType< i >::HessianRangeType hessianCoefficient ( const Point &x ) const
    {
      typename CoefficientFunctionSpaceType< i >::HessianRangeType result;
      std::get< i >( coefficients_ ).hessian( x, result );;
      return result;
    }
    EntityType entity_;
    IntersectionType intersection_;
    ConstantTupleType constants_;
    std::tuple< Dune::Fem::ConstLocalFunction< CoeffoldSolution > > coefficients_;

  public:
    template< class Point, class DRangeType >
    RRangeType hessStabilization ( const Point &x, const DRangeType &u ) const
    {
      RRangeType result;
      double tmp0 = constant< 1 >();
      const auto tmp1 = tmp0 * tmp0;
      double tmp2 = constant< 0 >();
      const auto tmp3 = tmp2 * tmp1;
      result[ 0 ] = tmp3;
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
      double tmp0 = constant< 0 >();
      const auto tmp1 = 2 * tmp0;
      result[ 0 ] = tmp1;
      return result;
    }

    template< class Point, class DRangeType >
    RRangeType linGradStabilization ( const Point &x, const DRangeType &u ) const
    {
      RRangeType result;
      result[ 0 ] = 0.0;
      return result;
    }

    template< class Point, class DRangeType >
    RRangeType massStabilization ( const Point &x, const DRangeType &u ) const
    {
      RRangeType result;
      result[ 0 ] = 1;
      return result;
    }

    template< class Point, class DRangeType >
    RRangeType linMassStabilization ( const Point &x, const DRangeType &u ) const
    {
      RRangeType result;
      result[ 0 ] = 0.0;
      return result;
    }
  };

} // namespace Integrands_4d19ca914e77f95c67985c41862adf39v1_3

PYBIND11_MODULE( vemintegrands_4d19ca914e77f95c67985c41862adf39v1_3_d5631e01dd7008453673985fc163d5bb, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::Vem::CubeGrid<2,2>::LeafGridView > GridPart;
  typedef Integrands_4d19ca914e77f95c67985c41862adf39v1_3::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Vem::CubeGrid<2,2>::LeafGridView>, Dune::FieldVector< double, 1 > > > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_4d19ca914e77f95c67985c41862adf39v1_3::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Vem::CubeGrid<2,2>::LeafGridView>, Dune::FieldVector< double, 1 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/vemintegrands_4d19ca914e77f95c67985c41862adf39v1_3_d5631e01dd7008453673985fc163d5bb.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] ( const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::Vem::CubeGrid<2,2>::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffoldSolution ) { return new Integrands( coeffoldSolution ); } ), pybind11::keep_alive< 1, 2 >() );
      cls.def_property( "dt", [] ( Integrands &self ) -> Integrands::Condt { return self.condt(); }, [] ( Integrands &self, const Integrands::Condt &v ) { self.condt() = v; } );
      cls.def_property( "eps", [] ( Integrands &self ) -> Integrands::Coneps { return self.coneps(); }, [] ( Integrands &self, const Integrands::Coneps &v ) { self.coneps() = v; } );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return true;});
  }
}
#endif // GuardIntegrands_4d19ca914e77f95c67985c41862adf39v1_3
