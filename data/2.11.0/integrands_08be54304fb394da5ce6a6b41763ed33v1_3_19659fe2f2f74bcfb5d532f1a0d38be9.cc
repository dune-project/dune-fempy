#ifndef GuardIntegrands_08be54304fb394da5ce6a6b41763ed33v1_3
#define GuardIntegrands_08be54304fb394da5ce6a6b41763ed33v1_3
#define USING_DUNE_PYTHON 1
#include <config.h>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/python/grid/hierarchical.hh>
#include <tuple>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/extensions.h>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/geometry/edgelength.hh>
#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/function/simplegridfunction.hh>
#include <dune/fem/misc/gridfunctionview.hh>
#include <dune/fempy/py/integrands.hh>

namespace Integrands_08be54304fb394da5ce6a6b41763ed33v1_3
{




  // Integrands
  // ----------

  template< class GridPart, class Coeffanalytical_solution_x, class Coeffanalytical_solution_y, class Coeffanalytical_solution_xA, class Coeffanalytical_solution_yA, class Coeffanalytical_solution_xAA, class Coeffanalytical_solution_yAA, class Coeffanalytical_solution_xAAA, class Coeffanalytical_solution_yAAA >
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
    typedef std::tuple< Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > >, Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > >, Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > >, Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > >, Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > >, Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > >, Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > >, Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< double, 1 > > > CoefficientFunctionSpaceTupleType;
    typedef std::tuple< Coeffanalytical_solution_x, Coeffanalytical_solution_y, Coeffanalytical_solution_xA, Coeffanalytical_solution_yA, Coeffanalytical_solution_xAA, Coeffanalytical_solution_yAA, Coeffanalytical_solution_xAAA, Coeffanalytical_solution_yAAA > CoefficientTupleType;
    template< std::size_t i >
    using CoefficientFunctionSpaceType = std::tuple_element_t< i, CoefficientFunctionSpaceTupleType >;
    template< std::size_t i >
    using CoefficientRangeType = typename CoefficientFunctionSpaceType< i >::RangeType;
    template< std::size_t i >
    using CoefficientJacobianRangeType = typename CoefficientFunctionSpaceType< i >::JacobianRangeType;
    static constexpr bool gridPartValid = Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffanalytical_solution_x>>() && Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffanalytical_solution_y>>() && Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffanalytical_solution_xA>>() && Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffanalytical_solution_yA>>() && Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffanalytical_solution_xAA>>() && Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffanalytical_solution_yAA>>() && Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffanalytical_solution_xAAA>>() && Dune::Fem::checkGridPartValid<GridPartType,Dune::Fem::ConstLocalFunction<Coeffanalytical_solution_yAAA>>();
    template< std::size_t i >
    using CoefficientType = std::tuple_element_t< i, CoefficientTupleType >;
    template< std::size_t i >
    using ConstantType = typename std::tuple_element_t< i, ConstantTupleType >::element_type;

    Integrands ( const Coeffanalytical_solution_x &coeffanalytical_solution_x, const Coeffanalytical_solution_y &coeffanalytical_solution_y, const Coeffanalytical_solution_xA &coeffanalytical_solution_xA, const Coeffanalytical_solution_yA &coeffanalytical_solution_yA, const Coeffanalytical_solution_xAA &coeffanalytical_solution_xAA, const Coeffanalytical_solution_yAA &coeffanalytical_solution_yAA, const Coeffanalytical_solution_xAAA &coeffanalytical_solution_xAAA, const Coeffanalytical_solution_yAAA &coeffanalytical_solution_yAAA, const Dune::Fem::ParameterReader &parameter = Dune::Fem::Parameter::container() )
      : coefficients_( Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_x >( coeffanalytical_solution_x ), Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_y >( coeffanalytical_solution_y ), Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_xA >( coeffanalytical_solution_xA ), Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_yA >( coeffanalytical_solution_yA ), Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_xAA >( coeffanalytical_solution_xAA ), Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_yAA >( coeffanalytical_solution_yAA ), Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_xAAA >( coeffanalytical_solution_xAAA ), Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_yAAA >( coeffanalytical_solution_yAAA ) )
    {}

    bool init ( const EntityType &entity )
    {
      entity_ = entity;
      std::get< 0 >( coefficients_ ).bind( this->entity() );
      std::get< 1 >( coefficients_ ).bind( this->entity() );
      std::get< 2 >( coefficients_ ).bind( this->entity() );
      std::get< 3 >( coefficients_ ).bind( this->entity() );
      std::get< 4 >( coefficients_ ).bind( this->entity() );
      std::get< 5 >( coefficients_ ).bind( this->entity() );
      std::get< 6 >( coefficients_ ).bind( this->entity() );
      std::get< 7 >( coefficients_ ).bind( this->entity() );
      return true;
    }

    void unbind ()
    {
      std::get< 0 >( coefficients_ ).unbind();
      std::get< 1 >( coefficients_ ).unbind();
      std::get< 2 >( coefficients_ ).unbind();
      std::get< 3 >( coefficients_ ).unbind();
      std::get< 4 >( coefficients_ ).unbind();
      std::get< 5 >( coefficients_ ).unbind();
      std::get< 6 >( coefficients_ ).unbind();
      std::get< 7 >( coefficients_ ).unbind();
    }

    bool init ( const IntersectionType &intersection )
    {
      intersection_ = intersection;
      return (intersection.boundary() && init( intersection.inside() ));
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
      const auto tmp0 = std::get< 0 >( u );
      const auto tmp1 = std::get< 1 >( u );
      return RangeValueType{ { 0.0, 0.0, tmp0[ 2 ] }, { { (tmp1[ 0 ])[ 0 ], (tmp1[ 0 ])[ 1 ] }, { (tmp1[ 1 ])[ 0 ], (tmp1[ 1 ])[ 1 ] }, { 0.0, 0.0 } } };
    }

    template< class Point >
    auto linearizedInterior ( const Point &x, const DomainValueType &u ) const
    {
      return [] ( const DomainValueType &phi ) {
          return RangeValueType{ { 0, 0, (std::get< 0 >( phi ))[ 2 ] }, { { ((std::get< 1 >( phi ))[ 0 ])[ 0 ], ((std::get< 1 >( phi ))[ 0 ])[ 1 ] }, { ((std::get< 1 >( phi ))[ 1 ])[ 0 ], ((std::get< 1 >( phi ))[ 1 ])[ 1 ] }, { 0, 0 } } };
        };
    }
    typedef Dune::FieldVector< double, 3 >  RRangeType;
    typedef Dune::FieldMatrix< double, 3, GridPartType::dimension >  RJacobianRangeType;
    typedef Dune::Fem::BoundaryIdProvider< typename GridPartType::GridType > BoundaryIdProviderType;
    typedef std::array<int,3> DirichletComponentType;

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
              const auto tmp0 = Dune::Fem::boundaryId<GridPartType>( (intersection_) );
              domainId = tmp0 == 1;
              if (domainId)
              {
                std::fill( dirichletComponent.begin(), dirichletComponent.end(), 2 );
                dirichletComponent[2] = 0;
                return true;
              }
          }
          {
              const auto tmp0 = Dune::Fem::boundaryId<GridPartType>( (intersection_) );
              domainId = tmp0 == 2;
              if (domainId)
              {
                std::fill( dirichletComponent.begin(), dirichletComponent.end(), 3 );
                dirichletComponent[2] = 0;
                return true;
              }
          }
          {
              const auto tmp0 = Dune::Fem::boundaryId<GridPartType>( (intersection_) );
              domainId = tmp0 == 3;
              if (domainId)
              {
                std::fill( dirichletComponent.begin(), dirichletComponent.end(), 4 );
                dirichletComponent[2] = 0;
                return true;
              }
          }
          {
              const auto tmp0 = Dune::Fem::boundaryId<GridPartType>( (intersection_) );
              domainId = tmp0 == 4;
              if (domainId)
              {
                std::fill( dirichletComponent.begin(), dirichletComponent.end(), 5 );
                dirichletComponent[2] = 0;
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
          typename CoefficientFunctionSpaceType< 0 >::RangeType tmp0 = evaluateCoefficient< 0 >( x );
          typename CoefficientFunctionSpaceType< 1 >::RangeType tmp1 = evaluateCoefficient< 1 >( x );
          result[ 0 ] = tmp0[ 0 ];
          result[ 1 ] = tmp1[ 0 ];
          result[ 2 ] = 0.0;
        }
        break;
      case 3:
        {
          typename CoefficientFunctionSpaceType< 2 >::RangeType tmp0 = evaluateCoefficient< 2 >( x );
          typename CoefficientFunctionSpaceType< 3 >::RangeType tmp1 = evaluateCoefficient< 3 >( x );
          result[ 0 ] = tmp0[ 0 ];
          result[ 1 ] = tmp1[ 0 ];
          result[ 2 ] = 0.0;
        }
        break;
      case 4:
        {
          typename CoefficientFunctionSpaceType< 4 >::RangeType tmp0 = evaluateCoefficient< 4 >( x );
          typename CoefficientFunctionSpaceType< 5 >::RangeType tmp1 = evaluateCoefficient< 5 >( x );
          result[ 0 ] = tmp0[ 0 ];
          result[ 1 ] = tmp1[ 0 ];
          result[ 2 ] = 0.0;
        }
        break;
      case 5:
        {
          typename CoefficientFunctionSpaceType< 6 >::RangeType tmp0 = evaluateCoefficient< 6 >( x );
          typename CoefficientFunctionSpaceType< 7 >::RangeType tmp1 = evaluateCoefficient< 7 >( x );
          result[ 0 ] = tmp0[ 0 ];
          result[ 1 ] = tmp1[ 0 ];
          result[ 2 ] = 0.0;
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
          typename CoefficientFunctionSpaceType< 0 >::JacobianRangeType tmp0 = jacobianCoefficient< 0 >( x );
          typename CoefficientFunctionSpaceType< 1 >::JacobianRangeType tmp1 = jacobianCoefficient< 1 >( x );
          (result[ 0 ])[ 0 ] = (tmp0[ 0 ])[ 0 ];
          (result[ 0 ])[ 1 ] = (tmp0[ 0 ])[ 1 ];
          (result[ 1 ])[ 0 ] = (tmp1[ 0 ])[ 0 ];
          (result[ 1 ])[ 1 ] = (tmp1[ 0 ])[ 1 ];
          (result[ 2 ])[ 0 ] = 0.0;
          (result[ 2 ])[ 1 ] = 0.0;
        }
        break;
      case 3:
        {
          typename CoefficientFunctionSpaceType< 2 >::JacobianRangeType tmp0 = jacobianCoefficient< 2 >( x );
          typename CoefficientFunctionSpaceType< 3 >::JacobianRangeType tmp1 = jacobianCoefficient< 3 >( x );
          (result[ 0 ])[ 0 ] = (tmp0[ 0 ])[ 0 ];
          (result[ 0 ])[ 1 ] = (tmp0[ 0 ])[ 1 ];
          (result[ 1 ])[ 0 ] = (tmp1[ 0 ])[ 0 ];
          (result[ 1 ])[ 1 ] = (tmp1[ 0 ])[ 1 ];
          (result[ 2 ])[ 0 ] = 0.0;
          (result[ 2 ])[ 1 ] = 0.0;
        }
        break;
      case 4:
        {
          typename CoefficientFunctionSpaceType< 4 >::JacobianRangeType tmp0 = jacobianCoefficient< 4 >( x );
          typename CoefficientFunctionSpaceType< 5 >::JacobianRangeType tmp1 = jacobianCoefficient< 5 >( x );
          (result[ 0 ])[ 0 ] = (tmp0[ 0 ])[ 0 ];
          (result[ 0 ])[ 1 ] = (tmp0[ 0 ])[ 1 ];
          (result[ 1 ])[ 0 ] = (tmp1[ 0 ])[ 0 ];
          (result[ 1 ])[ 1 ] = (tmp1[ 0 ])[ 1 ];
          (result[ 2 ])[ 0 ] = 0.0;
          (result[ 2 ])[ 1 ] = 0.0;
        }
        break;
      case 5:
        {
          typename CoefficientFunctionSpaceType< 6 >::JacobianRangeType tmp0 = jacobianCoefficient< 6 >( x );
          typename CoefficientFunctionSpaceType< 7 >::JacobianRangeType tmp1 = jacobianCoefficient< 7 >( x );
          (result[ 0 ])[ 0 ] = (tmp0[ 0 ])[ 0 ];
          (result[ 0 ])[ 1 ] = (tmp0[ 0 ])[ 1 ];
          (result[ 1 ])[ 0 ] = (tmp1[ 0 ])[ 0 ];
          (result[ 1 ])[ 1 ] = (tmp1[ 0 ])[ 1 ];
          (result[ 2 ])[ 0 ] = 0.0;
          (result[ 2 ])[ 1 ] = 0.0;
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
    std::tuple< Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_x >, Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_y >, Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_xA >, Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_yA >, Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_xAA >, Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_yAA >, Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_xAAA >, Dune::Fem::ConstLocalFunction< Coeffanalytical_solution_yAAA > > coefficients_;
  };

} // namespace Integrands_08be54304fb394da5ce6a6b41763ed33v1_3

PYBIND11_MODULE( integrands_08be54304fb394da5ce6a6b41763ed33v1_3_19659fe2f2f74bcfb5d532f1a0d38be9, module )
{
  typedef typename Dune::FemPy::GridPart< Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView > GridPart;
  typedef Integrands_08be54304fb394da5ce6a6b41763ed33v1_3::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > > Integrands;
  if constexpr( Integrands::gridPartValid )
  {
      auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("Integrands_08be54304fb394da5ce6a6b41763ed33v1_3::Integrands< GridPart, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > >, Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > >"), Dune::Python::IncludeFiles({"python/dune/generated/integrands_08be54304fb394da5ce6a6b41763ed33v1_3_19659fe2f2f74bcfb5d532f1a0d38be9.cc"})).first;
      Dune::FemPy::registerIntegrands< Integrands >( module, cls );
      cls.def( pybind11::init( [] ( const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffanalytical_solution_x, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffanalytical_solution_y, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffanalytical_solution_xA, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffanalytical_solution_yA, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffanalytical_solution_xAA, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffanalytical_solution_yAA, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffanalytical_solution_xAAA, const Dune::FemPy::VirtualizedGridFunction< Dune::FemPy::GridPart<Dune::YaspGrid< 2, Dune::EquidistantOffsetCoordinates< double, 2 > >::LeafGridView>, Dune::FieldVector< double, 1 > > &coeffanalytical_solution_yAAA ) { return new Integrands( coeffanalytical_solution_x, coeffanalytical_solution_y, coeffanalytical_solution_xA, coeffanalytical_solution_yA, coeffanalytical_solution_xAA, coeffanalytical_solution_yAA, coeffanalytical_solution_xAAA, coeffanalytical_solution_yAAA ); } ), pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 4 >(), pybind11::keep_alive< 1, 5 >(), pybind11::keep_alive< 1, 6 >(), pybind11::keep_alive< 1, 7 >(), pybind11::keep_alive< 1, 8 >(), pybind11::keep_alive< 1, 9 >() );
      cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return true;});
      cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return true;});
  }
}
#endif // GuardIntegrands_08be54304fb394da5ce6a6b41763ed33v1_3
