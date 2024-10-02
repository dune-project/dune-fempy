import hashlib

import sys
import logging
logger = logging.getLogger(__name__)

from dune.generator import Constructor, Method
from dune.generator.generator import SimpleGenerator
from dune.generator.generator import SimpleGenerator
from dune.fem.operator import load
generator = SimpleGenerator("Operator", "Dune::FemPy")

# limiter can be ScalingLimiter or FV based limiter with FV type reconstructions for troubled cells
def createLimiter(domainSpace, rangeSpace=None,
                  bounds = [1e-12,1.], limiter='scaling'):
    """
    Parameters:
        domainSpace  discrete space the domain function belongs to
        rangeSpace   discrete space the range function belongs to (default is None which means same as domainSpace)
        bounds       list of lists containing lower and upper bounds for each component. If a entry is none the component will not be limited.
        limiter      type of limiter, i.e. 'fv' or 'scaling'.

    Returns:
        Limiter object as fem operator.
    """

    if rangeSpace is None:
        rangeSpace = domainSpace

    domainSpaceType = domainSpace.cppTypeName
    rangeSpaceType = rangeSpace.cppTypeName

    # check for old parameters
    if len(bounds) == 2:
        foundList = False
        for b in bounds:
            if isinstance(b, (list,tuple)):
                foundList = True
        if not foundList:
            print(f"Deprecated parameters {bounds}: Use list of lists selecting bounds for each component, None if not limited!")
            # old style one list and last component is limited (see twophaseflow in tutorial)
            newBounds = [None]*domainSpace.dimRange
            newBounds[-1] = bounds
            bounds = newBounds
    else:
        raise Exception("createLimiter: bounds needs to be an instance of either list or tuple.")

    # obtain limited component numbers and corresponding bounds
    components = [ i for i, comp in enumerate(bounds) if comp is not None]
    newBounds = [(-1e308, 1e308) if comp is None else comp for comp in bounds]

    _, domainFunctionIncludes, domainFunctionType, _, _, _ = domainSpace.storage
    _, rangeFunctionIncludes, rangeFunctionType, _, _, _ = rangeSpace.storage

    includes = ["dune/fem-dg/operator/limiter/limiter.hh"]
    includes += domainSpace.cppIncludes + domainFunctionIncludes
    includes += rangeSpace.cppIncludes + rangeFunctionIncludes

    typeName = 'Dune::Fem::ScalingLimiter< ' + domainFunctionType + ', ' + rangeFunctionType + ' >'
    # FV type limiter where FV based reconstructions are done
    if limiter == 'fv':
        typeName = 'Dune::Fem::Limiter< ' + domainFunctionType + ', ' + rangeFunctionType + ' >'

    constructor = Constructor(['const '+domainSpaceType + ' &dSpace, const '+rangeSpaceType + ' &rSpace, const std::vector<int>& components, const std::vector< std::vector<double> >& bounds'],
                              ['return new ' + typeName + '(dSpace, rSpace, components, bounds);'],
                              ['"dSpace"_a', '"rSpace"_a', '"components"_a', '"bounds"_a',
                               'pybind11::keep_alive< 1, 2 >()', 'pybind11::keep_alive< 1, 3 >()'])

    # add method activated to inspect limited cells.
    activated = Method('activated', '&'+typeName+'::activated')

    return load(includes, typeName, constructor, activated).Operator( domainSpace, rangeSpace, components, newBounds )


def createOrderRedcution(domainSpace):

    domainSpaceType = domainSpace.cppTypeName

    _, domainFunctionIncludes, domainFunctionType, _, _, _ = domainSpace.storage

    includes = ["dune/fem-dg/operator/common/orderreduction.hh"]
    includes += domainSpace.cppIncludes + domainFunctionIncludes

    typeName = 'Dune::Fem::OrderReduction< ' + domainFunctionType + ' >'

    constructor = Constructor(['const '+domainSpaceType + ' &dSpace'],
                              ['return new ' + typeName + '(dSpace);'],
                              ['"dSpace"_a','pybind11::keep_alive< 1, 2 >()'])

    # add method maxRelevantOrder to get max order that is relevant per cell
    # maxRelevantOrder = Method('maxRelevantOrder', '&'+typeName+'::maxRelevantOrder')

    # return load(includes, typeName, constructor, maxRelevantOrder).Operator( domainSpace )
    return load(includes, typeName, constructor).Operator( domainSpace )
