

set(
		libqhullcpp_HEADERS
		Coordinates.h
		functionObjects.h
		PointCoordinates.h
		Qhull.h
		QhullError.h
		QhullFacet.h
		QhullFacetList.h
		QhullFacetSet.h
		QhullHyperplane.h
		QhullIterator.h
		QhullLinkedList.h
		QhullPoint.h
		QhullPoints.h
		QhullPointSet.h
		QhullQh.h
		QhullRidge.h
		QhullSet.h
		QhullSets.h
		QhullStat.h
		QhullUser.h
		QhullVertex.h
		QhullVertexSet.h
		RboxPoints.h
		RoadError.h
		RoadLogEvent.h
		#src/qhulltest/RoadTest.h
)

set(
		libqhullcpp_SOURCES
		Coordinates.cpp
		PointCoordinates.cpp
		Qhull.cpp
		QhullFacet.cpp
		QhullFacetList.cpp
		QhullFacetSet.cpp
		QhullHyperplane.cpp
		QhullPoint.cpp
		QhullPointSet.cpp
		QhullPoints.cpp
		QhullQh.cpp
		QhullRidge.cpp
		QhullSet.cpp
		QhullStat.cpp
		QhullUser.cpp
		QhullVertex.cpp
		QhullVertexSet.cpp
		RboxPoints.cpp
		RoadError.cpp
		RoadLogEvent.cpp
		${libqhullcpp_HEADERS}
)



add_library(libqhullcpp ${libqhullcpp_SOURCES})