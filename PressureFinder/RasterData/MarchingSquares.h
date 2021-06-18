#ifndef MARCHING_SQUARES_H
#define MARCHING_SQUARES_H

#include <unordered_set>
#include <array>

#include "../Math/Vector2.h"
#include "../Graphics/2d/PolyLine.h"

#include "./Image2d.h"

#ifdef MYMATH_NAMESPACE 
#	define MY_MATH_NS MYMATH_NAMESPACE 
#else
#	define MY_MATH_NS 
#endif

class MarchingSquares 
{
public:
	MarchingSquares();
	~MarchingSquares();
	
	std::vector<d2::PolyLine> & GetAllContours();
	const std::vector<d2::PolyLine> & GetAllContours() const;
		
	void SetBorderValue(float value);

	template <typename T>
	void Run(const Image2d<T> & data, float threshold);

	template <typename T>
	void Run(const Image2d<T> & data, int x, int y, int w, int h, float threshold);

protected:

	struct MarchingEdge
	{
		MY_MATH_NS::Vector2d pt;
		uint64_t edgeKey;
	};

	// 
	//  0 --- 1
	//  |     |
	//  |     |
	//  3 --- 2
	//
	// Edges: 0 = 01, 1 = 12, 2 = 23, 3 = 30
	//

	//indices of vertices for edges
	static const std::array<std::array< int, 2>, 4> edgesVertices;
	
	//offset to next cell over edge[id]
	//{dx, dy}
	static const std::array<std::array< int, 2>, 4> edgesNextCell;
	
	//list of edges for each case of corner values
	//[corners] = {edgesIds... -1} (-1 = no edge)
	//there can be max 4 edges for a cell
	static const std::array<std::array<int, 4>, 16> edges;
	
	std::unordered_set<uint64_t> visitedEdges;
	std::vector<d2::PolyLine> linesAll;

	float borderValue;

	inline uint64_t CalculateEdgeKey(uint32_t v0, uint32_t v1);

	inline int GetIndexFromPosition(int x, int y, int w) const;

	template <typename T>
	float GetValue(int x, int y, const Image2d<T> & data) const;

	template <typename T>
	void TrackContour(int x, int y, const Image2d<T> & data, float threshold);
};


#endif