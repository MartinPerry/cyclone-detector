#include "./MarchingSquares.h"

#include <vector>
#include <algorithm>

#include "../Math/MathUtils.h"

#ifdef MYMATH_NAMESPACE
using namespace MYMATH_NAMESPACE;
#endif

// 
//  0 --- 1
//  |     |
//  |     |
//  3 --- 2
//
// Edges: 0 = 01, 1 = 12, 2 = 23, 3 = 30
//

//indices of vertices for edges
const std::array<std::array<int, 2>, 4> MarchingSquares::edgesVertices = 
{
	std::array<int, 2>{0, 1},
	std::array<int, 2>{1, 2},
	std::array<int, 2>{2, 3},
	std::array<int, 2>{3, 0}
};

//offset to next cell over edge[id]
//{dx, dy}
const std::array<std::array< int, 2>, 4> MarchingSquares::edgesNextCell = 
{
	std::array<int, 2>{0, -1},
	std::array<int, 2>{1, 0},
	std::array<int, 2>{0, 1},
	std::array<int, 2>{-1, 0}
};

//list of edges for each case of corner values
//[corners] = {edgesIds... -1} (-1 = no edge)
//there can be max 4 edges for a cell
const std::array<std::array<int, 4>, 16> MarchingSquares::edges = 
[]() ->std::array<std::array<int, 4>, 16> 
{

	std::array<std::array<int, 4>, 16> edges;
	edges[0b0000] = { -1, -1, -1, -1 };

	edges[0b0001] = { 0, 3, -1, -1 };
	edges[0b0010] = { 0, 1, -1, -1 };
	edges[0b0100] = { 1, 2, -1, -1 };
	edges[0b1000] = { 2, 3, -1, -1 };

	edges[0b0011] = { 1, 3, -1, -1 };
	edges[0b0101] = { 0, 1, 2, 3 };
	edges[0b1001] = { 0, 2, -1, -1 };
	edges[0b0110] = { 0, 2, -1, -1 };
	edges[0b1010] = { 0, 1, 2, 3 };
	edges[0b1100] = { 1, 3, -1, -1 };

	edges[0b1110] = { 0, 3, -1, -1 };
	edges[0b1101] = { 0, 1, -1, -1 };
	edges[0b1011] = { 1, 2, -1, -1 };
	edges[0b0111] = { 2, 3, -1, -1 };

	edges[0b1111] = { -1, -1, -1, -1 };

	return edges;
}();


MarchingSquares::MarchingSquares() : 
	borderValue(static_cast<float>(std::numeric_limits<float>::max()))
{
}

MarchingSquares::~MarchingSquares()
{
}

void MarchingSquares::SetBorderValue(float value)
{
	this->borderValue = value;
}

std::vector<d2::PolyLine> & MarchingSquares::GetAllContours()
{
	return this->linesAll;
}

const std::vector<d2::PolyLine> & MarchingSquares::GetAllContours() const
{
	return this->linesAll;
}

/// <summary>
/// Calculate key for edge
/// This key is used to prevent multiple processing of the same edges
/// </summary>
/// <param name="v0"></param>
/// <param name="v1"></param>
/// <returns></returns>
uint64_t MarchingSquares::CalculateEdgeKey(uint32_t v0, uint32_t v1)
{
	if (v0 > v1) std::swap(v0, v1);

	return ((static_cast<uint64_t>(v0) << 32) | static_cast<uint64_t>(v1));
}

/// <summary>
/// Create vitual index for enlarged image
/// - Add one row / column from each side
/// => minimal top left is [-1, -1] and we "move it" to [0, 0]
/// => we also need to increase width by w + 2 (one column from each side)
/// </summary>
/// <param name="x"></param>
/// <param name="y"></param>
/// <param name="w"></param>
/// <returns></returns>
int MarchingSquares::GetIndexFromPosition(int x, int y, int w) const
{	
	return (x + 1) + (y + 1) * (w + 2);
}

template <typename T>
void MarchingSquares::Run(const Image2d<T> & data, float threshold)
{
	for (int y = 0; y < data.GetHeight() - 1; y++)
	{
		for (int x = 0; x < data.GetWidth() - 1; x++)
		{
			this->TrackContour(x, y, data, threshold);
		}
	}
}

/// <summary>
/// Runf for sub-image  from [x, y] to [x + w, y + h]
/// </summary>
/// <param name="data"></param>
/// <param name="xx"></param>
/// <param name="yy"></param>
/// <param name="w"></param>
/// <param name="h"></param>
/// <param name="threshold"></param>
template <typename T>
void MarchingSquares::Run(const Image2d<T> & data, int xx, int yy, int w, int h, float threshold)
{
	for (int y = yy; y < yy + h; y++)
	{
		for (int x = xx; x < xx + w; x++)
		{
			this->TrackContour(x, y, data, threshold);
		}
	}
}


/// <summary>
/// Get value from data
/// If input index is outside image data, return minimal value of T
/// to enclose contour
/// </summary>
/// <param name="index"></param>
/// <param name="data"></param>
/// <returns></returns>
template <typename T>
float MarchingSquares::GetValue(int x, int y, const Image2d<T> & data) const
{			
	if ((x < 0) || 
		(y < 0) || 
		(x >= data.GetWidth()) ||
		(y >= data.GetHeight()))
	{
		return this->borderValue;
	}	

	return static_cast<float>(*data.GetPixelStart(x, y));
}

template <typename T>
void MarchingSquares::TrackContour(int x, int y, const Image2d<T> & data, float t)
{
	int xx = x;
	int yy = y;
	bool canContinue = true;

	d2::PolyLine pl;
	float v[4]; 
	std::array<int, 3> indexes[4];

	while (canContinue)
	{		
		canContinue = false;
		
		/*
		//nejde, proc
		//get pixel values		
		if ((x >= 0) &&
			(y >= 0) &&
			(x < data.GetWidth() - 1) &&
			(y < data.GetHeight() - 1))
		{
			//if we are inside picture - no need for border conditions
			//also we can read neighbor pixels together (assuming that data are one channel)

			const T * tmp0 = data.GetPixelStart(xx, yy);
			const T * tmp1 = data.GetPixelStart(xx, yy + 1);
			v[0] = static_cast<float>(tmp0[0]);
			v[1] = static_cast<float>(tmp0[1]);
			v[2] = static_cast<float>(tmp1[1]);
			v[3] = static_cast<float>(tmp1[0]);
		}
		else
		{
			v[0] = this->GetValue(xx, yy, data);
			v[1] = this->GetValue(xx + 1, yy, data);
			v[2] = this->GetValue(xx + 1, yy + 1, data);
			v[3] = this->GetValue(xx, yy + 1, data);
		}
		*/

		v[0] = this->GetValue(xx, yy, data);
		v[1] = this->GetValue(xx + 1, yy, data);
		v[2] = this->GetValue(xx + 1, yy + 1, data);
		v[3] = this->GetValue(xx, yy + 1, data);

		//obtain threshold based key
		//if value is same as threshold - consider it part of key
		uint32_t key = 0;
		if (v[0] <= t) key |= 1;
		if (v[1] <= t) key |= 2;
		if (v[2] <= t) key |= 4;
		if (v[3] <= t) key |= 8;

		//early "termination" - no contour cross the square
		if ((key == 0) || (key == 15))
		{			
			continue;
		}

		//precache indexes of pixels and their 2D positions		
		indexes[0] = { this->GetIndexFromPosition(xx, yy, data.GetWidth()), xx, yy };
		indexes[1] = { this->GetIndexFromPosition(xx + 1, yy, data.GetWidth()), xx + 1, yy };
		indexes[2] = { this->GetIndexFromPosition(xx + 1, yy + 1, data.GetWidth()), xx + 1, yy + 1 };
		indexes[3] = { this->GetIndexFromPosition(xx, yy + 1, data.GetWidth()), xx, yy + 1 };


		//based on key, get edges that are crossed within the square
		auto edgesInfo = edges[key];


		//iterate all crossed edges and
		//for each edge calculate contour point location		
		for (size_t i = 0; i < 4; i += 2)
		{
			auto ei = edgesInfo[i];
			if (ei == -1) break;


			auto edge = edgesVertices[ei];

			//calculate edge index - neighbors have same edge index and we dont
			//need to calculate point twice
			uint64_t edgeKey = CalculateEdgeKey(indexes[edge[0]][0], indexes[edge[1]][0]);

			auto it = visitedEdges.find(edgeKey);
			if (it != visitedEdges.end())
			{
				//first edge already processed 
				//take the second edge on line
				ei = edgesInfo[i + 1];
				edge = edgesVertices[ei];
				edgeKey = CalculateEdgeKey(indexes[edge[0]][0], indexes[edge[1]][0]);
			}

			it = visitedEdges.find(edgeKey);
			if (it == visitedEdges.end())
			{
				auto nextCell = edgesNextCell[ei];
				xx += nextCell[0];
				yy += nextCell[1];

				//calculate percents of threshold within the edge
				auto mm = std::minmax(v[edge[0]], v[edge[1]]);
				double p = double(t - mm.first) / (mm.second - mm.first);

				//Lerp position of new point
				double px = MathUtils::Lerp((double)indexes[edge[0]][1], (double)indexes[edge[1]][1], p);
				double py = MathUtils::Lerp((double)indexes[edge[0]][2], (double)indexes[edge[1]][2], p);
								
				pl.AddPoint({ px, py });
				visitedEdges.insert(edgeKey);

				canContinue = true;
				break;
			}
		}
	}
	
	
	if (pl.Size() > 1)
	{
		linesAll.push_back(std::move(pl));
	}
}


template void MarchingSquares::Run(const Image2d<uint8_t> & data, float threshold);
template void MarchingSquares::Run(const Image2d<float> & data, float threshold);

template void MarchingSquares::Run(const Image2d<uint8_t> & data, int xx, int yy, int w, int h, float threshold);
template void MarchingSquares::Run(const Image2d<float> & data, int xx, int yy, int w, int h, float threshold);

template float MarchingSquares::GetValue(int x, int y, const Image2d<uint8_t> & data) const;
template float MarchingSquares::GetValue(int x, int y, const Image2d<float> & data) const;

template void MarchingSquares::TrackContour(int x, int y, const Image2d<uint8_t> & data, float threshold);
template void MarchingSquares::TrackContour(int x, int y, const Image2d<float> & data, float threshold);
