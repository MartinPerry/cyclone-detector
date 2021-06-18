#ifndef PROJECTION_DEBUGGER_H
#define PROJECTION_DEBUGGER_H

#include <cstdint>
#include <vector>
#include <string>
#include <unordered_map>

#include "./ProjectionInfo.h"
#include "./Reprojection.h"

namespace Projections
{
	
	class ProjectionRenderer
	{
	public:

		enum class RenderImageType 
		{
			GRAY = 1,
			RGB = 3,
			RGBA = 4
		};

		template <typename ReprojType>
		static void ReprojectImage(const uint8_t * fromData, RenderImageType fromType, 
			uint8_t * toData, RenderImageType toType, const Reprojection<ReprojType> & reproj);

		template <typename Proj>
		ProjectionRenderer(Proj * proj, RenderImageType type = RenderImageType::GRAY);
		~ProjectionRenderer();

		const uint8_t * GetRawData() const;

		template <typename Proj>
		void SetProjection(Proj * proj);
		void Clear();

		void SetPixelVal(uint8_t v);
		void SetRawDataTarget(uint8_t * target, RenderImageType targetType);

		void AddBorders(const char * fileName, int useEveryNthPoint = 1);
		void DrawBorders();
		void DrawParalells();
		void DrawParalells(MyRealType lonStep, MyRealType latStep);

		void DrawLine(Coordinate start, Coordinate end, int stepCount = 20);

		void DrawPoint(Coordinate p, int size = 5);

		void DrawLines(const std::vector<Coordinate> & points);

		template <typename Proj>
		void DrawImage(const uint8_t * imData, RenderImageType imType, int w, int h, Proj * imProj);

		template <typename ReprojType>
		void DrawImage(const uint8_t * imData, RenderImageType imType, const Reprojection<ReprojType> & reproj);

		void SetPixel(const Pixel<int> & p, uint8_t val);

		void FillData(std::vector<uint8_t> & output);
		void SaveToFile(const char * fileName);

	private:
		static const int INSIDE = 0; // 0000
		static const int LEFT = 1;   // 0001
		static const int RIGHT = 2;  // 0010
		static const int BOTTOM = 4; // 0100
		static const int TOP = 8;    // 1000


		uint8_t * rawData;
		RenderImageType type;
		bool externalData;

		uint8_t pixelVal;

		ProjectionFrame frame;
		std::function<Pixel<int>(const Coordinate & c)> projectCallback;
		std::function<Coordinate(const Pixel<int> & p)> projectInverseCallback;
		std::function<void(const Pixel<int> & start, const Pixel<int> & end, std::function<void(int x, int y)> callback)> lineBresenhamCallback;
	
		std::unordered_map<std::string, std::vector<Coordinate> > debugBorder;

		int ComputeOutCode(MyRealType x, MyRealType y);
		void CohenSutherlandLineClipAndDraw(MyRealType x0, MyRealType y0, MyRealType x1, MyRealType y1);

		void DrawLine(Pixel<int> pp1, Pixel<int> pp2);

		std::vector<std::string> Split(const std::string &s, char delim);
		std::string LoadFromFile(const char * filePath);

	};


	/// <summary>
	/// ctor
	/// </summary>
	/// <param name="proj"></param>
	template <typename Proj>
	ProjectionRenderer::ProjectionRenderer(Proj * proj, RenderImageType type)
		: rawData(nullptr), externalData(false), type(type), pixelVal(255)
	{
		this->SetProjection(proj);
	};

	/// <summary>
	/// Set current projection - rewrites projection from ctor
	/// It also clears current data, because new projection can have
	/// different frame size
	/// </summary>
	/// <param name="proj"></param>
	template <typename Proj>
	void ProjectionRenderer::SetProjection(Proj * proj)
	{
		frame = proj->GetFrame();
		projectCallback = [proj](const Coordinate & c) -> Pixel<int> {
            return proj->template Project<int>(c);
		};
		projectInverseCallback = [proj](const Pixel<int> & p) -> Coordinate {
			return proj->ProjectInverse(p);
		};

		lineBresenhamCallback = [proj](const Pixel<int> & start,
			const Pixel<int> & end,
			std::function<void(int x, int y)> callback) -> void {
			proj->LineBresenham(start, end, callback);
		};


		if (!externalData)
		{
			delete[] rawData;
			rawData = new uint8_t[static_cast<int>(frame.w) * static_cast<int>(frame.h) * static_cast<int>(type)];
		}
		this->Clear();
	};

	/// <summary>
	/// Draw image based on input projection
	/// For each pixel [x, y] -> calculate inverse based on this projection -> [lat, lon]
	/// Use [lat, lon] on imProj to get pixel coordinates [xx, yy]
	/// Map imData[xx, yy] to currentData[x, y]
	/// </summary>
	/// <param name="imData"></param>
	/// <param name="imType></param>
	/// <param name="w"></param>
	/// <param name="h"></param>
	/// <param name="imProj"></param>
	template <typename Proj>
	void ProjectionRenderer::DrawImage(const uint8_t * imData, RenderImageType imType, int w, int h, Proj * imProj)
	{
		if (imType != this->type)
		{
			return;
		}

		for (int y = 0; y < static_cast<int>(frame.h); y++)
		{
			for (int x = 0; x < static_cast<int>(frame.w); x++)
			{

				Coordinate cc = this->projectInverseCallback({ x,y });
                Pixel<int> p = imProj->template Project<int>(cc);

				if (p.x < 0) continue;
				if (p.y < 0) continue;
				if (p.x >= w) continue;
				if (p.y >= h) continue;

				int inIndex = (p.x + p.y * w) * static_cast<int>(type);
				int outIndex = x + y * static_cast<int>(frame.w) * static_cast<int>(type);

				for (int k = 0; k < static_cast<int>(type); k++)
				{					
					this->rawData[outIndex + k] = imData[inIndex + k];
				}

			}
		}

	};

	/// <summary>
	/// Draw image based on re-projection pixels
	/// e.g.: currentData[index] = imData[reprojection[index]]
	/// </summary>
	/// <param name="imData"></param>
	/// <param name="imType"></param>
	/// <param name="reproj"></param>
	template <typename ReprojType>
	void ProjectionRenderer::DrawImage(const uint8_t * imData, RenderImageType imType,
		const Reprojection<ReprojType> & reproj)
	{
		ProjectionRenderer::ReprojectImage(imData, imType, this->rawData, this->type, reproj);
	}

	template <typename ReprojType>
	void ProjectionRenderer::ReprojectImage(const uint8_t * fromData, RenderImageType fromType,
		uint8_t * toData, RenderImageType toType, const Reprojection<ReprojType> & reproj)
	{

		if (fromType != toType)
		{
			return;
		}

		for (int y = 0; y < reproj.outH; y++)
		{
			int yw = y * reproj.outW;
			for (int x = 0; x < reproj.outW; x++)
			{
				int index = x + yw;

				const auto& px = reproj.pixels[index];

				if ((px.x == -1) || (px.y == -1))
				{
					continue;
				}

				
				int origIndex = (px.x + px.y * reproj.inW) * static_cast<int>(fromType);
				int outIndex = index * static_cast<int>(toType);

				for (int k = 0; k < static_cast<int>(toType); k++)
				{
					toData[outIndex + k] = fromData[origIndex + k];
				}
			}
		}
	}
};

#endif

