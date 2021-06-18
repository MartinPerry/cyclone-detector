#include "ProjectionRenderer.h"

#if __has_include("./lodepng.h")
#	include "./lodepng.h"
#else
#	include "../Compression/3rdParty/lodepng.h"
#endif
#include "MapProjectionStructures.h"

using namespace Projections;

//=======================================================================
// Debug renderer
//
//=======================================================================


/// <summary>
/// dtor
/// </summary>
ProjectionRenderer::~ProjectionRenderer()
{
	if (!externalData)
	{
		delete[] rawData;
		rawData = nullptr;
	}
}

const uint8_t * ProjectionRenderer::GetRawData() const
{
	return this->rawData;
}


/// <summary>
/// Clear data
/// </summary>
void ProjectionRenderer::Clear()
{
	memset(rawData, 0, static_cast<int>(frame.w) * static_cast<int>(frame.h) * static_cast<int>(type));
}

void ProjectionRenderer::SetPixelVal(uint8_t v)
{
	this->pixelVal = v;
}

void ProjectionRenderer::SetRawDataTarget(uint8_t * target, RenderImageType targetType)
{
	if (targetType != this->type)
	{
		return;
	}

	if (target == nullptr)
	{
		this->rawData = nullptr;
		this->rawData = new uint8_t[static_cast<int>(frame.w) * static_cast<int>(frame.h) * static_cast<int>(type)];

		this->externalData = false;
		return;
	}
	if (!externalData)
	{
		delete[] rawData;
	}
	this->rawData = target;
	this->externalData = true;
}

/// <summary>
/// Save data to *.png file
/// </summary>
/// <param name="fileName"></param>
void ProjectionRenderer::SaveToFile(const char * fileName)
{
	if (type == RenderImageType::GRAY)
	{
		lodepng::encode(fileName, this->rawData,
			static_cast<int>(frame.w), static_cast<int>(frame.h),
			LodePNGColorType::LCT_GREY);
	}
	else if (type == RenderImageType::RGB)
	{
		lodepng::encode(fileName, this->rawData,
			static_cast<int>(frame.w), static_cast<int>(frame.h),
			LodePNGColorType::LCT_RGB);
	}
	else if (type == RenderImageType::RGBA)
	{
		lodepng::encode(fileName, this->rawData,
			static_cast<int>(frame.w), static_cast<int>(frame.h),
			LodePNGColorType::LCT_RGBA);
	}
}

void ProjectionRenderer::FillData(std::vector<uint8_t> & output)
{
	output.resize(static_cast<int>(frame.w) * static_cast<int>(frame.h) * static_cast<int>(type));
	memcpy(&output[0], this->rawData, output.size());
}


/// <summary>
/// Load input data from file
/// </summary>
/// <param name="filePath"></param>
/// <returns></returns>
std::string ProjectionRenderer::LoadFromFile(const char * filePath)
{
	FILE * f = NULL;  //pointer to file we will read in
	my_fopen(&f, filePath, "rb");
	if (f == NULL)
	{
		printf("Failed to open file: \"%s\"\n", filePath);
		return "";
	}

	fseek(f, 0L, SEEK_END);
	long size = ftell(f);
	fseek(f, 0L, SEEK_SET);

	char * data = new char[size + 1];
	fread(data, sizeof(char), size, f);
	fclose(f);

	data[size] = 0;
	std::string tmp = data;
	delete[] data;

	return tmp;
};

std::vector<std::string> ProjectionRenderer::Split(const std::string &s, char delim)
{
	std::stringstream ss;
	ss.str(s);
	std::string item;
	std::vector<std::string> elems;
	while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	}

	return elems;
}

/// <summary>
/// Load borders from file in format "lat;lon\n"
/// </summary>
/// <param name="fileName"></param>
/// <param name="useEveryNthPoint"></param>
void ProjectionRenderer::AddBorders(const char * fileName, int useEveryNthPoint)
{
	std::string content = LoadFromFile(fileName);


	std::vector<std::string> lines = Split(content, '\n');

	std::string keyName;
	int tmp = 0;

	for (size_t i = 0; i < lines.size(); i++)
	{
		std::vector<std::string> line = Split(lines[i], ';');
		if (line.size() <= 2)
		{
			continue;
		}
		if (line.size() > 5)
		{
			keyName = line[7];
			tmp = 0;
		}
		else
		{
			if (tmp % useEveryNthPoint != 0)
			{
				tmp++;
				continue;
			}

			std::string key = keyName;
			key += "_";
			key += line[3];


			Coordinate point;
			point.lon = Longitude::deg(atof(line[0].c_str()));
			point.lat = Latitude::deg(atof(line[1].c_str()));
			this->debugBorder[key].push_back(point);

			tmp++;
		}

	}


}

void ProjectionRenderer::DrawLine(Pixel<int> pp1, Pixel<int> pp2)
{
	this->CohenSutherlandLineClipAndDraw(pp1.x, pp1.y, pp2.x, pp2.y);

	//process world possible wrapping around and repeting
	{
		int p1x = pp1.x;
		int p2x = pp2.x;

		MyRealType nc = this->frame.repeatNegCount;
		int offset = static_cast<int>(frame.ww);
		while (nc > 0)
		{
			pp1.x -= offset;
			pp2.x -= offset;

			this->CohenSutherlandLineClipAndDraw(pp1.x, pp1.y, pp2.x, pp2.y);
			nc--;
		}

		pp1.x = p1x;
		pp2.x = p2x;

		MyRealType pc = this->frame.repeatPosCount;
		offset = static_cast<int>(frame.ww);
		while (pc > 0)
		{
			pp1.x += offset;
			pp2.x += offset;

			this->CohenSutherlandLineClipAndDraw(pp1.x, pp1.y, pp2.x, pp2.y);
			pc--;
		}
	}
}

/// <summary>
/// Draw borders
/// </summary>
void ProjectionRenderer::DrawBorders()
{	
	for (const auto & it : this->debugBorder)
	{
		const std::vector<Coordinate> & b = it.second;

		for (size_t i = 0; i < b.size(); i++)
		{
			Coordinate p = b[i % b.size()];
			Coordinate p1 = b[(i + 1) % b.size()];

			Pixel<int> pp1 = this->projectCallback(p);
			Pixel<int> pp2 = this->projectCallback(p1);
			
			this->DrawLine(pp1, pp2);
		}


	}
}

/// <summary>
/// Draw parallels
/// </summary>
void ProjectionRenderer::DrawParalells()
{	
	this->DrawParalells(10, 5);
}

void ProjectionRenderer::DrawParalells(MyRealType lonStep, MyRealType latStep)
{

	for (MyRealType lat = -90; lat <= 90; lat += latStep)
	{
		for (MyRealType lon = -180; lon <= 180 - lonStep; lon += lonStep)
		{
			Coordinate p;
			p.lat = Latitude::deg(lat);
			p.lon = Longitude::deg(lon);

			Coordinate p1;
			p1.lat = p.lat;
			p1.lon = Longitude::deg(lon + lonStep);

			Pixel<int> pp1 = this->projectCallback(p);
			Pixel<int> pp2 = this->projectCallback(p1);
			this->DrawLine(pp1, pp2);

			p1.lat = Latitude::deg(lat + latStep);
			p1.lon = p.lon;

			pp1 = this->projectCallback(p);
			pp2 = this->projectCallback(p1);
			this->DrawLine(pp1, pp2);
		}
	}
}

/// <summary>
/// Draw line from [start] to [end] with a stepCount
/// </summary>
/// <param name="start"></param>
/// <param name="end"></param>
/// <param name="stepCount"></param>
void ProjectionRenderer::DrawLine(Coordinate start, Coordinate end, int stepCount)
{
	MyRealType difLat = (end.lat.rad() - start.lat.rad());
	MyRealType difLon = (end.lon.rad() - start.lon.rad());

	MyRealType lonStep = difLon / stepCount;
	MyRealType latStep = difLat / stepCount;


	Coordinate p = start;

	for (int i = 0; i < stepCount; i++)
	{
		Coordinate p1 = p;
		p1.lat = Latitude::rad(p1.lat.rad() + latStep);
		p1.lon = Longitude::rad(p1.lon.rad() + lonStep);

		Pixel<int> pp1 = this->projectCallback(p);
		Pixel<int> pp2 = this->projectCallback(p1);

		this->DrawLine(pp1, pp2);

		p = p1;
	}
}

/// <summary>
/// Draw point - point is represented by small square
/// </summary>
/// <param name="p"></param>
void ProjectionRenderer::DrawPoint(Coordinate p, int size)
{	
	Pixel<int> center = this->projectCallback(p);

	Pixel<int> a = center;
	Pixel<int> b = center;
	Pixel<int> c = center;
	Pixel<int> d = center;

	a.x -= size; a.y -= size;
	b.x += size; b.y -= size;
	c.x += size; c.y += size;
	d.x -= size; d.y += size;

	this->DrawLine(a, b);
	this->DrawLine(b, c);
	this->DrawLine(c, d);
	this->DrawLine(d, a);

}

/// <summary>
/// Draw lines created by points
/// </summary>
/// <param name="points"></param>
void ProjectionRenderer::DrawLines(const std::vector<Coordinate> & points)
{
	if (points.size() <= 1)
	{
		return;
	}

	for (size_t i = 0; i < points.size() - 1; i++)
	{
		this->DrawLine(points[i], points[i + 1]);
	}

}




/// <summary>
/// Set pixel value
/// </summary>
/// <param name="p"></param>
/// <param name="val"></param>
void ProjectionRenderer::SetPixel(const Pixel<int> & p, uint8_t val)
{	
	int index = (p.x + p.y * static_cast<int>(frame.w)) * static_cast<int>(type);
	for (int k = 0; k < static_cast<int>(type); k++)
	{
		this->rawData[index + k] = val;
	}	
}

/// <summary>
/// Codes for Cohen-Sutherlan
/// </summary>
/// <param name="x"></param>
/// <param name="y"></param>
/// <returns></returns>
int ProjectionRenderer::ComputeOutCode(MyRealType x, MyRealType y)
{
	int code = INSIDE;

	if (x < 0) code |= LEFT;
	else if (x >= static_cast<int>(frame.w)) code |= RIGHT;

	if (y < 0) code |= BOTTOM;
	else if (y >= static_cast<int>(frame.h)) code |= TOP;

	return code;
}

/// <summary>
/// Cohen-Sutherland line clipping
/// </summary>
/// <param name="x0"></param>
/// <param name="y0"></param>
/// <param name="x1"></param>
/// <param name="y1"></param>
void ProjectionRenderer::CohenSutherlandLineClipAndDraw(MyRealType x0, MyRealType y0, MyRealType x1, MyRealType y1)
{
	// compute outcodes for P0, P1, and whatever point lies outside the clip rectangle
	int outcode0 = ComputeOutCode(x0, y0);
	int outcode1 = ComputeOutCode(x1, y1);
	bool accept = false;

	MyRealType xmin = 0;
	MyRealType xmax = static_cast<int>(frame.w) - 1;

	MyRealType ymin = 0;
	MyRealType ymax = static_cast<int>(frame.h) - 1;

	while (true) {
		if (!(outcode0 | outcode1)) { // Bitwise OR is 0. Trivially accept and get out of loop
			accept = true;
			break;
		}
		else if (outcode0 & outcode1) { // Bitwise AND is not 0. (implies both end points are in the same region outside the window). Reject and get out of loop
			break;
		}
		else {
			// failed both tests, so calculate the line segment to clip
			// from an outside point to an intersection with clip edge
			MyRealType x = 0;
			MyRealType y = 0;

			// At least one endpoint is outside the clip rectangle; pick it.
			int outcodeOut = outcode0 ? outcode0 : outcode1;

			// Now find the intersection point;
			// use formulas y = y0 + slope * (x - x0), x = x0 + (1 / slope) * (y - y0)
			if (outcodeOut & TOP) {           // point is above the clip rectangle
				x = x0 + (x1 - x0) * (ymax - y0) / (y1 - y0);
				y = ymax;
			}
			else if (outcodeOut & BOTTOM) { // point is below the clip rectangle
				x = x0 + (x1 - x0) * (ymin - y0) / (y1 - y0);
				y = ymin;
			}
			else if (outcodeOut & RIGHT) {  // point is to the right of clip rectangle
				y = y0 + (y1 - y0) * (xmax - x0) / (x1 - x0);
				x = xmax;
			}
			else if (outcodeOut & LEFT) {   // point is to the left of clip rectangle
				y = y0 + (y1 - y0) * (xmin - x0) / (x1 - x0);
				x = xmin;
			}

			// Now we move outside point to intersection point to clip
			// and get ready for next pass.
			if (outcodeOut == outcode0) {
				x0 = x;
				y0 = y;
				outcode0 = ComputeOutCode(x0, y0);
			}
			else {
				x1 = x;
				y1 = y;
				outcode1 = ComputeOutCode(x1, y1);
			}
		}
	}
	if (accept)
	{
		// Following functions are left for implementation by user based on
		// their platform (OpenGL/graphics.h etc.)
		//DrawRectangle(xmin, ymin, xmax, ymax);
		//LineSegment(x0, y0, x1, y1);
		//DrawLine(x0, y0, x1, y1);
		Pixel<int> start, end;
		start.x = static_cast<int>(x0); start.y = static_cast<int>(y0);
		end.x = static_cast<int>(x1); end.y = static_cast<int>(y1);

		lineBresenhamCallback(start, end,
			[&](int x, int y) -> void {			
			this->SetPixel({ x, y }, this->pixelVal);
		});
	}
}


//=====

