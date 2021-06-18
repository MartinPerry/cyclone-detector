#ifndef PRESSURE_FINDER_H
#define PRESSURE_FINDER_H


class PressureFinderSaver;

#include <functional>
#include <vector>
#include <list>

#include "./MapProjections/GeoCoordinate.h"


#include "./RasterData/Image2d.h"
#include "./RasterData/ImageFilters.h"

#include "./PressureExtrema.h"


class PressureFinder 
{
public:
	
	static std::string debugName;
	static std::string debugDirName;


	PressureFinder(const char * fileName, int w, int h);
	PressureFinder(const std::vector<float> & data, int w, int h);
	PressureFinder(Image2d<float>&& input);

	~PressureFinder();

	double GetContoursStep() const;
	int GetSystemsCount() const;

	const std::vector<PressureExtrema>& GetExtremas() const;
	const std::vector<Contour>& GetContours() const;
	const Image2d<float>& GetRawData() const;

	float GetDataMin() const;
	float GetDataMax() const;
	
	void SetAreaLatitudeCorrectionFactor(double f);
	void SetMinCentersDistanceThresholdPixel(int v);
	void SetSmallAreaThreshold(AreaThreshold v);
	void SetExtremaCorrectRatio(double v);
	void SetExtremaMaskRadius(int r);
	void SetContoursStep(double step);
	void Run();

	PressureFinderSaver GetSaver();

	friend class PressureFinderSaver;
	
protected:
	struct Derivatives 
	{
		std::vector<float> signDx;
		std::vector<float> signDy;
	};

	Image2d<float> rawData;
	float minValue;
	float maxValue;

	std::vector<Contour> contours; //all contours
	std::vector<PressureExtrema> extremas; //only contours enclosing extrema
	
	double latCorrectionFactor;
	AreaThreshold smallAreaThreshold;

	int maskRadius; //radius of mask for extrema detection
	double contourStep; //step between contours in Pa

	double extremaCorrectRatio;

	int minCentersDistanceThresholdPixel;

	void FindContours();
	void BuildContoursHierarchy();

	double GetLatitudeAreaCorrection(double y);
	double CalcArea(const Contour& c, Projections::IProjectionInfo* proj);
	
	template <typename Proj>
	double CalcAreaProj(const Contour& c, Proj* proj);

	void InitExtremaMask(ConvolutionFilter2d & maskX,
		ConvolutionFilter2d & maskY);
	Derivatives CreateDerivatives();

	void FindExtrema();

	void RenderContours(const std::vector<Contour> & c) const;
	void RenderContours();
};



class PressureFinderSaver
{
public:
	PressureFinderSaver(const PressureFinder * pf);
	~PressureFinderSaver();


	void Save(const char * modelName, struct tm t,
		const char * dbFileName,
		std::function<void(int, int, Latitude &, Longitude &)> positionToGps);

protected:

	const PressureFinder * pf;

};

#endif
