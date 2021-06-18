//===================================================================================
// main
//===================================================================================


#include "./PressureFinder.h"

#include "./Utils/CmdParser.h"
#include "./Utils/Logger.h"

#include "./MapProjections/Projections/Equirectangular.h"
#include "./MapProjections/Projections/Mercator.h"

int main(int argc, char ** argv)
{
#if defined (_DEBUG) || defined (DEBUG)
	printf("Running DEBUG build => tests speed will be incorrect.\n");
#endif
	
	std::string fileName = "D://20160421_1200_world.png";
	std::string debugOutputDir = "D://";
	std::string debugOutputName = "debug";
	bool verbose = true;
	double contourStep = 10;
	int extremaMaskRadius = 0;
	int minCenterDist = 0;
	int smallAreaThreshold = 10'000;	
	
	double latMin = -90.0;
	double latMax = 90.0;
	double lonMin = -180.0;
	double lonMax = 180.0;
	std::string projType = "eq";

	double extremaCorRatio = 0.6;
	double areaLatCorFactor = 0.0;

	CmdParser cmd(argc, argv);	
	cmd.AddTarget({ "file", "Name of input file with raw pressure data (single byte)", true }, &fileName);
	
	cmd.AddTarget({ "verbose", "Verbose mode (default: false)", true }, &verbose);
		
	cmd.AddTarget({ "cont_step", "Contour step 'C' (in units based on input file - grayscale-level values in this demo) (default: 10)", true }, &contourStep);
		
	cmd.AddTarget({ "small_area_threshold", 
					"Small area threshold 'T'. Extrema must have occupy at least area of this size in km^2. "
					"Note: If the projection of input data is unknown, this value represents area in image pixels. "
					"However, in thats case, area_lat_correction should be set as well. "
					"(default: 10'000)", 
		true }, &smallAreaThreshold);
	
	cmd.AddTarget({ "mask_radius", "Extrema mask radius 'r'. Size of area where to look for extrema (default: 0 = automatic)", true }, &extremaMaskRadius);

	cmd.AddTarget({ "min_dist", "Minimal distance 'D' between extrema. To disable this feature, set value to 0. (default: 0)", true }, &minCenterDist);

	cmd.AddTarget({ "min_lat", "Minimal latitude of the image area AABB. (default: -90)", true }, &latMin);
	cmd.AddTarget({ "max_lat", "Maximal latitude of the image area AABB. (default: 90)", true }, &latMax);
	cmd.AddTarget({ "min_lon", "Minimal longitude of the image area AABB. (default: -180)", true }, &lonMin);
	cmd.AddTarget({ "max_lon", "Maximal longitude of the image area AABB. (default: 180)", true }, &lonMax);
	cmd.AddTarget({ "proj_type", "Type of input data projection. ? = unknown, eq = Equirectangular, me = Meractor (default: eq)", true }, &projType);

	//=======================================================

	
	//if area is set in pixels instead of km2, use this to correct projection distortion
	//this is only estimate of the correction
	cmd.AddTarget({ "area_lat_correction", 
					"!!! USE ONLY IF small_area_threshold is in pixels, not ik km^2 !!!! "
					"Correction factor for area latitude. Usefull for whole world inputs. "
					"On equator, area is taken directly from small_area_threshold. On poles, it is multipled by this factor "
					"and between, the factor is linearly interpolated "
					"To disable this feature, set factor to 0.0 "
					"(default: 0.0)", 
		true }, &areaLatCorFactor);
	

	cmd.AddTarget({ "extrema_correction_ratio", 
					"Pixel is considered extrema if ratio of neighborhood "
					"pixels that correspond to mask against 'incorrect' pixels is above extrema_correction_ratio "
					"(default: 0.6)", 
		true }, &extremaCorRatio);


#ifndef DISABLE_DEBUG_OUTPUT
	cmd.AddTarget({ "debug_dir", "Set directory for debug output", true }, &debugOutputDir);
	cmd.AddTarget({ "debug_name", "Set base name for debug output", true }, &debugOutputName);
#endif


	if (cmd.Parse() == false)
	{
		MyUtils::Logger::GetInstance()->EnableInfoLogging(MyUtils::Logger::LOG_STDOUT);
		cmd.PrintHelp();
		return 0;
	}

	if (verbose)
	{
		MyUtils::Logger::GetInstance()->EnableInfoLogging(MyUtils::Logger::LOG_STDOUT);
	}


	Image2d<uint8_t> rawDataUint(fileName.c_str());
	if (rawDataUint.GetPixelsCount() == 0)
	{
		MY_LOG_ERROR("Input image %s is empty", fileName.c_str());
		return 0;
	}

	//If you want to load data from different source (GRIB, NetCDF)
	//change this
	Image2d<float> rawDataFloat = rawDataUint.CreateAs<float>();


	PressureFinder pf(std::move(rawDataFloat));

#ifndef DISABLE_DEBUG_OUTPUT	
	pf.debugName = debugOutputName;
	pf.debugDirName = debugOutputDir;
#endif

	//pf.SetAreaLatitudeCorrectionFactor(areaLatCorFactor);

	AreaThreshold th;

	Projections::Coordinate bbMin, bbMax;
	bbMin.lat = Latitude::deg(latMin); bbMin.lon = Longitude::deg(lonMin);
	bbMax.lat = Latitude::deg(latMax); bbMax.lon = Longitude::deg(lonMax);

	if (projType == "eq")
	{
		Projections::Equirectangular* eq = new Projections::Equirectangular();
		eq->SetFrame(bbMin, bbMax,
			rawDataFloat.GetWidth(), rawDataFloat.GetHeight(),
			Projections::STEP_TYPE::PIXEL_BORDER, false);

		th = AreaThreshold::CreateKm2(smallAreaThreshold, eq);
	}
	else if (projType == "me")
	{
		Projections::Mercator* me = new Projections::Mercator();
		me->SetFrame(bbMin, bbMax,
			rawDataFloat.GetWidth(), rawDataFloat.GetHeight(),
			Projections::STEP_TYPE::PIXEL_BORDER, false);

		th = AreaThreshold::CreateKm2(smallAreaThreshold, me);
	}
	else if (projType == "?")
	{		
		th = AreaThreshold::CreatePixels(smallAreaThreshold);
		pf.SetAreaLatitudeCorrectionFactor(areaLatCorFactor);
	}

	pf.SetMinCentersDistanceThresholdPixel(minCenterDist);	
	pf.SetSmallAreaThreshold(th);
	pf.SetExtremaCorrectRatio(extremaCorRatio);

	pf.SetContoursStep(contourStep);
	pf.SetExtremaMaskRadius(extremaMaskRadius);
	pf.Run();


	auto extremas = pf.GetExtremas();

	int count = 0;
	for (const auto& c : extremas)
	{
		if (c.IsValid() == false)
		{
			continue;
		}
		count++;
	}
	MY_LOG("Extrema count: %d\n", count);
	

	for (const auto& c : extremas)
	{
		if (c.IsValid() == false)
		{
			continue;
		}
		MY_LOG("%f %f\n", c.center.x, c.center.y);
	}
	
	return 0;	
}
