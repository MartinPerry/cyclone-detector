# Cyclone detector

Run program as \<program name\> -\<paramy> -\<param\> ...

*Note:
For the test program the input file is a single-channel grayscale image in range 0 - 255. 
The real pressure values are therefore scaled to be in this range.*

# Main input parameters:

- **file** 
	- default: "" 
	- Name of input file with raw pressure data (single byte)
- **vebrose**
	- default: false 
	- verbose mode
- **cont_step**
	- default: 10
	- Contour step **C** (in units based on the input file)
	 *Note: The contour step is the greyscale value step. The greyscale values can be converted to Pa / hPa based on interval mappings*
- **small_area_threshold**
	- default: 10000
	- Small area threshold **T**. Extrema must occupy at least area of this size in km^2. 
	  *Note: If the projection of input data is unknown, this value represents area in image pixels. However, in that case, `area_lat_correction` should be set as well.* 

## Additional input parameters:

- **proj_type**
	- default: eq
	- Type of input data projection. ? = unknown, eq = Equirectangular, me = Meractor. Input projection is ussed to calculate area in km^2.
- **min_lat**
	- default: -90
	- Minimal latitude of the image area AABB
- **max_lat**
	- default: 90
	- Maximal latitude of the image area AABB
- **min_lon**
	- default: -180
	- Minimal longitude of the image area AABB
- **max_lon**
	- default: 180
	- Maximal longitude of the image area AABB

- **extrema_correction_ratio**
	- default: 0.6
	- Pixel is considered as extrema if the ratio of extrema vs. non-extrema pixels inside the mask is above this value
- **min_dist**
	- default: 0  (disabled)
	- Minimal distance **D** between extrema. To disable this feature, set value to 0.
- **mask_radius**
	- default: 0 (auto-calculated)
	- Extrema mask radius **r**. Pixel size of the area where to look for extrema
- **area_lat_correction**
	- default: 0.0 (disabled)
	- !!! USE ONLY IF `small_area_threshold` is in pixels, not in km^2 !!!!
	- Correction factor for the latitude of the area. Usefull for whole world inputs.
	  On the Equator, the area is taken directly from `small_area_threshold`. 
	  On poles, `small_area_threshold` is multipled by this factor "
	  and between, the factor is linearly interpolated
	  To disable this feature, set factor to 0.0

If program is compiled without `DISABLE_DEBUG_OUTPUT` (default state) we can output intermediate results from the program

- **debug_dir**
	- default: ""
	- Directory where we want to store the intermediate results. If empty, disable this feature
- **debug_name**
	- default: "debug" 
	- Intermediate file name prefix


# Compilation

There are two ways how to compile the source codes:

## Using makefile
We use a simple makefile script. 
If you add `DISABLE_DEBUG_OUTPUT`, program will be compiled without the availability to generate intermediate results. 
This option will slightly improve performance, since some checks are ommited.

## Using Visual Studio 2019 solution
We have also provided Visual Studio Solution files with x64 Debug and Release build.
In this case in `VC++ Directories` add or remove include and library directory for OpenCV library. By default, these path are set to
`D:\opencv\build\include` and `D:\opencv\build\x64\vc15\lib` respectively. If you do not want OpenCV support, simple remove the variables.

If you add `DISABLE_DEBUG_OUTPUT` under `C/C++ -> Preprocessor -> Preprocessor Definition`, program will be compiled without the availability to generate intermediate results. 
This option will slightly improve performance, since some checks are ommited.


# Description

The program does not require to add any any 3rd party dependencies 
(however, you have the opportunity to enable OpenCV for some debug outputs). 
It is written in C++17. 
All source files required for the compialtion are included in the repository. 
The main algorithm is located in files `PressureExtrema` and `CycloneDetector`.

# Input files
This demo uses only single-channel grayscale images for the input. 

## Edit input files
However, as you can see in the code, the image is converted to float array. 
If you edit the source code, you may load the input file for example directly from GRIB or NetCDF file and use a float array.
In this case, assign your data to an image `Image2d<float> rawDataFloat` in `main.cpp`.
In the code of `main.cpp`, this can be seen around line 102

```C++
	// ======================================================================
	// NOTE =================================================================
	// ======================================================================
	 
	//If you want to load data from different source (GRIB, NetCDF)
	//change this to eg:	
	//Image2d<float> rawDataFloat = GetDataFromGrib(fileName);
	//where GetDataFromGrib is your method to load data

	Image2d<uint8_t> rawDataUint(fileName.c_str());
	if (rawDataUint.GetPixelsCount() == 0)
	{
		MY_LOG_ERROR("Input image %s is empty", fileName.c_str());
		return 0;
	}
	
	Image2d<float> rawDataFloat = rawDataUint.CreateAs<float>();
	// ======================================================================
	// ======================================================================
	// ======================================================================
```
We have choosed this because we did not want to add any 3rd party dependencies. Usually GRIB or NetCDF libraries
are large projects and precompiled libraries may not be compatible with all configurations.
This way, the demo program is independent on external libraries and can be used as it is.



https://guides.github.com/features/mastering-markdown/