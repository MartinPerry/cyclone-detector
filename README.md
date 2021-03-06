# Cyclone detector

Run program as \<program name\> -\<param name\> -\<param name\> ...

*Note:
For the test program the input file is a single-channel grayscale image in range 0 - 255. 
The real pressure values are therefore scaled to be in this range.*

# Main input parameters:

- **file** 
	- default: "" 
	- Name of input file with the raw pressure data (single byte)
- **vebrose**
	- default: false 
	- verbose mode
- **cont_step**
	- default: 10
	- Contour step **C** (in units based on the input file)
	 *Note: The contour step is the grayscale value step for single-channel grayscale image. The grayscale values can be converted to Pa / hPa based on interval mappings*
- **small_area_threshold**
	- default: 10000
	- Small area threshold **T**. Extrema must occupy at least area of this size in km^2. 
	  *Note: If the projection of input data is unknown, this value represents area in image pixels. However, in that case, `area_lat_correction` should be set as well.* 

## Additional input parameters:

- **proj_type**
	- default: eq
	- Type of input data projection. ? = unknown, eq = Equirectangular, me = Meractor. Input projection is used to calculate area in km^2.
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
	- Extrema mask radius **r**. Pixel size of the area where to look for extrema.
- **area_lat_correction**
	- default: 0.0 (disabled)
	- !!! USE ONLY IF `small_area_threshold` is in pixels, not in km^2 !!!!
	- Correction factor for the latitude of the area. Useful for the whole world inputs.
	  On the Equator, the area is taken directly from `small_area_threshold`. 
	  On poles, `small_area_threshold` is multiplied by this factor and between, the factor is linearly interpolated. To disable this feature, set factor to 0.0.

If program is compiled without `DISABLE_DEBUG_OUTPUT` (default state) we can output intermediate results from the program.

- **debug_dir**
	- default: ""
	- Directory where we want to store the intermediate results. If empty, disable this feature.
- **debug_name**
	- default: "debug" 
	- Intermediate file name prefix. Each debug file name will have this prefix and the rest of the name is taken form the input name.


# Compilation

There are two ways how to compile the source codes:

## Using makefile
We use a simple makefile script. 
If you add `DISABLE_DEBUG_OUTPUT`, program will be compiled without the availability to generate intermediate results. 
This option will slightly improve performance since some checks are omitted.

## Using Visual Studio 2019 solution
We have also provided Visual Studio Solution files with x64 Debug and Release build.
In this case in `VC++ Directories` add or remove include and library directory for OpenCV library. By default, these paths are set to
`D:\opencv\build\include` and `D:\opencv\build\x64\vc15\lib` respectively. If you do not want OpenCV support, simple remove the variables.

If you add `DISABLE_DEBUG_OUTPUT` under `C/C++ -> Preprocessor -> Preprocessor Definition`, program will be compiled without the availability to generate intermediate results. 
This option will slightly improve performance since some checks are omitted.


# Description

The program does not require to add any 3rd party dependencies 
(however, you have the opportunity to enable OpenCV for some debug outputs). 
It is written in C++17. 
All source files required for the compilation are included in the repository. 
The main algorithm is located in files `PressureExtrema` and `CycloneDetector`.

# Input files
This demo uses only single-channel grayscale images for the input. 

Since most of the NWP models have outputs in GRIB or NetCDF, 
we have provided two simple Python script `grib_to_png.py` and `netcdf_to_png.py` to convert single-channel GRIB/NetCDF file to the grayscale image.

*Usage:* 
Specify input file `--input_file`. Other parameters are optional. 
If `--out_w` and/or `--out_h` is set, output image is resized to this resolution. 
If only one of them is set, the second dimension is auto-calculated to keep aspect ratio of the input data.
If `--output_file` is set, image is stored in the file with this name (required extension is `png`, since the main program expects data as `png`).

*Note:*
If data are flipped, modify `values = np.flipud(values)`.


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
//or
//Image2d<float> rawDataFloat = GetDataFromNetCDF(fileName);
//where GetDataFromGrib/GetDataFromNetCDF is your method to load data

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
We have chosen this because we did not want to add any 3rd party dependencies. 
Usually, GRIB or NetCDF libraries are large projects and precompiled libraries may not be compatible with all configurations.
This way, the demo program is independent on external libraries and can be used as it is.

