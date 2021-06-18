#ifndef COLOR_SPACE_H
#define COLOR_SPACE_H

template <typename T>
class Image2d;

#include <cstdint>
#include <cmath>

class ColorSpace 
{
public:

	enum class PixelFormat
	{
		NONE,
		GRAY,
		RGB,
		RGBA,
		XYZ,
		CIE_LUV,
		HSV,
		RG //special format - only 2 chanels (cannot be saved) - for compatibility with OpenGL GL_RG
	};

	
		
	static size_t GetChannelsCount(PixelFormat pf);	

	template <typename T>
	static PixelFormat GetFormatFromChannelsCount(size_t c)
	{
		if constexpr (std::is_same<T, uint8_t>::value)
		{
			switch (c)
			{
			case 1: return PixelFormat::GRAY;
			case 2: return PixelFormat::RG;
			case 3: return PixelFormat::RGB;
			case 4: return PixelFormat::RGBA;
			default: return PixelFormat::NONE;
			}
		}
		else 
		{
			switch (c)
			{
			case 1: return PixelFormat::GRAY;			
			case 3: return PixelFormat::XYZ;			
			default: return PixelFormat::NONE;
			}
		}
	}

	static Image2d<uint8_t> ConvertGrayToRgb(const Image2d<uint8_t>& input);


};

#endif

