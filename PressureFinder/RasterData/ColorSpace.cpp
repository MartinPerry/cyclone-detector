#include "./ColorSpace.h"

#include "./Image2d.h"

//=============================================================================================
// Helper methods
//=============================================================================================


/// <summary>
/// Get number of channels based on PixelFormat
/// </summary>
/// <param name="pf"></param>
/// <returns></returns>
size_t ColorSpace::GetChannelsCount(PixelFormat pf)
{
	switch (pf)
	{
	case PixelFormat::GRAY: return 1;
	case PixelFormat::RGB: return 3;
	case PixelFormat::RG: return 2;
	case PixelFormat::RGBA: return 4;
	case PixelFormat::XYZ: return 3;	
	case PixelFormat::CIE_LUV: return 3;
	case PixelFormat::HSV: return 3;
	default: return 0;
	}
}


/// <summary>
/// Convert Grayscale image to RGB
/// (Each channel has the same value)
/// </summary>
/// <param name="input"></param>
/// <returns></returns>
Image2d<uint8_t> ColorSpace::ConvertGrayToRgb(const Image2d<uint8_t>& input)
{
	size_t len = input.GetPixelsCount();
	std::vector<uint8_t> data;
	data.resize(len * 3);

	size_t outIndex = 0;
	for (size_t i = 0; i < len; i++)
	{
		const uint8_t gray = input.GetData()[i];

		data[outIndex++] = gray;
		data[outIndex++] = gray;
		data[outIndex++] = gray;

	}

	return Image2d<uint8_t>(input.GetWidth(),
		input.GetHeight(),
		std::move(data),
		PixelFormat::RGB);
}
