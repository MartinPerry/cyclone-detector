#ifndef PNG_LOADER_H
#define PNG_LOADER_H

//#define HAVE_LIBPNG

struct IFile;

#include <vector>
#include <cstdint>

#if __has_include (<libpng/png.h>)
#    include <libpng/png.h>
#    define HAVE_LIBPNG 1
#elif __has_include (<libpng15/png.h>)
#    include <libpng15/png.h>
#    define HAVE_LIBPNG 1
#endif

class PNGLoader
{
public:
	enum class USED_LIBRARY 
	{ 
#ifdef HAVE_LIBPNG
		LIBPNG = 0, 
#endif
		LODEPNG = 1
	};

	typedef struct RGBA
	{
		union
		{
			struct 
			{
				uint8_t r;
				uint8_t g;
				uint8_t b;
				uint8_t a;
			};
			uint8_t _rgba[4];
		};
	} RGBA;

	typedef struct DecompressedImage
	{
		unsigned w;
		unsigned h;
		unsigned channelsCount;
		unsigned bitDepth;

		std::vector<uint8_t> data;
		std::vector<RGBA> palette;
		bool grayScalePallete;

	} DecompressedImage;

	PNGLoader();
	PNGLoader(USED_LIBRARY lib);
	~PNGLoader();

	void SetKeepPalette(bool val);
	
	DecompressedImage DecompressFromMemory(uint8_t * mem, size_t memSize);
	DecompressedImage DecompressFromFile(const char * fileName);
	DecompressedImage DecompressFromFile(IFile * file);

private:

	typedef struct PngRawData
	{
		uint8_t * data;
		size_t offset;

	} PngRawData;

	typedef struct PngFileData
	{
		IFile * file;
		size_t offset;

	} PngFileData;


	const USED_LIBRARY lib;
	const size_t PNG_SIG_SIZE = 8;

#ifdef HAVE_LIBPNG
	png_structp pngPtr;
	png_infop infoPtr;
	png_bytep * rowPtrs;
#endif

	bool keepPalette;

	void Release();

	DecompressedImage DecompressWithLodePNG(uint8_t * mem, size_t memSize);
#ifdef HAVE_LIBPNG
	DecompressedImage DecompressWithLibPNG(uint8_t * mem, size_t memSize);
	DecompressedImage DecompressWithLibPNG(IFile * file);

	bool InitLibPNG();
	void LibPNGReadHeader(DecompressedImage & d);
	void LibPNGReadData(DecompressedImage & d);
	void LibPNGReadPalette(DecompressedImage & d);		
	bool LibPNGIsPalleteGreyScale();


    static void UserWarningFn(png_structp png_ptr, png_const_charp warning_msg);
	static void UserReadDataMemory(png_structp pngPtr, png_bytep data, png_size_t length);
	static void UserReadDataFile(png_structp pngPtr, png_bytep data, png_size_t length);
#endif
};

#endif
