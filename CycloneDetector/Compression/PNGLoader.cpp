#include "./PNGLoader.h"


#ifdef HAVE_LIBPNG
#	ifdef _WIN32
#		pragma comment(lib, "libpng16.lib")
#		pragma comment(lib, "zlib.lib")
#	endif
#endif

#include <stdexcept>

#include "./3rdParty/lodepng.h"

#include "../Utils/Logger.h"

//#include "../VFS/VFS.h"
#include "../Utils/FileWrapper.h"


PNGLoader::PNGLoader()  :
#ifdef HAVE_LIBPNG
	PNGLoader(USED_LIBRARY::LIBPNG)
#else
	PNGLoader(USED_LIBRARY::LODEPNG)
#endif
{
}

PNGLoader::PNGLoader(USED_LIBRARY lib) : 
	lib(lib), 
#ifdef HAVE_LIBPNG
	pngPtr(nullptr), infoPtr(nullptr), rowPtrs(nullptr),
#endif
	keepPalette(false)
{

}

PNGLoader::~PNGLoader()
{
	this->Release();	
}

void PNGLoader::Release()
{
#ifdef HAVE_LIBPNG
	delete[](png_bytep)rowPtrs;
	rowPtrs = nullptr;

	if (pngPtr != nullptr)
	{
		if (infoPtr != nullptr)
		{
			png_destroy_info_struct(pngPtr, &infoPtr);
			infoPtr = nullptr;
		}
		png_destroy_read_struct(&pngPtr, (png_infopp)0, (png_infopp)0);
		pngPtr = nullptr;
	}
#endif
}

/// <summary>
/// Keep data in palette
/// If set to false, data are unpacked from palette
/// and returned image is real colors
/// If set to true, palette is kept and returned
/// image is created from indexes to palette
/// </summary>
/// <param name="val"></param>
void PNGLoader::SetKeepPalette(bool val)
{
	keepPalette = val;
}

PNGLoader::DecompressedImage PNGLoader::DecompressFromMemory(uint8_t * mem, size_t memSize)
{
	this->Release();

	if (lib == USED_LIBRARY::LODEPNG)
	{
		return this->DecompressWithLodePNG(mem, memSize);
	}

#ifdef HAVE_LIBPNG
	if (lib == USED_LIBRARY::LIBPNG)
	{
		return this->DecompressWithLibPNG(mem, memSize);
	}
#endif	
    
    throw std::runtime_error("Unknown decompression library");
    
}

PNGLoader::DecompressedImage PNGLoader::DecompressFromFile(const char * fileName)
{
	RawFile rw(fileName);
	return this->DecompressFromFile(&rw);
}

PNGLoader::DecompressedImage PNGLoader::DecompressFromFile(IFile * file)
{
	this->Release();

	if (lib == USED_LIBRARY::LODEPNG)
	{
		uint8_t * buf;
		size_t bufSize = file->ReadAll((void **)&buf);
		auto dec = this->DecompressWithLodePNG(buf, bufSize);
		free(buf);
		return dec;
	}

#ifdef HAVE_LIBPNG
	if (lib == USED_LIBRARY::LIBPNG)
	{
		return this->DecompressWithLibPNG(file);
	}	
#endif
	

	DecompressedImage dec;
	dec.w = 0;
	dec.h = 0;
	dec.channelsCount = 0;
	dec.bitDepth = 0;

	return dec;
}


PNGLoader::DecompressedImage PNGLoader::DecompressWithLodePNG(uint8_t * mem, size_t memSize)
{
	DecompressedImage dec;

	lodepng::State pngState;
	pngState.decoder.color_convert = 0;
	unsigned error = lodepng::decode(dec.data, dec.w, dec.h, pngState, mem, memSize);

	if (error != 0)
	{
		dec.w = 0;
		dec.h = 0;
		dec.channelsCount = 0;
		dec.bitDepth = 0;
		dec.data.clear();
		return dec;
	}

	if (dec.data.size() == dec.w * dec.h * sizeof(uint8_t)) dec.channelsCount = 1;
	else if (dec.data.size() == dec.w * dec.h * 3 * sizeof(uint8_t)) dec.channelsCount = 3;
	else if (dec.data.size() == dec.w * dec.h * 4 * sizeof(uint8_t)) dec.channelsCount = 4;

	dec.bitDepth = dec.channelsCount * sizeof(uint8_t);

	return dec;
}


#ifdef HAVE_LIBPNG

void PNGLoader::UserWarningFn(png_structp png_ptr, png_const_charp warning_msg)
{
    MY_LOG_INFO("LibPNG > %s", warning_msg);
}

void PNGLoader::UserReadDataMemory(png_structp pngPtr, png_bytep data, png_size_t length)
{	
	PngRawData * raw = (PngRawData *)png_get_io_ptr(pngPtr);

	memcpy(data, raw->data + raw->offset, length);
	raw->offset += length;
}


void PNGLoader::UserReadDataFile(png_structp pngPtr, png_bytep data, png_size_t length)
{	
	IFile * file = (IFile *)png_get_io_ptr(pngPtr);
	
	file->Read(data, sizeof(char), length);
	//raw->offset += length * sizeof(png_bytep);
}

bool PNGLoader::InitLibPNG()
{
	//Here we create the png read struct. The 3 NULL's at the end can be used
	//for your own custom error handling functions, but we'll just use the default.
	//if the function fails, NULL is returned. Always check the return values!
    this->pngPtr = png_create_read_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, PNGLoader::UserWarningFn);

	if (!pngPtr)
	{

		MY_LOG_ERROR("ERROR: Couldn't initialize png read struct");
		return false; //Do your own error recovery/handling here
	}

	//Here we create the png info struct.
	//Note that this time, if this function fails, we have to clean up the read struct!
	this->infoPtr = png_create_info_struct(pngPtr);
	if (!infoPtr)
	{
		MY_LOG_ERROR("ERROR: Couldn't initialize png info struct");
		png_destroy_read_struct(&pngPtr, (png_infopp)0, (png_infopp)0);
		return false; //Do your own error recovery/handling here

	}
	
	return true;
}

PNGLoader::DecompressedImage PNGLoader::DecompressWithLibPNG(uint8_t * mem, size_t memSize)
{
	//http://www.piko3d.net/tutorials/libpng-tutorial-loading-png-files-from-streams/

	

	DecompressedImage dec;
	dec.w = 0;
	dec.h = 0;
	dec.channelsCount = 0;
	dec.bitDepth = 0;

	if (this->InitLibPNG() == false)
	{
		return dec;
	}
	

	if (setjmp(png_jmpbuf(pngPtr)))
	{
		//An error occurred, so clean up what we have allocated so far...
		png_destroy_read_struct(&pngPtr, &infoPtr, (png_infopp)0);

		if (rowPtrs != nullptr)
		{
			delete[] rowPtrs;
		}
	
		MY_LOG_ERROR("ERROR: An error occured while reading the PNG file");

		//Make sure you return here. libPNG will jump to here if something
		//goes wrong, and if you continue with your normal code, you might
		//End up with an infinite loop.				
		return dec;
	}

	PngRawData raw;
	raw.data = mem;
	raw.offset = PNG_SIG_SIZE;

	png_set_read_fn(pngPtr, (png_voidp)&raw,PNGLoader::UserReadDataMemory);

	png_set_sig_bytes(pngPtr, static_cast<int>(PNG_SIG_SIZE));

	this->LibPNGReadHeader(dec);
	


	this->LibPNGReadData(dec);

	return dec;
}

PNGLoader::DecompressedImage PNGLoader::DecompressWithLibPNG(IFile * file)
{
	DecompressedImage dec;
	dec.w = 0;
	dec.h = 0;
	dec.channelsCount = 0;
	dec.bitDepth = 0;

	if (this->InitLibPNG() == false)
	{
		return dec;
	}


	if (setjmp(png_jmpbuf(pngPtr))) 
	{
		// If we get here, we had problem reading the file
		// Free all of the memory associated with the png_ptr and info_ptr

		png_destroy_read_struct(&pngPtr, &infoPtr, nullptr);
		return dec;
	}

	
	png_set_read_fn(pngPtr, (png_voidp)file, PNGLoader::UserReadDataFile);

	//https://blog.nobel-joergensen.com/2010/11/07/loading-a-png-as-texture-in-opengl-using-libpng/
	//png_init_io(pngPtr, file);	

	png_set_sig_bytes(pngPtr, 0);
	
	this->LibPNGReadHeader(dec);

	this->LibPNGReadData(dec);

	return dec;
}

void PNGLoader::LibPNGReadHeader(DecompressedImage & d)
{
	//Now call png_read_info with our pngPtr as image handle, and infoPtr to receive the file info.
	png_read_info(pngPtr, infoPtr);

	png_uint_32 imgWidth = png_get_image_width(pngPtr, infoPtr);
	png_uint_32 imgHeight = png_get_image_height(pngPtr, infoPtr);

	//bits per CHANNEL! note: not per pixel!
	png_uint_32 bitdepth = png_get_bit_depth(pngPtr, infoPtr);
	
	//Number of channels
	//png_uint_32 channels = png_get_channels(pngPtr, infoPtr);

	//Color type. (RGB, RGBA, Luminance, luminance alpha... palette... etc)
	png_uint_32 color_type = png_get_color_type(pngPtr, infoPtr);
	
	
	if (keepPalette)
	{		
		if (color_type == PNG_COLOR_TYPE_PALETTE)
		{
			this->LibPNGReadPalette(d);
		}
	}
	else 
	{
		if (bitdepth < 8)
		{
			png_set_expand_gray_1_2_4_to_8(pngPtr);
		}

		if (color_type == PNG_COLOR_TYPE_PALETTE)
		{
			png_set_palette_to_rgb(pngPtr);

			if (this->LibPNGIsPalleteGreyScale())
			{
				png_set_rgb_to_gray(pngPtr, 1, 0.0, 0.0);
			}			
		}

		/*if the image has a transparency set.. convert it to a full Alpha channel..*/
		if (png_get_valid(pngPtr, infoPtr, PNG_INFO_tRNS))
		{
			png_set_tRNS_to_alpha(pngPtr);
		}
	}

		

	//We don't support 16 bit precision.. so if the image Has 16 bits per channel
	//precision... round it down to 8.
	if (bitdepth == 16)
	{
		png_set_strip_16(pngPtr);				
	}

	png_read_update_info(pngPtr, infoPtr);

	d.w = imgWidth;
	d.h = imgHeight;
	d.channelsCount = png_get_channels(pngPtr, infoPtr);
	d.bitDepth = png_get_bit_depth(pngPtr, infoPtr);	
}

bool PNGLoader::LibPNGIsPalleteGreyScale()
{
	png_colorp palette;
	int num_palette;
	png_get_PLTE(pngPtr, infoPtr, &palette, &num_palette);

	png_bytep trans;
	int num_trans;
	png_color_16p col;
	png_uint_32 retTRNS = png_get_tRNS(pngPtr, infoPtr, &trans, &num_trans, &col);

	if (retTRNS == PNG_INFO_tRNS)
	{		
		return false;
		/*
		for (int i = 0; i < num_palette; i++)
		{
			RGBA rgba;
			rgba.r = palette[i].red;
			rgba.g = palette[i].green;
			rgba.b = palette[i].blue;

			if (rgba.r != rgba.g) return false;
			else if (rgba.r != rgba.b) return false;
			else if (rgba.g != rgba.b) return false;

			rgba.a = (i < num_trans) ? trans[i] : 255;
			d.palette.push_back(std::move(rgba));
		}
		*/
	}
	else
	{		
		for (int i = 0; i < num_palette; i++)
		{
			RGBA rgb;
			rgb.r = palette[i].red;
			rgb.g = palette[i].green;
			rgb.b = palette[i].blue;

			if (rgb.r != rgb.g) return false;
			else if (rgb.r != rgb.b) return false;
			else if (rgb.g != rgb.b) return false;
		}
	}

	return true;
}

void PNGLoader::LibPNGReadPalette(DecompressedImage & d)
{
	png_colorp palette;
	int num_palette;
	png_get_PLTE(pngPtr, infoPtr, &palette, &num_palette);

	png_bytep trans;
	int num_trans;
	png_color_16p col;
	png_uint_32 retTRNS = png_get_tRNS(pngPtr, infoPtr, &trans, &num_trans, &col);

	if (retTRNS == PNG_INFO_tRNS)
	{
		d.grayScalePallete = true;
		d.palette.reserve(num_palette);
		for (int i = 0; i < num_palette; i++)
		{
			RGBA rgba;
			rgba.r = palette[i].red;
			rgba.g = palette[i].green;
			rgba.b = palette[i].blue;

			if (rgba.r != rgba.g) d.grayScalePallete = false;
			else if (rgba.r != rgba.b) d.grayScalePallete = false;
			else if (rgba.g != rgba.b) d.grayScalePallete = false;

			rgba.a = (i < num_trans) ? trans[i] : 255;
            d.palette.push_back(std::move(rgba));
		}
	}
	else
	{
		d.grayScalePallete = true;
		d.palette.reserve(num_palette);
		for (int i = 0; i < num_palette; i++)
		{
			RGBA rgb;
			rgb.r = palette[i].red;
			rgb.g = palette[i].green;
			rgb.b = palette[i].blue;

			if (rgb.r != rgb.g) d.grayScalePallete = false;
			else if (rgb.r != rgb.b) d.grayScalePallete = false;
			else if (rgb.g != rgb.b) d.grayScalePallete = false;

			rgb.a = 255;
            d.palette.push_back(std::move(rgb));
		}
	}
}

void PNGLoader::LibPNGReadData(DecompressedImage & dec)
{
	
	//Here's one of the pointers we've defined in the error handler section:
	//Array of row pointers. One for every row.
	rowPtrs = new png_bytep[dec.h];

	//Alocate a buffer with enough space.
	//(Don't use the stack, these blocks get big easilly)
	//This pointer was also defined in the error handling section, so we can clean it up on error.
	//data = new char[dec.w * dec.h * bitdepth * channels / 8];
		
	//This is the length in bytes, of one row.
	//const unsigned int stride = dec.w * dec.bitDepth * dec.channelsCount / 8;
	const png_size_t stride = png_get_rowbytes(pngPtr, infoPtr);


	dec.data.resize(dec.h * stride);

	//A little for-loop here to set all the row pointers to the starting
	//Adresses for every row in the buffer
	
	for (unsigned i = 0; i < dec.h; i++)
	{
		//Set the pointer to the data pointer + i times the row stride.
		
		//Notice that the row order is reversed with q.
		//This is how at least OpenGL expects it,
		//and how many other image loaders present the data.
		//png_uint_32 q = (dec.h - i - 1) * stride;

		//original order
		rowPtrs[i] = (png_bytep)dec.data.data() + i * stride;
	}
		
	//And here it is! The actuall reading of the image!
	//Read the imagedata and write it to the adresses pointed to
	//by rowptrs (in other words: our image databuffer)
	png_read_image(pngPtr, rowPtrs);
}

#endif
