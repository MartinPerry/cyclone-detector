#include "./ImageLoader.h"

//#include "../../Strings/MyStringAnsi.h"

#include "../Compression/PNGLoader.h"

#include "../Utils/FileWrapper.h"
#include "../Utils/Logger.h"

#include "../Macros.h"

//================================================================================================
// ctor
//================================================================================================

/// <summary>
/// Image loader
/// - supposed to run in background thread
/// </summary>
/// <param name="textureName">name of texture that is being loaded</param>
ImageLoader::ImageLoader(const char * textureName) :	
	IDataLoader(textureName), 	
	hasAdditionalData(false),
	corrupted(false),	
	joinFiles(false),
	channelMapping(false),
	optionalAlpha(false),
	storeAlpha(true)
{
	
}

//================================================================================================
// Setters & getters
//================================================================================================

bool ImageLoader::IsOptionalAlphaEnabled() const
{
	return this->optionalAlpha;
}


/// <summary>
/// If we load multiple files, we want to join them to one
/// Example:
/// Load two files 8bit -> join them to a single RG file
/// </summary>
/// <param name="val"></param>
void ImageLoader::EnableFilesJoin(bool val)
{
	this->joinFiles = val;
}

/// <summary>
/// Enable channel mapping
/// If channel mapping is disabled - decompressed output is
/// used as it is 
/// If channel mapping is enabled - decompressed output
/// is mapped to channels specified by SetChannelMapping
/// </summary>
/// <param name="val"></param>
void ImageLoader::EnableChannelMapping(bool val)
{
	this->channelMapping = val;
}

/// <summary>
/// Enable / disable optional alpha
/// If optional alpha is enabled and image is pallete,
/// test if all values are solid (a = 255)
/// If yes, alpha channel is not stored and channels count is decreased
/// </summary>
/// <param name="val"></param>
void ImageLoader::EnableOptionalAlpha(bool val)
{
	this->optionalAlpha = val;
}

/// <summary>
/// Set mapping of channels
/// -> input channel is mapped to output channel
/// E.g. input RGBA => mapped to RG
/// R => R
/// G => NONE
/// B => NONE
/// A => G
/// 
/// Note: If channel mapping has been disabled, calling
/// this method will automatically enable it
/// </summary>
/// <param name="input"></param>
/// <param name="output"></param>
void ImageLoader::SetChannelMapping(size_t fileIndex, CHANNEL input, CHANNEL output)
{
	if (input == CHANNEL::NONE)
	{
		return;
	}

	while (outMapping.size() <= fileIndex)
	{
		outMapping.push_back({ CHANNEL::NONE, CHANNEL::NONE, CHANNEL::NONE, CHANNEL::NONE });
	}

	outMapping[fileIndex][input] = output;
	
	this->channelMapping = true;
}

/// <summary>
/// If we know file type, set it
/// It it is not set, file type is deduced from file itself
/// </summary>
/// <param name="fileIndex"></param>
/// <param name="type"></param>
void ImageLoader::SetFileType(size_t fileIndex, FILE_TYPE type)
{
	while (fileTypes.size() <= fileIndex)
	{
		fileTypes.push_back(FILE_TYPE::UNKNOWN);
	}

	fileTypes[fileIndex] = type;
}

/// <summary>
/// Get number of all channels summed from all files
/// If file join is enabled, this is number of channels in joined file
/// </summary>
/// <returns></returns>
int ImageLoader::GetSumChannelsCount() const
{
	int sum = 0;
	for (size_t i = 0; i < this->files.size(); i++)
	{
		sum += this->outputChannelsCount[i];
	}

	return sum;
}

//================================================================================================
// Processing
//================================================================================================

/// <summary>
/// 
/// </summary>
void ImageLoader::Start()
{	
	if (this->files.empty())
	{
		//no files available
		return;
	}

	if ((this->channelMapping) && (this->outMapping.empty()))
	{
		//channels mapping is enabled
		//but there is no mapping available

        MY_LOG_ERROR("No output mapping count");
        for (FileHandle & fh : this->files)
        {
            if (fh.closeAtFinish)
            {
                SAFE_DELETE(fh.f);
            }
        }
        
		return;
	}
	
	//prepare file types for files that have not filetype pre-set
	while (fileTypes.size() < this->files.size())
	{
		fileTypes.push_back(FILE_TYPE::UNKNOWN);
	}

	
	if ((this->channelMapping) && (this->files.size() != this->outMapping.size()))
	{
		MY_LOG_ERROR("Number of files (%i) is not same as output mapping count (%i)", 
			this->files.size(), this->outMapping.size());
        
        for (FileHandle & fh : this->files)
        {
            if (fh.closeAtFinish)
            {
                SAFE_DELETE(fh.f);
            }
        }
        
		return;
	}

	this->corrupted = false;

	//Calculate number of channels for output files
	//based on mapping
	//If mapping is disabled, result is not important - 
	//result from decompressor is used directly	
	this->outputChannelsCount.clear();

	for (size_t i = 0; i < this->files.size(); i++)
	{
		this->outputChannelsCount.push_back(0);

		if (this->channelMapping)
		{
			const auto & tmp = this->outMapping[i];
			for (int j = 0; j < 4; j++)
			{
				if (tmp[j] != CHANNEL::NONE)
				{
					this->outputChannelsCount[i]++;
				}
			}
		}
	}

	//Decompress all files
	for (size_t i = 0; i < this->files.size(); i++)
	{
		FileHandle & fh = this->files[i];

		FILE_TYPE ft = this->fileTypes[i];
		if (ft == FILE_TYPE::UNKNOWN)
		{
			ft = this->GetFileType(fh.f);
		}

		//check header
		if (ft == FILE_TYPE::PNG)
		{
			this->LoadPNG(fh.f, i);
		}		
		else
		{
			MY_LOG_ERROR("UNKNOWN file format for %s", dataName.c_str());
			this->corrupted = true;
		}
	}

	//Finished decompression
	//Callbacks & cleanup

	if (this->onFinishCallback)
	{
		this->onFinishCallback(this);
	}
    
	
	for (size_t i = 0; i < this->files.size(); i++)
	{
		FileHandle & fh = this->files[i];

        if (fh.closeAtFinish)
        {
            SAFE_DELETE(fh.f);
        }
	}

	if (this->joinFiles)
	{
		this->JoinAllToOneImage();
	}

    
	finished.store(true);
}


void ImageLoader::JoinAllToOneImage()
{
	if (this->files.size() == 1)
	{
		//nothing to join if there was only one input file
		return;
	}

	if (this->loaded.empty())
	{
		return;
	}

	LoadedData ldJoin;
	ldJoin.channelsCount = this->loaded[0].channelsCount;
	ldJoin.w = this->loaded[0].w;
	ldJoin.h = this->loaded[0].h;
    
	//check if join can be done
	for (auto & l : this->loaded)
	{
		if (ldJoin.w != l.w)
		{		
			MY_LOG_ERROR("Unable to join - different width");
			return;
		}
		if (ldJoin.h != l.h)
		{
			MY_LOG_ERROR("Unable to join - different height");
			return;
		}
		if (ldJoin.channelsCount != l.channelsCount)
		{
			MY_LOG_ERROR("Unable to join - different channels count");
			return;
		}
	}

	//actually join the files
	ldJoin.rawData.reserve(this->loaded.size() * this->loaded[0].rawData.size());
	
	for (size_t i = 0; i < this->loaded[0].rawData.size(); i++)
	{
		for (size_t j = 0; j < this->loaded.size(); j++)
		{
            ldJoin.rawData.push_back(std::move(this->loaded[j].rawData[i]));			
		}
	}
	this->loaded.clear();
    this->loaded.push_back(std::move(ldJoin));
}

/// <summary>
/// Get file type from its starting bytes
/// </summary>
/// <returns></returns>
ImageLoader::FILE_TYPE ImageLoader::GetFileType(IFile * f)
{
	uint8_t header[8];
	f->Read(header, sizeof(uint8_t), 8);
	f->Seek(0, SEEK_SET);

	if ((header[0] == 137) && (header[1] == 'P')) // N G \r \n \032 \n
	{
		return FILE_TYPE::PNG;
	}

	return FILE_TYPE::UNKNOWN;
}

//================================================================================================
// Image loading section
//================================================================================================

/// <summary>
/// Read PNG file
/// </summary>
/// <param name="outChannelsCount">number of output channels (1 = R or 2 = RA)</param>
void ImageLoader::LoadPNG(IFile * f, size_t fileIndex)
{
	PNGLoader pngLoad;
    if (this->channelMapping)
    {
		//with channel mapping, keep palette
		//because we build the result image manually
        pngLoad.SetKeepPalette(true);
    }
    else
    {
		//no channel mapping - no need for palette
        pngLoad.SetKeepPalette(false);
    }
	
	auto dec = pngLoad.DecompressFromFile(f);

	if ((dec.w == 0) || (dec.h == 0))
	{
		//file is corrupted
		this->corrupted = true;	
		return;
	}

	if ((this->channelMapping == false) && (dec.bitDepth == 2))
	{		
		//2 bits not supported without pallete !!

		//reset stream to beginning !!!
		//and reload image
		f->Seek(0L, SEEK_SET);
		pngLoad.SetKeepPalette(false);		
		dec = pngLoad.DecompressFromFile(f);
	}

    if ((dec.w == 0) || (dec.h == 0))
    {
        //file is corrupted
        this->corrupted = true;
        return;
    }
	
	
	LoadedData l;	
	l.w = dec.w;
	l.h = dec.h;	

	if (this->channelMapping == false)
	{
		//no channels mapping, return image as it is
		this->outputChannelsCount[fileIndex] = dec.channelsCount;
		l.channelsCount = dec.channelsCount;
		l.rawData = std::move(dec.data);
		this->loaded.push_back(std::move(l));
		return;
	}

	int outChannelsCount = this->outputChannelsCount[fileIndex];
	const std::array<char, 4> & mapping = this->outMapping[fileIndex];
	
	
	if (dec.palette.empty())
	{
		//no palette

		l.channelsCount = outChannelsCount;
		l.rawData.resize(dec.w * dec.h * outChannelsCount, 255);

		if (dec.bitDepth == 1)
		{
			dec.data = this->Convert1BitTo8Bit(dec.data, dec.w, dec.h);
			dec.bitDepth = 8;
		}
		else if (dec.bitDepth == 4)
		{
			dec.data = this->Convert4BitTo8Bit(dec.data, dec.w, dec.h);
			dec.bitDepth = 8;
		}
		this->ColorMapping(fileIndex, dec.w, dec.h, dec.channelsCount, dec.data, l);
	}
	else
	{
		//palette is kept in image
		//unpack palette and do color mapping directly
	
		this->storeAlpha = true;
		if ((this->optionalAlpha) && (outChannelsCount > 1))
		{
			this->storeAlpha = false;

			//alpha channel is optional
			//it can be omited if there is non (all alpha are 255)
			for (auto & v : dec.palette)
			{
				if (v.a != 255)
				{
					this->storeAlpha = true;
					break;
				}
			}

			if (this->storeAlpha == false)
			{
				outChannelsCount--;
				this->outputChannelsCount[fileIndex] = outChannelsCount;
			}
		}

		l.channelsCount = outChannelsCount;
		l.rawData.resize(dec.w * dec.h * outChannelsCount, 255);

		if (dec.bitDepth == 1)
		{	
			this->UnpackPallete1Bit(dec, outChannelsCount, mapping,
				l.rawData);
		}
		else if (dec.bitDepth == 2)
		{
			this->UnpackPallete2Bit(dec, outChannelsCount, mapping,
				l.rawData);

		}
		else if (dec.bitDepth == 4)
		{
			this->UnpackPallete4Bit(dec, outChannelsCount, mapping,
				l.rawData);

		}
		else
		{
			int index = 0;
			for (unsigned int i = 0; i < dec.w * dec.h; i++)
			{
				auto rgba = dec.palette[dec.data[i]];

				WriteToTarget(index, rgba, mapping, l.rawData);

				index += outChannelsCount;								
			}
		}
	}
	

    this->loaded.push_back(std::move(l));
			
}

//================================================================================================
// PNG-based pallete and unpacking manipulation
//================================================================================================

/// <summary>
/// Convert 1bit input data to 8 bit output
/// -> unpack data from 1 bit no pallete PNG
/// </summary>
/// <param name="data"></param>
/// <returns></returns>
std::vector<uint8_t> ImageLoader::Convert1BitTo8Bit(const std::vector<uint8_t> & data, size_t w, size_t h)
{
	
	const uint64_t MASK = 0x00000000000000ff;

	std::vector<uint8_t> unpacked;
	unpacked.resize(w * h);

	size_t index = 0;


	//end size must be multiple of sizeof(uint64_t)
	//rest of data is processed after that
	size_t remain = data.size() % sizeof(uint64_t);
	size_t endSize = data.size();
	endSize -= remain;
	
	//process multiples of sizeof(uint64_t)
	for (size_t i = 0; i < endSize; i += sizeof(uint64_t))
	{				
		const uint64_t val = *(reinterpret_cast<const uint64_t *>(data.data() + i));

		const uint8_t v[sizeof(uint64_t)] = {
				static_cast<uint8_t>(val & MASK),
				static_cast<uint8_t>((val & MASK << 8) >> 8),
				static_cast<uint8_t>((val & MASK << 16) >> 16),
				static_cast<uint8_t>((val & MASK << 24) >> 24),
				static_cast<uint8_t>((val & MASK << 32) >> 32),
				static_cast<uint8_t>((val & MASK << 40) >> 40),
				static_cast<uint8_t>((val & MASK << 48) >> 48),
				static_cast<uint8_t>((val & MASK << 56) >> 56)
		};

		for (size_t j = 0; j < sizeof(uint64_t); j++)
		{
			uint8_t t = v[j];
			//calculate last "column" width
			int endB = static_cast<int>(w - (index % w));
			if (endB > 7)
			{
				endB = 0;
			}
			else
			{
				endB = sizeof(uint64_t) - endB;
			}

			//endian !!
			for (int b = 7; b >= endB; b--)			
			{
				uint32_t bit = (t >> b) & 1;
				unpacked[index] = static_cast<uint8_t>(bit * 255);
				index++;
			}
		}		
	}
	
	for (size_t i = endSize; i < data.size(); i++)
	{
		uint8_t val = data[i];

		//calculate last "column" width
		int endB = static_cast<int>(w - (index % w));
		if (endB > 7)
		{
			endB = 0;
		}
		else
		{
			endB = sizeof(uint64_t) - endB;
		}

		//endian !!
		for (int b = 7; b >= endB; b--)
		{
			uint32_t bit = (val >> b) & 1;
			unpacked[index] = static_cast<uint8_t>(bit * 255);
			index++;
		}
	}



	return unpacked;
	
}


/// <summary>
/// Convert 4bit input data to 8 bit output
/// -> unpack data from 4 bit no pallete PNG
/// </summary>
/// <param name="data"></param>
/// <returns></returns>
std::vector<uint8_t> ImageLoader::Convert4BitTo8Bit(const std::vector<uint8_t> & data, size_t w, size_t h)
{
	const uint64_t MASK = 0x00000000000000ff;

	std::vector<uint8_t> unpacked;
	unpacked.resize(w * h);

	int index = 0;


	//end size must be multiple of sizeof(uint64_t)
	//rest of data is processed after that
	size_t remain = data.size() % sizeof(uint64_t);
	size_t endSize = data.size();
	endSize -= remain;

	//process multiples of sizeof(uint64_t)
	for (size_t i = 0; i < endSize; i += sizeof(uint64_t))
	{		
		const uint64_t val = *(reinterpret_cast<const uint64_t *>(data.data() + i));

		const uint8_t v[sizeof(uint64_t)] = {
				static_cast<uint8_t>(val & MASK),
				static_cast<uint8_t>((val & MASK << 8) >> 8),
				static_cast<uint8_t>((val & MASK << 16) >> 16),
				static_cast<uint8_t>((val & MASK << 24) >> 24),
				static_cast<uint8_t>((val & MASK << 32) >> 32),
				static_cast<uint8_t>((val & MASK << 40) >> 40),
				static_cast<uint8_t>((val & MASK << 48) >> 48),
				static_cast<uint8_t>((val & MASK << 56) >> 56)
		};

		for (size_t j = 0; j < sizeof(uint64_t); j++)
		{
			uint8_t t = v[j];
			//calculate last "column" width
			int endB = static_cast<int>(w - (index % w));
			if (endB > 2)
			{
				uint32_t fourBit1 = (t >> 4) & 0xF;
				uint32_t fourBit2 = (t) & 0x0F;

				unpacked[index] = static_cast<uint8_t>(fourBit1 * 16);
				index++;

				unpacked[index] = static_cast<uint8_t>(fourBit2 * 16);
				index++;

			}
			else
			{
				uint32_t fourBit1 = (t >> 4) & 0xF;

				unpacked[index] = static_cast<uint8_t>(fourBit1 * 16);
				index++;
			}

		}
	}

	for (size_t i = endSize; i < data.size(); i++)
	{
		uint8_t val = data[i];

		//calculate last "column" width
		int endB = static_cast<int>(w - (index % w));
		if (endB > 2)
		{
			uint32_t fourBit1 = (val >> 4) & 0xF;
			uint32_t fourBit2 = (val) & 0x0F;

			unpacked[index] = static_cast<uint8_t>(fourBit1 * 16);
			index++;

			unpacked[index] = static_cast<uint8_t>(fourBit2 * 16);
			index++;
		}
		else
		{
			uint32_t fourBit1 = (val >> 4) & 0xF;

			unpacked[index] = static_cast<uint8_t>(fourBit1 * 16);
			index++;
		}

	}


	return unpacked;
}

/// <summary>
/// Unpack 1bit pallete (2 colors) to output target
/// with color mapping
/// </summary>
/// <param name="dec"></param>
/// <param name="outChannelsCount"></param>
/// <param name="mapping"></param>
/// <param name="target"></param>
void ImageLoader::UnpackPallete1Bit(const PNGLoader::DecompressedImage & dec,
	int outChannelsCount,
	const std::array<char, 4> & mapping,
	std::vector<uint8_t> & target)
{
	int index = 0;

	const uint64_t MASK = 0x00000000000000ff;

	//end size must be multiple of sizeof(uint64_t)
	//rest of data is processed after that
	size_t remain = dec.data.size() % sizeof(uint64_t);
	size_t endSize = dec.data.size();
	endSize -= remain;


	//process multiples of sizeof(uint64_t)
	for (size_t i = 0; i < endSize; i += sizeof(uint64_t))
	{
		const uint64_t val = *(reinterpret_cast<const uint64_t *>(dec.data.data() + i));

		const uint8_t v[sizeof(uint64_t)] = {
			static_cast<uint8_t>(val & MASK),
			static_cast<uint8_t>((val & MASK << 8) >> 8),
			static_cast<uint8_t>((val & MASK << 16) >> 16),
			static_cast<uint8_t>((val & MASK << 24) >> 24),
			static_cast<uint8_t>((val & MASK << 32) >> 32),
			static_cast<uint8_t>((val & MASK << 40) >> 40),
			static_cast<uint8_t>((val & MASK << 48) >> 48),
			static_cast<uint8_t>((val & MASK << 56) >> 56)
		};

		for (size_t j = 0; j < sizeof(uint64_t); j++)
		{
			uint8_t t = v[j];

			//calculate last "column" width
			int endB = (dec.w - (index % dec.w));
			if (endB > 7)
			{
				endB = 0;
			}
			else
			{
				endB = sizeof(uint64_t) - endB;
			}

			//endian !!
			for (int b = 7; b >= endB; b--)
			{
				uint32_t bit = (t >> b) & 1;
				const auto & rgba = dec.palette[bit];


				WriteToTarget(index, rgba, mapping, target);

				index += outChannelsCount;
			}
		}
	}

	//process rest of data
	for (size_t i = endSize; i < dec.data.size(); i++)
	{
		uint8_t val = dec.data[i];

		//calculate last "column" width
		int endB = (dec.w - (index % dec.w));
		if (endB > 7)
		{
			endB = 0;
		}
		else
		{
			endB = sizeof(uint64_t) - endB;
		}

		//endian !!
		for (int b = 7; b >= endB; b--)
		{
			uint32_t bit = (val >> b) & 1;
			const auto & rgba = dec.palette[bit];

			WriteToTarget(index, rgba, mapping, target);

			index += outChannelsCount;
		}
	}
}

/// <summary>
/// Unpack 2bit pallete (4 colors) to output target
/// with color mapping
/// </summary>
/// <param name="dec"></param>
/// <param name="outChannelsCount"></param>
/// <param name="mapping"></param>
/// <param name="target"></param>
void ImageLoader::UnpackPallete2Bit(const PNGLoader::DecompressedImage & dec,
	int outChannelsCount,
	const std::array<char, 4> & mapping,
	std::vector<uint8_t> & target)
{
	int index = 0;	
	const int lastByteCount = dec.w % 4;

	if (lastByteCount == 0)
	{
		//even width - all bytes are used

		size_t count = (dec.w * dec.h) / 4;
		for (size_t i = 0; i < count; i++)
		{
			const uint8_t val = dec.data[i];
			
			WriteToTarget(index, dec.palette[(val & 0xC0) >> 6], mapping, target);
			index += outChannelsCount;

			WriteToTarget(index, dec.palette[(val & 0x30) >> 4], mapping, target);
			index += outChannelsCount;

			WriteToTarget(index, dec.palette[(val & 0x0C) >> 2], mapping, target);
			index += outChannelsCount;

			WriteToTarget(index, dec.palette[(val & 0x03)], mapping, target);
			index += outChannelsCount;
		}
	}
	else
	{
		uint8_t indices[4];

		//width not multiple of 4 - last byte is only used from part
		size_t i = 0;
		for (size_t y = 0; y < dec.h; y++)
		{
			for (size_t x = 0; x < dec.w - lastByteCount; x += 4)
			{
				const uint8_t val = dec.data[i];
				i++;
								
				WriteToTarget(index, dec.palette[(val & 0xC0) >> 6], mapping, target);
				index += outChannelsCount;

				WriteToTarget(index, dec.palette[(val & 0x30) >> 4], mapping, target);
				index += outChannelsCount;

				WriteToTarget(index, dec.palette[(val & 0x0C) >> 2], mapping, target);
				index += outChannelsCount;

				WriteToTarget(index, dec.palette[(val & 0x03)], mapping, target);
				index += outChannelsCount;				
			}

			const uint8_t val = dec.data[i];
			i++;
			
			indices[0] = (val & 0xC0) >> 6;
			indices[1] = (val & 0x30) >> 4;
			indices[2] = (val & 0x0C) >> 2;
			indices[3] = (val & 0x03);			

			for (int b = 0; b < lastByteCount; b++)
			{				
				WriteToTarget(index, dec.palette[indices[b]], mapping, target);
				index += outChannelsCount;
			}
		}
	}

}


/// <summary>
/// Unpack 4bit pallete (16 colors) to output target
/// with color mapping
/// </summary>
/// <param name="dec"></param>
/// <param name="outChannelsCount"></param>
/// <param name="mapping"></param>
/// <param name="target"></param>
void ImageLoader::UnpackPallete4Bit(const PNGLoader::DecompressedImage & dec,
	int outChannelsCount,
	const std::array<char, 4> & mapping,
	std::vector<uint8_t> & target)
{
	int index = 0;
	
	if (dec.w % 2 == 0)
	{
		//even width - all bytes are used

		size_t count = (dec.w * dec.h) / 2;
		for (size_t i = 0; i < count; i++)
		{
			const uint8_t val = dec.data[i];
			
			WriteToTarget(index, dec.palette[(val & 0xF0) >> 4], mapping, target);
			index += outChannelsCount;
			
			WriteToTarget(index, dec.palette[(val & 0x0F)], mapping, target);
			index += outChannelsCount;
			
		}
	}
	else
	{
		//odd width - last byte is only used from half
		size_t i = 0;
		for (size_t y = 0; y < dec.h; y++)
		{
			for (size_t x = 0; x < dec.w - 1; x += 2)
			{
				const uint8_t val = dec.data[i];				
				i++;

				WriteToTarget(index, dec.palette[(val & 0xF0) >> 4], mapping, target);
				index += outChannelsCount;

				WriteToTarget(index, dec.palette[(val & 0x0F)], mapping, target);
				index += outChannelsCount;				
			}

			const uint8_t val = dec.data[i];
			i++;

			WriteToTarget(index, dec.palette[(val & 0xF0) >> 4], mapping, target);
			index += outChannelsCount;						
		}
	}

}

void ImageLoader::WriteToTarget(int index, const PNGLoader::RGBA & val,
	const std::array<char, 4> & mapping,
	std::vector<uint8_t> & target)
{	
	if (mapping[0] != CHANNEL::NONE)
	{
		target[index + mapping[0]] = val.r;
	}

	if (mapping[1] != CHANNEL::NONE)
	{
		target[index + mapping[1]] = val.g;
	}

	if (mapping[2] != CHANNEL::NONE)
	{
		target[index + mapping[2]] = val.b;
	}

	if ((mapping[3] != CHANNEL::NONE) && (this->storeAlpha))
	{
		target[index + mapping[3]] = val.a;
	}
}

//================================================================================================
// Final color mapping
//================================================================================================

/// <summary>
/// copy "raw" data loaded from file to our output buffer using defined color mapping
/// </summary>
/// <param name="fileIndex"></param>
/// <param name="w"></param>
/// <param name="h"></param>
/// <param name="channelsCount"></param>
/// <param name="data"></param>
/// <param name="l"></param>
void ImageLoader::ColorMapping(size_t fileIndex, size_t w, size_t h, int channelsCount,
	const std::vector<uint8_t> & data, LoadedData & l)
{

	int outChannelsCount = this->outputChannelsCount[fileIndex];
	const auto & mapping = this->outMapping[fileIndex];

	size_t imgSize = w * h * channelsCount;


	int index = 0;

	if (channelsCount == 4)  //PNG_COLOR_TYPE_RGB_ALPHA
	{
		for (size_t i = 0; i < imgSize; i += channelsCount)
		{
			for (int j = 0; j < 4; j++)
			{
				if (mapping[j] != CHANNEL::NONE)
				{
					l.rawData[index + mapping[j]] = data[i + j];
				}
			}

			index += outChannelsCount;
		}
	}
	else if (channelsCount == 2) //PNG_COLOR_TYPE_GRAY_ALPHA 
	{
		for (size_t i = 0; i < imgSize; i += channelsCount)
		{

			if (mapping[0] != CHANNEL::NONE) //RED
			{
				l.rawData[index + mapping[0]] = data[i + 0];
			}
			if (mapping[3] != CHANNEL::NONE) //ALPHA
			{
				l.rawData[index + mapping[3]] = data[i + 1];
			}

			//added channels
			/*
			if (mapping[1] != CHANNEL::NONE) //GREEN
			{
				l.rawData[index + mapping[1]] = 255;
			}
			if (mapping[2] != CHANNEL::NONE) //BLUE
			{
				l.rawData[index + mapping[2]] = 255;
			}
			*/

			index += outChannelsCount;
		}
	}
	else if (channelsCount == 3) //PNG_COLOR_TYPE_RGB
	{
		for (size_t i = 0; i < imgSize; i += channelsCount)
		{

			if (mapping[0] != CHANNEL::NONE) //RED
			{
				l.rawData[index + mapping[0]] = data[i + 0];
			}

			if (mapping[1] != CHANNEL::NONE) //GREEN
			{
				l.rawData[index + mapping[1]] = data[i + 1];
			}
			if (mapping[2] != CHANNEL::NONE) //BLUE
			{
				l.rawData[index + mapping[2]] = data[i + 2];
			}

			//added channels
			/*
			if (mapping[3] != CHANNEL::NONE) //ALPHA
			{
				l.rawData[index + mapping[3]] = 255;
			}
			*/

			index += outChannelsCount;
		}
	}
	else if (channelsCount == 1) //PNG_COLOR_TYPE_GRAY
	{
		//data.size() == imgSize
		
		const uint64_t MASK = 0x00000000000000ff;

		//end size must be multiple of sizeof(uint64_t)
		//rest of data is processed after that
		size_t remain = data.size() % sizeof(uint64_t);
		size_t endSize = data.size();
		endSize -= remain;
		
		for (size_t i = 0; i < endSize; i += sizeof(uint64_t))
		{

			//uint64_t val = 0;
			//memcpy(&val, data.data() + i, sizeof(uint64_t));

			/*
			uint8_t v[sizeof(uint64_t)];
			v[0] = static_cast<uint8_t>(val & MASK);
			v[1] = static_cast<uint8_t>((val & MASK << 8) >> 8);
			v[2] = static_cast<uint8_t>((val & MASK << 16) >> 16);
			v[3] = static_cast<uint8_t>((val & MASK << 24) >> 24);
			v[4] = static_cast<uint8_t>((val & MASK << 32) >> 32);
			v[5] = static_cast<uint8_t>((val & MASK << 40) >> 40);
			v[6] = static_cast<uint8_t>((val & MASK << 48) >> 48);
			v[7] = static_cast<uint8_t>((val & MASK << 56) >> 56);
			*/

			const uint64_t val = *(reinterpret_cast<const uint64_t *>(data.data() + i));

			const uint8_t v[sizeof(uint64_t)] = {
				static_cast<uint8_t>(val & MASK),
				static_cast<uint8_t>((val & MASK << 8) >> 8),
				static_cast<uint8_t>((val & MASK << 16) >> 16),
				static_cast<uint8_t>((val & MASK << 24) >> 24),
				static_cast<uint8_t>((val & MASK << 32) >> 32),
				static_cast<uint8_t>((val & MASK << 40) >> 40),
				static_cast<uint8_t>((val & MASK << 48) >> 48),
				static_cast<uint8_t>((val & MASK << 56) >> 56)
			};
		
			for (size_t j = 0; j < sizeof(uint64_t); j++)
			{
				if (mapping[0] != CHANNEL::NONE) //RED
				{
					//l.rawData[index + mapping[0]] = data[i];
					l.rawData[index + mapping[0]] = v[j];
				}

				index += outChannelsCount;
			}
		}
		
		
		for (size_t i = endSize; i < data.size(); i++)
		{
			if (mapping[0] != CHANNEL::NONE) //RED
			{
				l.rawData[index + mapping[0]] = data[i];
			}
			index += outChannelsCount;
		}

	}
}
