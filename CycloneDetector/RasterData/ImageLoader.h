#ifndef IMAGE_LOADER_H
#define IMAGE_LOADER_H

struct IFile;

#include <array>

#include "../Utils/IDataLoader.h"

#include "../Compression/PNGLoader.h"

class ImageLoader : public MyUtils::IDataLoader
{
public:
	enum class FILE_TYPE { PNG = 0, UNKNOWN = 2 };
	typedef enum CHANNEL { RED = 0, GREEN = 1, BLUE = 2, ALPHA = 3, NONE = 4 } CHANNEL;
	
					
	int hasAdditionalData;
	bool corrupted;	
	
	ImageLoader(const char * textureName);
	~ImageLoader() = default;
		
	bool IsOptionalAlphaEnabled() const;

	void EnableFilesJoin(bool val);	
	void EnableChannelMapping(bool val);
	void EnableOptionalAlpha(bool val);
	void SetChannelMapping(size_t fileIndex, CHANNEL input, CHANNEL output);	
	void SetFileType(size_t fileIndex, FILE_TYPE type);

	int GetSumChannelsCount() const;

	virtual void Start() override;
	
	
	
private:
	bool joinFiles;
	bool channelMapping;
	bool optionalAlpha;
	bool storeAlpha;

	std::vector<FILE_TYPE> fileTypes;
	std::vector<std::array<char, 4>> outMapping;
	std::vector<int> outputChannelsCount;

	FILE_TYPE GetFileType(IFile * f);

	void LoadPNG(IFile * f, size_t fileIndex);
	
	void JoinAllToOneImage();

	std::vector<uint8_t> Convert1BitTo8Bit(const std::vector<uint8_t> & data, size_t w, size_t h);
	std::vector<uint8_t> Convert4BitTo8Bit(const std::vector<uint8_t> & data, size_t w, size_t h);

	void UnpackPallete1Bit(const PNGLoader::DecompressedImage & dec,
		int outChannelsCount,
		const std::array<char, 4> & mapping,
		std::vector<uint8_t> & target);
	void UnpackPallete2Bit(const PNGLoader::DecompressedImage & dec,
		int outChannelsCount,
		const std::array<char, 4> & mapping,
		std::vector<uint8_t> & target);
	void UnpackPallete4Bit(const PNGLoader::DecompressedImage & dec,
		int outChannelsCount,
		const std::array<char, 4> & mapping,
		std::vector<uint8_t> & target);

	void WriteToTarget(int index, const PNGLoader::RGBA & val, 
		const std::array<char, 4> & mapping,
		std::vector<uint8_t> & target);

	void ColorMapping(size_t fileIndex, size_t w, size_t, int channelsCount, const std::vector<uint8_t> & data, 
		LoadedData & l);
};


#endif
