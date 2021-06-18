#include "./FileWrapper.h"

#include <cstdlib>

#include "../Macros.h"

size_t IFile::ReadAll(void ** buffer)
{
	size_t fs = this->GetSize();
	*buffer = malloc(fs);

	return this->Read(*buffer, sizeof(char), fs);
}

//==========================================================================

RawFile::RawFile(const char * path, const char * mode) :
	size(0),
	fp(nullptr)
{
	my_fopen(&fp, path, mode);
}

RawFile::RawFile(FILE * fp, size_t size) :
	size(size),
	fp(fp)
{
}


RawFile::~RawFile()
{
	this->Close();
}

size_t RawFile::GetSize() const
{
	if (this->size == 0)
	{
		fseek(fp, 0L, SEEK_END);
		this->size = static_cast<size_t>(ftell(fp));
		fseek(fp, 0L, SEEK_SET);
	}

	return this->size;
}

size_t RawFile::Read(void * buffer, size_t elementSize, size_t elementCount)
{
	return fread(buffer, elementSize, elementCount, fp);
}

void RawFile::Seek(long  offset, int origin)
{
	fseek(fp, offset, origin);
}

void * RawFile::GetRawFilePtr()
{
	return this->fp;
}

void RawFile::Close()
{
	if (fp)
	{
		fclose(fp);
		fp = nullptr;
	}
}
