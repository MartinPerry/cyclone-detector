#ifndef FILE_WRAPPER_H
#define FILE_WRAPPER_H

#include <cstdio>

struct IFile
{
	virtual ~IFile() = default;

	virtual size_t GetSize() const = 0;
	virtual size_t Read(void *  buffer, size_t elementSize, size_t elementCount) = 0;
	virtual void Seek(long  offset, int origin) = 0;
	virtual void Close() = 0;
	virtual void * GetRawFilePtr() = 0;


	size_t ReadAll(void ** buffer);
};

struct RawFile : public IFile
{
	RawFile(const char * path, const char * mode = "rb");
	RawFile(FILE * fp, size_t size = 0);
	virtual ~RawFile();


	size_t GetSize() const override;
	size_t Read(void * buffer, size_t elementSize, size_t elementCount) override;
	void Seek(long  offset, int origin) override;
	void Close() override;
	void * GetRawFilePtr() override;

protected:
	mutable size_t size;
	FILE * fp;
};

#endif

