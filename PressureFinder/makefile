CC = /opt/rh/devtoolset-7/root/usr/bin/g++

INCLUDE_DIRS = -I/usr/include/

LIB_DIRS = -L./ -L/usr/lib64/ -L/usr/lib64/mysql/           

LIBS = -lm -lpthread -lstdc++ -ldl

LINKER_OPTS = -Wall -O2

CFLAGS = -g -mavx -std=c++17 $(INCLUDE_DIRS)

BIN = pressure_finder

#=========================================================================================
COMPRESSION_FILES = $(wildcard ./Compression/*.cpp)
COMPRESSION_3RDPARTY_FILES = $(wildcard ./Compression/3rdParty/*.cpp)

GRAPHICS_FILES = $(wildcard ./Graphics/*.cpp)
GRAPHICS_2D_FILES = $(wildcard ./Graphics/2d/*.cpp)

RASTER_DATA_FILES = $(wildcard ./RasterData/*.cpp)

ADT_FILES = $(wildcard ./ADT/*.cpp)
MATH_FILES = $(wildcard ./Math/*.cpp)
SIMD_FILES = $(wildcard ./SIMD/*.cpp)
UTILS_FILES = $(wildcard ./Utils/*.cpp)
PROJECTION_FILES = $(wildcard ./MapProjections/*.cpp)
ROOT_FILES = $(wildcard ./*.cpp)

#=========================================================================================

OBJ = $(addprefix ,$(notdir $(ADT_FILES:.cpp=.obj))) \
	  $(addprefix ,$(notdir $(COMPRESSION_FILES:.cpp=.obj))) \	  
	  $(addprefix ,$(notdir $(COMPRESSION_3RDPARTY_FILES:.cpp=.obj))) \
	  $(addprefix ,$(notdir $(GRAPHICS_FILES:.cpp=.obj))) \
	  $(addprefix ,$(notdir $(GRAPHICS_2D_FILES:.cpp=.obj))) \	  	  
	  $(addprefix ,$(notdir $(MATH_FILES:.cpp=.obj))) \
	  $(addprefix ,$(notdir $(RASTER_DATA_FILES:.cpp=.obj))) \	  	  
	  $(addprefix ,$(notdir $(SIMD_FILES:.cpp=.obj))) \
	  $(addprefix ,$(notdir $(UTILS_FILES:.cpp=.obj))) \
	  $(addprefix ,$(notdir $(PROJECTION_FILES:.cpp=.obj))) \
	  $(addprefix ,$(notdir $(ROOT_FILES:.cpp=.obj)))

#=========================================================================================

.PHONY: build
.PHONY: clean

#-----------------------------------------

clean:
	rm -f *.obj
	rm -f *.o


#-----------------------------------------
build: ${BIN}


$(BIN): $(OBJ)	
	$(CC) $(LINKER_OPTS) $^ $(LIB_DIRS) $(LIBS) -o $@
#-----------------------------------------

%.obj: ./ADT/%.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

%.obj: ./Compression/%.cpp
	$(CC) $(CFLAGS) -c -o $@ $<
	
%.obj: ./Compression/3rdParty/%.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

%.obj: ./Graphics/%.cpp
	$(CC) $(CFLAGS) -c -o $@ $< 

%.obj: ./Graphics/2d/%.cpp
	$(CC) $(CFLAGS) -c -o $@ $< 

%.obj: ./Math/%.cpp
	$(CC) $(CFLAGS) -c -o $@ $< 

%.obj: ./RasterData/%.cpp
	$(CC) $(CFLAGS) -c -o $@ $< 

%.obj: ./SIMD/%.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

%.obj: ./Utils/%.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

%.obj: ./MapProjections/%.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

%.obj: ./%.cpp
	$(CC) $(CFLAGS) -c -o $@ $<
