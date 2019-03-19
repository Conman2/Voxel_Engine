# Find OS 
ifeq '$(findstring ;,$(PATH))' ';'
	detected_OS := Windows
else
	detected_OS := $(shell uname 2>/dev/null || echo Unknown)
	detected_OS := $(patsubst CYGWIN%,Cygwin,$(detected_OS))
	detected_OS := $(patsubst MSYS%,MSYS,$(detected_OS))
	detected_OS := $(patsubst MINGW%,MSYS,$(detected_OS))
endif

# Print OS
$(info Compiling on $(detected_OS))

# Set flags depending on OS
ifeq ($(detected_OS),MSYS)
	#$(info hi richard)
	#COMPILER_FLAGS specifies the additional compilation options we're using
	# -w suppresses all warnings
	# -Wl,-subsystem,windows gets rid of the console window
	COMPILER_FLAGS += -w -Wl,-subsystem,windows
	#INCLUDE_PATHS specifies the additional include paths we'll need
	INCLUDE_PATHS += -IC:\mingw_dev_lib\include\SDL2
	#LIBRARY_PATHS specifies the additional library paths we'll need
	LIBRARY_PATHS += -LC:\mingw_dev_lib\lib
	#LINKER_FLAGS specifies the libraries we're linking against
	LINKER_FLAGS += -lmingw32 -lSDL2main -lSDL2
else ifeq ($(detected_OS),Linux)
	LINKER_FLAGS += -lSDL2
else
	$(error Unsupported OS)
endif

#OBJS specifies which files to compile as part of the project
OBJS = 3D_Engine.cpp
#CC specifies which compiler we're using
CC = g++
# Because of richards shitty code
#COMPILER_FLAGS += -fpermissive
#OBJ_NAME specifies the name of our exectuable
OBJ_NAME = Engine

#This is the target that compiles our executable
all : $(OBJS)
	$(CC) $(OBJS) $(INCLUDE_PATHS) $(LIBRARY_PATHS) $(COMPILER_FLAGS) $(LINKER_FLAGS) -o $(OBJ_NAME)

