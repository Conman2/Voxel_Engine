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
ifeq ($(detected_OS), MSYS)
	COMPILER_FLAGS += -Ofast
	LINKER_FLAGS += -lmingw32 -lSDL2main -lSDL2 -lOpenGL32 -lglew32
else ifeq ($(detected_OS), Linux)
	LINKER_FLAGS += -lSDL2 -lOpenGL32 -lglew32 
else
	$(error Unsupported OS)
endif

#OBJS specifies which files to compile as part of the project
OBJS = Game.cpp
#CC specifies which compiler we're using
CC = g++
# Because of richards great code
COMPILER_FLAGS += -g -fpermissive
#OBJ_NAME specifies the name of our exectuable
OBJ_NAME = Game

#This is the target that compiles our executable
all : $(OBJS)
	$(CC) $(OBJS) $(INCLUDE_PATHS) $(LIBRARY_PATHS) $(COMPILER_FLAGS) $(LINKER_FLAGS) -o $(OBJ_NAME)

