CC			=/usr/local/bin/g++-10
CFLAGS		=-c -Wall -O4 -fopenmp -std=c++17 -I./
LDFLAGS		=-fopenmp -std=c++17

SOURCES		=./testAIFMM.cpp
OBJECTS		=$(SOURCES:.cpp=.o)

# DTYPE_FLAG  = -DUSE_real # 1/R
DTYPE_FLAG  = -DUSE_Hankel # Hankel

EXECUTABLE	=./aifmm

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
		$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
		$(CC) $(CFLAGS) $(SOLVER) $(DTYPE_FLAG) $< -o $@

clean:
	rm a.out aifmm *.o
