CC			=/usr/local/bin/g++-10
# CFLAGS		=-c -Wall -fopenmp -std=c++17 #-I/usr/local/include #-I/usr/local/opt/icu4c/include #-ffast-math -Ofast -ffinite-math-only
CFLAGS		=-c -Wall -O4 -fopenmp -std=c++17 -I./#-I/Users/vaishnavi/Documents/GitHub/PhD_Vaishnavi/AIFMM_LS#-I/usr/local/include #-I/usr/local/opt/icu4c/include #-ffast-math -Ofast -ffinite-math-only
LDFLAGS		=-fopenmp -std=c++17 #-I/usr/local/include #-L/usr/local/opt/icu4c/lib #-Ofast

SOLVER		=-DAIFMM_Direct
# SOLVER		=-DGMRES
# SOLVER		=-DPRECOND_AIFMM
SOURCES		=./preconditioner.cpp
OBJECTS		=$(SOURCES:.cpp=.o)

# DTYPE_FLAG  = -DUSE_integralEqn # 1/R
# LSTYPE_FLAG  = -DnoToLS

# DTYPE_FLAG  = -DUSE_real # 1/R
# LSTYPE_FLAG  = -DnoToLS

DTYPE_FLAG  = -DUSE_Hankel # Hankel
LSTYPE_FLAG  = -DnoToLS

# DTYPE_FLAG  = -DUSE_LS # Lippmann-Schwinger
# LSTYPE_FLAG  = -DyestoLS

EXECUTABLE	=./convergence

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
		$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
		$(CC) $(CFLAGS) $(SOLVER) $(DTYPE_FLAG) $(LSTYPE_FLAG) $< -o $@

clean:
	rm a.out convergence *.o
