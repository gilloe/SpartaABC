CC = g++
 
# compiler flags:
#  -g     - this flag adds debugging information to the executable file
#  -Wall  - this flag is used to turn on most compiler warnings
CFLAGS  = -std=c++11 -g

# The build target 
TARGET = SpartaABC



all: $(TARGET)

$(TARGET): src/$(TARGET).cpp
	$(CC) -o $(TARGET) src/*.cpp $(CFLAGS)

clean:
	$(RM) $(TARGET)
