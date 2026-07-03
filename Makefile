CC      = gcc
CFLAGS  = -O3 -Wall
LDLIBS  = -lm
TARGET  = repeatnet
SRC     = repeatnet.c

.PHONY: all debug clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET) $(LDLIBS)

# Sanitizer build for debugging memory/UB issues
debug: $(SRC)
	$(CC) -O1 -g -Wall -fsanitize=address,undefined $(SRC) -o $(TARGET) $(LDLIBS)

clean:
	rm -f $(TARGET)
