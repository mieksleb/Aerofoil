CC = gcc
CFLAGS = -Wall -Wextra -DDEBUG
LIBS = -lm -luxhw

SRC = vector.c utils.c math_utils.c main.c
OBJ = $(SRC:.c=.o)
EXECUTABLE = aerofoil.exe

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(EXECUTABLE)
