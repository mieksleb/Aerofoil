CC = gcc
CFLAGS = -Wall -Wextra -ggdb -O0 -DDEBUG -v
LIBS += -lm

SRC = vector.c utils.c math_utils.c main.c
OBJ = $(SRC:.c=.o)
EXECUTABLE = aerofoil

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(EXECUTABLE)
