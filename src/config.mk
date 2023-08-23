CC = gcc
CFLAGS = -Wall -Wextra -g -DDEBUG

LIBS += -lm

SRC = math_utils.c utils.c vector.c main.c
OBJ = $(SRC:.c=.o)
EXECUTABLE = aerofoil

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(EXECUTABLE)
