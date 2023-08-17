CC = gcc
CFLAGS = -Wall -Wextra -g -DDEBUG

SRC = utils.c main.c
OBJ = $(SRC:.c=.o)
EXECUTABLE = aerofoil

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(EXECUTABLE)
