all:
	@gcc \
		-Wall -Wextra \
		-O3 -funroll-loops \
		-g \
		rng.c \
		u512.S fp.S \
		mont.c \
		csidh.c \
		main.c \
		-o main

bench:
	@gcc \
		-Wall -Wextra \
		-O3 -funroll-loops \
		-g -pg \
		rng.c \
		u512.S fp.S \
		mont.c \
		csidh.c \
		bench.c \
		-o main


debug:
	gcc \
		-Wall -Wextra \
		-g \
		rng.c \
		u512.S fp.S \
		mont.c \
		csidh.c \
		main.c \
		-o main

clean:
	rm -f main

