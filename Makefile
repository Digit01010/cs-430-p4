all: raycast.c
	gcc raycast.c -o raycast

clean:
	rm -rf raycast *~

