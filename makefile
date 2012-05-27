
all: projekt movie.mp4

movie.mp4: projekt
	rm -f image*.png movie.mp4
	./projekt
	ffmpeg -r 10 -b 9600 -i image%05d.png movie.mp4
	rm -f image*.png


projekt: projekt.c
	gcc -O3 -g -msse4 -std=gnu99 -Wall -Wextra projekt.c `pkg-config --cflags --libs fftw3 cairo` -o projekt

clean:
	rm -f projekt
	rm -f movie.mp4

