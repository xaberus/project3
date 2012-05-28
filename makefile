
all: projekt

movie.mp4: .projekt4
	rm -f movie.mp4
	ffmpeg -r 10 -b 9600 -i image%05d.png movie.mp4
	rm -f image*.png

.projekt4: projekt4
	./projekt4
	touch .projekt4

plotc.txt plotck.txt: projekt3
	./projekt3

plotc.png plotck.png: plotck.txt plots.gp
	gnuplot plots.gp


projekt3: projekt.c
	gcc -DP3 -O3 -g -msse4 -std=gnu99 -Wall -Wextra projekt.c `pkg-config --cflags --libs fftw3 cairo` -o projekt3

projekt4: projekt.c
	gcc -DP4 -O3 -g -msse4 -std=gnu99 -Wall -Wextra projekt.c `pkg-config --cflags --libs fftw3 cairo` -o projekt4

clean:
	rm -f projekt3
	rm -f projekt4
	rm -f movie.mp4
	rm -f image*.png
	rm -f plotck.txt

