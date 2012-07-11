#include "sdl_draw.h"

int main(int argc, char *argv[])
{
	int w = 222 * 2;
	int h = 222;
	int lebuf[w][h];

	sdl_open(w, h);

	for (int c = 0; c < 22; c++) {
		for (int x = 0; x < (w * h); x++)
			lebuf[0][x] = c;
		sdl_plot(w, h, lebuf);
		if (sdl_waitkey())
			return 1;
	}

	sdl_close();
	return 0;
}
