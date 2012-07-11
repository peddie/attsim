#include <SDL/SDL.h>
#include "sdl_draw.h"

void HSVtoRGB(int *r, int *g, int *b, int h, int s, int v);
void setpixel(SDL_Surface * screen, int x, int y, Uint8 r, Uint8 g, Uint8 b);

static SDL_Surface *screen;

void HSVtoRGB(int *r, int *g, int *b, int h, int s, int v)
{
	int f;
	long p, q, t;

	if (s == 0) {
		*r = *g = *b = v;
		return;
	}

	f = ((h % 60) * 255) / 60;
	h /= 60;
	p = (v * (256 - s)) / 256;
	q = (v * (256 - (s * f) / 256)) / 256;
	t = (v * (256 - (s * (256 - f)) / 256)) / 256;

	switch (h) {
	case 0:
		*r = v;
		*g = t;
		*b = p;
		break;
	case 1:
		*r = q;
		*g = v;
		*b = p;
		break;
	case 2:
		*r = p;
		*g = v;
		*b = t;
		break;
	case 3:
		*r = p;
		*g = q;
		*b = v;
		break;
	case 4:
		*r = t;
		*g = p;
		*b = v;
		break;
	default:
		*r = v;
		*g = p;
		*b = q;
		break;
	}
}

void setpixel(SDL_Surface * screen, int x, int y, Uint8 r, Uint8 g, Uint8 b)
{
	Uint32 *pixmem32;
	Uint32 colour;

	colour = SDL_MapRGB(screen->format, r, g, b);

	pixmem32 = (Uint32 *) screen->pixels + y * screen->w + x;
	*pixmem32 = colour;
}

int sdl_plot(int w, int h, const int *pix)
{

	if (SDL_MUSTLOCK(screen)) {
		if (SDL_LockSurface(screen) < 0)
			return 2;
	}

	for (int y = 0; y < h; y++) {
		for (int x = 0; x < w; x++) {
			int c = pix[x + y * w] % 255;
			int r, g, b;
			if (c == 0)
				r = g = b = 0;
			else {
				c--;
				HSVtoRGB(&r, &g, &b, (c % 10) * 255 / 9,
					 255 - c * 5, 255);
			}
			setpixel(screen, x, y, r, g, b);
		}
	}

	if (SDL_MUSTLOCK(screen))
		SDL_UnlockSurface(screen);

	SDL_Flip(screen);

	return 0;

}

int sdl_open(int w, int h)
{

	if (SDL_Init(SDL_INIT_VIDEO) < 0)
		return 1;

	if (!(screen = SDL_SetVideoMode(w, h, 32, SDL_HWSURFACE))) {
		SDL_Quit();
		return 1;
	}
	return 0;
}

void sdl_close()
{

	SDL_Quit();
}

int sdl_poll_key() {
  SDL_Delay(5); // 5 ms delay to allow some other tasks to run
	SDL_Event event;
  while (SDL_PollEvent(&event)) { // Retrieve all pending events
	  switch (event.type) {
		case SDL_QUIT:
			return 2;
			break;
		case SDL_KEYDOWN:
			return 1;
			break;
    default:
      // Not one we want, discard.
      break;
		}
	}
  return 0;
}


int sdl_waitkey()
{
  int a=0;
  while(!(a = sdl_poll_key()));
  if (a == 2)
    return 2;
  return 0;
}
