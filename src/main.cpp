#include <SDL.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>

uint32_t to_argb(uint8_t a, uint8_t r, uint8_t g, uint8_t b) {
  uint32_t A = a << 24;
  uint32_t R = r << 16;
  uint32_t G = g << 8;
  uint32_t B = b;

  return A + R + G + B;
}

uint32_t* cool_image(const int tx_w, const int tx_h) {
  uint32_t* bitmap = new uint32_t[tx_w * tx_h];
  
  for(int y = 0; y < tx_w; y++) {
    for(int x = 0; x < tx_h; x++) {
      // White edge
      auto shift = std::min({x, y, (tx_w-1) - x, (tx_h-1) - y, 31});
      int edge   = 0x1ff >> shift;

      // Darker circle in the center
      auto norm_x = 2 * (x / (float)tx_w) - 1.f;
      auto norm_y = 2 * (y / (float)tx_h) - 1.f;
      auto norm_r = std::hypot(norm_x, norm_y);

      auto max_darken = 64;
      auto r_reduce   = 2;
      auto circle = (int)(0xff - 0xff * r_reduce * norm_r);

      int darken = std::clamp(circle, 0, max_darken);

      // Cool pattern
      int r = (~x | ~y) & 0xff;
      int g = ( x | ~y) & 0xff;
      int b = (~x |  y) & 0xff;

      // Putting it all together
      int R = std::clamp(edge, r - darken, 0xff);
      int G = std::clamp(edge, g - darken, 0xff);
      int B = std::clamp(edge, b - darken, 0xff);
      
      bitmap[y*tx_w+x] = to_argb(0, R, G, B);
    }
  }

  return bitmap;
}

struct Camera {
  double position[3];
};

int main (int argc, char* argv[]) {
  // const int W = 1280/2, H = 720/2 ;
  const int W = 256, H = 256;

  SDL_Window* window = SDL_CreateWindow(
    "Software Renderer", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
    W*2, H*2, SDL_WINDOW_RESIZABLE
  );
  SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, 0);
  SDL_Texture* texture = SDL_CreateTexture(
    renderer, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STREAMING, W, H
  );

  const int tx_w = 256, tx_h = 256;
  uint32_t* bitmap = cool_image(tx_w, tx_h);

  // Main loop
  uint64_t curr_time = SDL_GetPerformanceCounter();
  uint64_t last_time = 0;
  uint64_t perf_freq = SDL_GetPerformanceFrequency();
  double dt_ms = 0.0;
  double dt    = 0.0;

  double point[2]    = { 0.0, 0.0 };
  double velocity[2] = { 10, 20 };

  std::map<int, bool> keys;
  while(!keys[SDLK_ESCAPE]) {
    last_time = curr_time;
    curr_time = SDL_GetPerformanceCounter();
    perf_freq = SDL_GetPerformanceFrequency();
    dt_ms = (double)((curr_time - last_time) * 1000 / (double)perf_freq);
    dt    = dt_ms / 1000.0;
    std::cout << dt << std::endl;

    SDL_Event ev;
    while(SDL_PollEvent(&ev)) {
      switch(ev.type) {
        case SDL_QUIT:    keys[SDLK_ESCAPE] = true; break;
        case SDL_KEYDOWN: keys[ev.key.keysym.sym] = true; break;
        case SDL_KEYUP:   keys[ev.key.keysym.sym] = false; break;
      }
    }

    uint32_t pixels[W*H] = { 0 };

    for (int i = 0, val = 0; i < H; i++) {
      for (int j = 0; j < W; j++) {
        if (i < tx_h && j < tx_w) pixels[W*i + j] = bitmap[tx_w*i + j];
        
      }
    }


    int loc = (int)point[1]*W + (int)point[0];
    // int loc = (int)point[0];
    std::cout << point[0] << ", " << point[1] << ": " << loc << std::endl;
    pixels[loc] = 0x00;

    point[0] = fmod(dt * velocity[0] + point[0], W);
    if (point[0] < 0) point[0] += W;
    point[1] = fmod(dt * velocity[1] + point[1], H);
    if (point[1] < 0) point[1] += H;


    SDL_UpdateTexture(texture, nullptr, pixels, 4*W);
    SDL_RenderCopy(renderer, texture, nullptr, nullptr);
    SDL_RenderPresent(renderer);

    SDL_Delay(1000 / 32);
  }

  delete[] bitmap;

  return 0;
}