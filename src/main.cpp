#include <SDL.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>

#include "sw_math.hpp"

struct PixelData {
  int w, h;
  uint32_t* px;

  PixelData(std::size_t w, std::size_t h) : w(w), h(h) {
    px = new uint32_t[w * h]();
  }

  ~PixelData() {
    std::cout << "DESTRUCTOR" << std::endl;
    delete[] px;
  }

  inline uint32_t& at(int x, int y) {
    return px[y*w + x];
  }
};

struct Camera {
  math::Vec<3> pos;
  math::Vec<3> rot;
};

uint32_t to_argb(uint8_t a, uint8_t r, uint8_t g, uint8_t b) {
  uint32_t A = a << 24;
  uint32_t R = r << 16;
  uint32_t G = g << 8;
  uint32_t B = b;

  return A + R + G + B;
}

PixelData cool_image() {
  PixelData tx(256, 256);
  
  for(int y = 0; y < tx.w; y++) {
    for(int x = 0; x < tx.h; x++) {
      // White edge
      auto shift = std::min({x, y, (tx.w-1) - x, (tx.h-1) - y, 31});
      int edge   = 0x1ff >> shift;

      // Darker circle in the center
      auto norm_x = 2 * (x / (float)tx.w) - 1.f;
      auto norm_y = 2 * (y / (float)tx.h) - 1.f;
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
      
      tx.at(x, y) = to_argb(0, R, G, B);
    }
  }

  return tx;
}

// Converts from normalized coordinates [-1, 1] to texture coordinates
math::Vec<2> unnormalize_coords(PixelData px, math::Vec<2> v) {
  return { px.w / 2.0 * (v[0] + 1.0), -px.h / 2.0 * (v[1] - 1.0) };
}

// Converts from texture coordinates to normalized coordinates [-1, 1]
math::Vec<2> normalize_coords(PixelData px, math::Vec<2> v) {
  return { v[0] * (2.0 / px.w) - 1.0, v[1] * (2.0 / -px.h) + 1.0 };
}

void draw_scanline(PixelData& pixels, int y, int x_start, int x_end, uint32_t c) {
  std::cout << pixels.px << std::endl;
  y       = std::clamp(y,       0, pixels.h);
  x_start = std::clamp(x_start, 0, pixels.w);
  x_end   = std::clamp(x_end  , 0, pixels.w);
  for (int x = x_start; x <= x_end; x++) {
    // std::cout << x << ", " << y << std::endl;
    pixels.at(x, y) = 0x00ffffff;
  }
}

/**
 * Given PixelData, color and 3 Vec<2> with values between -1.0 and 1.0,
 * rasterize_triangle will draw a triangle on the PixelData
 */
void rasterize_triangle(
  PixelData px, math::Vec<2> a, math::Vec<2> b, math::Vec<2> c, uint32_t color
) {
  a = unnormalize_coords(px, a);
  b = unnormalize_coords(px, b);
  c = unnormalize_coords(px, c);

  if (std::tie(b[1], b[0]) < std::tie(a[1], a[0])) { std::swap(a, b); }
  if (std::tie(c[1], c[0]) < std::tie(a[1], a[0])) { std::swap(a, c); }
  if (std::tie(c[1], c[0]) < std::tie(b[1], b[0])) { std::swap(b, c); }

  if ((int)a[0] == (int)c[0]) return;

  bool middle_on_right = 
    (b[1] - a[1]) * (c[0] - a[0]) < (b[0] - a[0]) * (c[1] - a[1]);
  
  double slopes[2];

  if ((int)a[1] < (int)b[1]) {

  }

  if ((int)b[1] < (int)c[1]) {

  }
}

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

  auto bitmap = cool_image();

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
    std::cout << "Nothing executed yet..." << std::endl;
    last_time = curr_time;
    curr_time = SDL_GetPerformanceCounter();
    perf_freq = SDL_GetPerformanceFrequency();
    dt_ms = (double)((curr_time - last_time) * 1000 / (double)perf_freq);
    dt    = dt_ms / 1000.0;
    // std::cout << dt << std::endl;

    SDL_Event ev;
    while(SDL_PollEvent(&ev)) {
      switch(ev.type) {
        case SDL_QUIT:    keys[SDLK_ESCAPE] = true; break;
        case SDL_KEYDOWN: keys[ev.key.keysym.sym] = true; break;
        case SDL_KEYUP:   keys[ev.key.keysym.sym] = false; break;
      }
    }

    std::cout << "Huh? " << std::flush;
    PixelData pixels(W, H);
    std::cout << "What?" << std::endl;

    for (int i = 0, val = 0; i < H; i++) {
      for (int j = 0; j < W; j++) {
        if (i < bitmap.h && j < bitmap.w) pixels.at(i, j) = bitmap.at(i, j);
        
      }
    }

    // pixels.at(point[0], point[1]) = 0x00;

    // draw_scanline(pixels, point[1],       point[0], point[0] + 2, 0xffffffff);
    // draw_scanline(pixels, point[1] + 1.0, point[0], point[0] + 2, 0xffffffff);
    // std::cout << "BEFORE" << std::endl;
    std::cout << pixels.px << " ";

    // uint32_t color = 0;
    // for (int asdf = 0; asdf < H; asdf++) {
    //   draw_scanline(pixels, asdf, 0, W - 1, 0xff00ff);
    //   color <<= 4;
    //   color += 1;
    // }
  //   auto y = 100, x_start = 100, x_end = 200;
  //     std::cout << pixels.px << std::endl;
  // y       = std::clamp(y,       0, pixels.h);
  // x_start = std::clamp(x_start, 0, pixels.w);
  // x_end   = std::clamp(x_end  , 0, pixels.w);
  // for (int x = x_start; x <= x_end; x++) {
  //   // std::cout << x << ", " << y << std::endl;
  //   pixels.at(x, y) = 0x00ffffff;
  // }

    // std::cout << "HERE" << std::endl;

    // point[0] = fmod(dt * velocity[0] + point[0], W);
    // if (point[0] < 0) point[0] += W;
    // point[1] = fmod(dt * velocity[1] + point[1], H);
    // if (point[1] < 0) point[1] += H;

    std::cout << "1 " << std::flush;
    SDL_UpdateTexture(texture, nullptr, pixels.px, 4*W);
    std::cout << "2 " << std::flush;
    SDL_RenderCopy(renderer, texture, nullptr, nullptr);
    std::cout << "3 " << std::flush;
    SDL_RenderPresent(renderer);
    std::cout << "4\n" << std::flush;

    SDL_Delay(1000 / 60);
    std::cout << "Post delay." << std::endl;
  }

  return 0;
}