#ifdef _WIN32
#include <SDL.h>
#endif

#ifdef __linux__ 
#include <SDL2/SDL.h>
#endif

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <iomanip>
#include <map>

#include "sw_math.hpp"

struct PixelData {
  int w, h;
  uint32_t* px;

  PixelData(std::size_t w, std::size_t h) : w(w), h(h) {
    px = new uint32_t[w * h]();
  }

  ~PixelData() {
    delete[] px;
  }

  inline uint32_t& at(int x, int y) {
    return px[y*w + x];
  }
};

struct Camera {
  math::Vec<4> pos;
  math::Matrix<4, 4> tform;
  // math::Rotor3 rot;

  math::Vec<4> pos_vel;
  math::BiVec3 rot_vel;


  Camera() {
    pos = { 0.0, 0.0, -5.0, 1.0 };
    // rot = math::Rotor3();
    pos_vel = { 0.0, 0.0, 0.0, 0.0 };
    rot_vel = { 0.0, 0.0, 0.0 };
    tform = rot.to_matrix4();
  }
};

struct Vertex {
  math::Vec<4> pos;
  math::Vec<2> uv;
  uint32_t color;
};

struct Triangle {
  std::array<Vertex, 3> v;
  PixelData* tx;
};

struct Controls {
  bool up, down, left, right, fwd, back;
  bool rup, rdown, rleft, rright, rfwd, rback;

  void update(std::map<int, bool>& keys) {
    up    = keys[SDLK_SPACE]  || keys[SDLK_e];
    down  = keys[SDLK_LSHIFT] || keys[SDLK_q];
    left  = keys[SDLK_a];
    right = keys[SDLK_d];
    fwd   = keys[SDLK_w];
    back  = keys[SDLK_s];

    rup    = keys[SDLK_i];
    rdown  = keys[SDLK_k];
    rleft  = keys[SDLK_j];
    rright = keys[SDLK_l];
    rfwd   = keys[SDLK_o];
    rback  = keys[SDLK_u];
  }
};

math::Vec<2> perspective_project(math::Vec<4> v) {
  // TODO: Make it work with non-square aspect ratios
  return { v[0] / v[2], v[1] / v[2] };
}

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
math::Vec<2> unnormalize_coords(PixelData& px, math::Vec<2> v) {
  return { px.w / 2.0 * (v[0] + 1.0), -px.h / 2.0 * (v[1] - 1.0) };
}

// Converts from texture coordinates to normalized coordinates [-1, 1]
math::Vec<2> normalize_coords(PixelData& px, math::Vec<2> v) {
  return { v[0] * (2.0 / px.w) - 1.0, v[1] * (2.0 / -px.h) + 1.0 };
}

void draw_scanline(PixelData& pixels, int y, int x_start, int x_end, uint32_t c) {
  // std::cout << pixels.px << std::endl;
  y       = std::clamp(y,       0, pixels.h - 1);
  x_start = std::clamp(x_start, 0, pixels.w - 1);
  x_end   = std::clamp(x_end  , 0, pixels.w - 1);
  for (int x = x_start; x <= x_end; x++) {
    pixels.at(x, y) = c;
  }
}

/**
 * Given PixelData, color and 3 Vec<2> with values between -1.0 and 1.0,
 * rasterize_triangle will draw a triangle on the PixelData
 */
void rasterize_triangle(
  PixelData& px, math::Vec<2> a, math::Vec<2> b, math::Vec<2> c, uint32_t color
) {
  // a = unnormalize_coords(px, a);
  // b = unnormalize_coords(px, b);
  // c = unnormalize_coords(px, c);

  if (std::tie(b[1], b[0]) < std::tie(a[1], a[0])) { std::swap(a, b); }
  if (std::tie(c[1], c[0]) < std::tie(a[1], a[0])) { std::swap(a, c); }
  if (std::tie(c[1], c[0]) < std::tie(b[1], b[0])) { std::swap(b, c); }

  if ((int)a[1] == (int)c[1]) return;

  bool bend_side = 
    (b[1] - a[1]) * (c[0] - a[0]) < (b[0] - a[0]) * (c[1] - a[1]);
  
  double sl[2];
  double  x[2];
  sl[bend_side]  = (b[0] - a[0]) / (b[1] - a[1]);
  sl[!bend_side] = (c[0] - a[0]) / (c[1] - a[1]);
  x[bend_side] = a[0];//+ ((int)a[1]) * sl[bend_side];
  x[!bend_side] = a[0];// + ((int)a[1]) * sl[!bend_side];


  for (double y = a[1], y_end = b[1]; ;y++) {
    if (y >= y_end) {
      if (y >= c[1]) break;

      sl[bend_side] = (c[0] - b[0]) / (c[1] - b[1]);
      x[bend_side]  = b[0];//y + ((int)y) * sl[bend_side];
      y_end = c[1];
    }

    draw_scanline(px, y, x[0], x[1], color);

    x[bend_side]  += sl[bend_side];
    x[!bend_side] += sl[!bend_side];
  }
}

int main (int argc, char* argv[]) {
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

  Camera camera;
  Controls controls;

  Triangle tri1;
  tri1.v[0].pos = {0.0, 0.0, 0.0, 1.0};
  tri1.v[1].pos = {1.0, 0.0, 0.0, 1.0};
  tri1.v[2].pos = {0.0, 1.0, 0.0, 1.0};

  Triangle tri2;
  tri2.v[0].pos = {0.0, 0.0, 0.0, 1.0};
  tri2.v[1].pos = {1.0, 0.0, 0.0, 1.0};
  tri2.v[2].pos = {0.0, 0.0, 1.0, 1.0};

  Triangle tri3;
  tri3.v[0].pos = {0.0, 0.0, 0.0, 1.0};
  tri3.v[1].pos = {0.0, 1.0, 0.0, 1.0};
  tri3.v[2].pos = {0.0, 0.0, 1.0, 1.0};

  int count = 20;

  double aa=0.094, ab=-0.11, ac=-0.021, ad=1.0;

  std::map<int, bool> keys;
  while(!keys[SDLK_ESCAPE]) {
    last_time = curr_time;
    curr_time = SDL_GetPerformanceCounter();
    perf_freq = SDL_GetPerformanceFrequency();
    dt_ms = (double)((curr_time - last_time) * 1000 / (double)perf_freq);
    dt    = dt_ms / 1000.0;

    SDL_Event ev;
    while(SDL_PollEvent(&ev)) {
      switch(ev.type) {
        case SDL_QUIT:    keys[SDLK_ESCAPE] = true; break;
        case SDL_KEYDOWN: keys[ev.key.keysym.sym] = true; break;
        case SDL_KEYUP:   keys[ev.key.keysym.sym] = false; break;
      }
    }

    PixelData px(W, H);

    controls.update(keys);

    math::Vec<3> delta_rot_vel = {
      (double)(controls.rfwd   - controls.rback),  // xy
      (double)(controls.rright - controls.rleft), // xz
      (double)(controls.rup    - controls.rdown), // yz
    };
    
    camera.rot_vel = 0.5 * camera.rot_vel + 0.5 * delta_rot_vel;

    if (auto rlen = camera.rot_vel.length(); rlen > 0.0001) {
      camera.tform = math::Rotor3(camera.rot_vel.normal(), rlen*0.03).to_matrix4() * camera.tform;
    }

    math::Vec<4> d_vel = {
      0.1 * (double)(controls.right - controls.left),
      0.1 * (double)(controls.up    - controls.down),
      0.1 * (double)(controls.fwd   - controls.back),
      0.0
    };

    // Back into world space
    camera.pos_vel = camera.tform.transpose() * d_vel;
    camera.pos += camera.pos_vel;

    auto a1 = unnormalize_coords(px, perspective_project(camera.tform * (tri1.v[0].pos - camera.pos)));
    auto b1 = unnormalize_coords(px, perspective_project(camera.tform * (tri1.v[1].pos - camera.pos)));
    auto c1 = unnormalize_coords(px, perspective_project(camera.tform * (tri1.v[2].pos - camera.pos)));
    if (!(a1[0] < 0 || a1[0] > W-1 || a1[1] < 0 || a1[1] > H-1 || b1[0] < 0 || b1[0] > W-1 || b1[1] < 0 || b1[1] > H-1 || c1[0] < 0 || c1[0] > W-1 || c1[1] < 0 || c1[1] > H-1))
      rasterize_triangle(px, a1, b1, c1, 0xff0000);

    auto a2 = unnormalize_coords(px, perspective_project(camera.tform * (tri2.v[0].pos - camera.pos)));
    auto b2 = unnormalize_coords(px, perspective_project(camera.tform * (tri2.v[1].pos - camera.pos)));
    auto c2 = unnormalize_coords(px, perspective_project(camera.tform * (tri2.v[2].pos - camera.pos)));
    if (!(a2[0] < 0 || a2[0] > W-1 || a2[1] < 0 || a2[1] > H-1 || b2[0] < 0 || b2[0] > W-1 || b2[1] < 0 || b2[1] > H-1 || c2[0] < 0 || c2[0] > W-1 || c2[1] < 0 || c2[1] > H-1))
      rasterize_triangle(px, a2, b2, c2, 0x0000ff);

    auto a3 = unnormalize_coords(px, perspective_project(camera.tform * (tri3.v[0].pos - camera.pos)));
    auto b3 = unnormalize_coords(px, perspective_project(camera.tform * (tri3.v[1].pos - camera.pos)));
    auto c3 = unnormalize_coords(px, perspective_project(camera.tform * (tri3.v[2].pos - camera.pos)));
    if (!(a3[0] < 0 || a3[0] > W-1 || a3[1] < 0 || a3[1] > H-1 || b3[0] < 0 || b3[0] > W-1 || b3[1] < 0 || b3[1] > H-1 || c3[0] < 0 || c3[0] > W-1 || c3[1] < 0 || c3[1] > H-1))
      rasterize_triangle(px, a3, b3, c3, 0x0dff00);

    SDL_UpdateTexture(texture, nullptr, px.px, 4*W);
    SDL_RenderCopy(renderer, texture, nullptr, nullptr);
    SDL_RenderPresent(renderer);

    SDL_Delay(1000/60);
  }

  return 0;
}