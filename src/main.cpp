#ifdef _WIN32
#include <SDL.h>
#endif

#ifdef __linux__ 
#include <SDL2/SDL.h>
#endif

#include <cmath>
#include <cassert>

#include <algorithm>
#include <functional>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <vector>

#include "sw_fast_math.hpp"

// x86 is little endian, God help you if you try to use this on anything else
union Color {
  uint32_t u32;
  struct { uint8_t b, g, r, a; };

  Color() : u32(0) { }
  Color(uint32_t u32) : u32(u32) { }
  Color(uint8_t a, uint8_t r, uint8_t g, uint8_t b) : a(a), r(r), g(g), b(b) { }

  operator uint32_t() { return u32; }
  operator int()      { return u32; }
};

struct TwoDimData {
  int w;
};

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

/**
 * NOTE: Camera faces the -z axis like in OpenGL. This makes it so that the
 * right-hand-rule is obeyed (x points ->, y points ^, z points x cross y).
 */
struct Camera {
  fmth::Vec pos; // Position
  fmth::Mat rmx; // Rotation Matrix
  fmth::Vec vel; // Velocity
  fmth::Vec ang; // Angular Velocity

  fmth::Rot rot; // Rotor (probably useless)

  Camera() {
    pos = { 0, 0, 150, 1 };
    vel = { 0, 0, 0, 0 };
    ang = { 0, 0, 0, 0 };
    rmx = fmth::IdentityMat;
  }
};

struct Vertex {
  fmth::Vec pos;
  fmth::Vec uv;
  Color color;
};

struct Triangle {
  std::array<Vertex, 3> v;
  PixelData* tx;
  Color temp_color;
};

struct Controls {
  bool up, down, left, right, fwd, back;
  bool xy, yx, xz, zx, zy, yz;
  bool speed_up;

  void update(std::map<int, bool>& keys) {
    speed_up = keys[SDLK_LSHIFT];

    up    = keys[SDLK_SPACE]  || keys[SDLK_e];
    down  = keys[SDLK_q];
    left  = keys[SDLK_a];
    right = keys[SDLK_d];
    fwd   = keys[SDLK_w];
    back  = keys[SDLK_s];

    yz = keys[SDLK_k];
    zy = keys[SDLK_i];

    zx = keys[SDLK_l];
    xz = keys[SDLK_j];

    xy = keys[SDLK_u];
    yx = keys[SDLK_o];
  }
};

// TODO: Make it work with non-square aspect ratios
// Negative x and y to make right hand coordinates
fmth::Vec perspective_project(fmth::Vec v) {
  return { -v.x/v.z, -v.y/v.z, -v.z, v.w };
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
      
      tx.at(x, y) = Color(0, R, G, B);
    }
  }

  return tx;
}

// Converts from normalized coordinates [-1, 1] to texture coordinates
fmth::Vec unnormalize_coords(PixelData* px, fmth::Vec v) {
  return { px->w / 2.0 * (v.x + 1.0), -px->h / 2.0 * (v.y - 1.0), v.z, v.w };
}

// Converts from texture coordinates to normalized coordinates [-1, 1]
fmth::Vec normalize_coords(PixelData* px, fmth::Vec v) {
  return { v.x * (2.0 / px->w) - 1.0, v.y * (2.0 / -px->h) + 1.0, v.z, v.w };
}

void draw_scanline(
  PixelData* px, double* zbf, int y, int x_start, int x_end, double z, Color c
) {
  auto W = px->w - 1;
  if (y < 0 || y > px->h - 1) return;
  if (x_start > x_end) return;
  if (x_start < 0 && x_end < 0) return;
  if (x_start > W && x_end > W) return; 
  if (x_start < 0) x_start = 0;
  if (x_end > W) x_end = W;

  for (int x = x_start; x <= x_end; x++) {
    if (zbf[y*px->w + x] == 0.0 || z < zbf[y*px->w + x]) {
      px->at(x, y) = c;
      zbf[y*px->w + x] = z;
    }
  }
}

/**
 * Given PixelData, color and 3 Vec<2> with values between -1.0 and 1.0,
 * rasterize_triangle will draw a triangle on the PixelData
 */
void rasterize_triangle(
  PixelData* px, double* zbf, fmth::Vec a, fmth::Vec b, fmth::Vec c, Color color
) {
  // a = unnormalize_coords(px, a);
  // b = unnormalize_coords(px, b);
  // c = unnormalize_coords(px, c);

  if (std::tie(b.y, b.x) < std::tie(a.y, a.x)) { std::swap(a, b); }
  if (std::tie(c.y, c.x) < std::tie(a.y, a.x)) { std::swap(a, c); }
  if (std::tie(c.y, c.x) < std::tie(b.y, b.x)) { std::swap(b, c); }

  if ((int)a.y == (int)c.y) return;

  bool bend_side = (b.y - a.y) * (c.x - a.x) < (b.x - a.x) * (c.y - a.y);
  
  double sl[2];
  double  x[2];
  sl[bend_side]  = (b.x - a.x) / (b.y - a.y);
  sl[!bend_side] = (c.x - a.x) / (c.y - a.y);
  x[bend_side]  = a.x;//+ ((int)a[1]) * sl[bend_side];
  x[!bend_side] = a.x;// + ((int)a[1]) * sl[!bend_side];
  // Totally wrong
  double Z = (a.z + b.z + c.z) / 3.0;

  for (double y = a.y, y_end = b.y; ;y++) {
    if (y >= y_end) {
      if (y >= c.y) break;

      sl[bend_side] = (c.x - b.x) / (c.y - b.y);
      x[bend_side]  = b.x;//y + ((int)y) * sl[bend_side];
      y_end = c.y;
    }

    draw_scanline(px, zbf, y, x[0], x[1], Z, color);

    x[bend_side]  += sl[bend_side];
    x[!bend_side] += sl[!bend_side];
  }
}

std::vector<Triangle> load_from_obj(std::string filename) {
  std::vector<Vertex> verts;
  std::vector<Triangle> tris;

  std::ifstream obj_file(filename);

  if (obj_file.is_open()) {
    std::string s;
    while (getline(obj_file, s)) {
      if (s.length() == 0) continue;
      std::stringstream ss(s);
      std::string id;
      double a, b, c;
      ss >> id >> a >> b >> c;

      if (id == "v") {
        Vertex v;
        v.pos.x = a, v.pos.y = b; v.pos.z = c;
        verts.push_back(v);
      } else if (id == "f") {
        Triangle tri;
        auto u = verts[(int)a - 1].pos - verts[(int)b - 1].pos;
        auto v = verts[(int)c - 1].pos - verts[(int)b - 1].pos;
        auto x = fmth::wedge(u, v).normalized().y;
        tri.v = { verts[(int)a - 1], verts[(int)b - 1],verts[(int)c - 1] };
        // tri.temp_color = Color(0, rand() % 0xff, rand() % 0xff, rand() % 0xff);
        tri.temp_color = Color(0, 0x9f + 0x5f*x, 0x7f + 0x7f*x, 0x7f + 0x7f*x);
        // std::cout << std::hex << (uint32_t)tri.temp_color << std::endl;
        // tri.temp_color = rand() % (0xffffff + 1);
        tris.push_back(tri);
      }
    }
  } else {
    std::cout << "unable to open file";
  }

  return tris;
}

int main (int argc, char* argv[]) {
  const double scale = 1.01;
  const int W = 512 / scale, H = 512 / scale;

  SDL_Window* window = SDL_CreateWindow(
    "Software Renderer", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
    W*scale, H*scale, SDL_WINDOW_RESIZABLE
  );
  SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, 0);
  SDL_Texture* texture = SDL_CreateTexture(
    renderer, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STREAMING, W, H
  );

  // auto bitmap = cool_image();

  // Main loop
  uint64_t curr_time = SDL_GetPerformanceCounter();
  uint64_t last_time = 0;
  uint64_t perf_freq = SDL_GetPerformanceFrequency();
  double dt_ms = 0.0;
  double dt    = 0.0;

  Camera camera;
  Controls controls;

  std::vector<Triangle> tris;

  if (argc != 1) {
    tris = load_from_obj(argv[1]);
  }

  // Triangle tri1; // xy - red
  // tri1.v[0].pos = {1.0, 0.0, 0.0, 1.0};
  // tri1.v[1].pos = {0.0, 1.0, 0.0, 1.0}; 
  // tri1.v[2].pos = {0.0, 0.0, 0.0, 1.0};

  // Triangle tri2; // xz - blue
  // tri2.v[0].pos = {1.0, 0.0, 0.0, 1.0};
  // tri2.v[1].pos = {0.0, 0.0, 1.0, 1.0}; 
  // tri2.v[2].pos = {0.0, 0.0, 0.0, 1.0};

  // Triangle tri3; // yz - green
  // tri3.v[0].pos = {0.0, 1.0, 0.0, 1.0};
  // tri3.v[1].pos = {0.0, 0.0, 1.0, 1.0}; 
  // tri3.v[2].pos = {0.0, 0.0, 0.0, 1.0};

  int count = 20;

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

    controls.update(keys);

    PixelData px(W, H);
    double zbf[W*H] = { 0.0 };

    // Camera transformations
    // 
    fmth::Vec d_ang = {
      (double)(controls.xy - controls.yx), // xy
      (double)(controls.xz - controls.zx), // xz
      (double)(controls.yz - controls.zy), // yz
      0.0
    };
    
    camera.ang = 0.5 * camera.ang + 0.5 * d_ang;

    // Through the magic of Geometric Algebra, rotation is a one-liner!
    if (auto rlen = camera.ang.length(); rlen > 0.0001) {
      auto m = fmth::Rot(camera.ang.normalized(), rlen*0.03).matrix();
      camera.rmx = m * camera.rmx;
    }

    fmth::Vec d_vel = {
      0.1 * (double)(controls.right - controls.left),
      0.1 * (double)(controls.up    - controls.down),
      0.1 * (double)(controls.back   - controls.fwd), // Into the screen is -z
      0.0
    };

    if (controls.speed_up) d_vel *= 10;

    // Back into world space
    camera.vel  = camera.rmx.transpose() * d_vel;
    camera.pos += camera.vel;

    const bool debug = false;
    if constexpr (debug) {
      if (++count >= 20) {
        count = 0;
        std::cout << "POS: " << camera.pos << std::endl;
        std::cout << "VEL: " << camera.vel << std::endl;
        std::cout << "RMX:\n";
        std::cout << camera.rmx << std::endl;
      }
    }

    for (Triangle tri : tris) {
      auto a1 = unnormalize_coords(&px, perspective_project(camera.rmx * (tri.v[0].pos - camera.pos)));
      auto b1 = unnormalize_coords(&px, perspective_project(camera.rmx * (tri.v[1].pos - camera.pos)));
      auto c1 = unnormalize_coords(&px, perspective_project(camera.rmx * (tri.v[2].pos - camera.pos)));
      rasterize_triangle(&px, zbf, a1, b1, c1, tri.temp_color);
    }

    // auto a1 = unnormalize_coords(&px, perspective_project(camera.rmx * (tri1.v[0].pos - camera.pos)));
    // auto b1 = unnormalize_coords(&px, perspective_project(camera.rmx * (tri1.v[1].pos - camera.pos)));
    // auto c1 = unnormalize_coords(&px, perspective_project(camera.rmx * (tri1.v[2].pos - camera.pos)));
    // // if (!(a1.x < 0 || a1.x > W-1 || a1.y < 0 || a1.y > H-1 || b1.x < 0 || b1.x > W-1 || b1.y < 0 || b1.y > H-1 || c1.x < 0 || c1.x > W-1 || c1.y < 0 || c1.y > H-1))
    // rasterize_triangle(&px, zbf, a1, b1, c1, 0xff0000);

    // auto a2 = unnormalize_coords(&px, perspective_project(camera.rmx * (tri2.v[0].pos - camera.pos)));
    // auto b2 = unnormalize_coords(&px, perspective_project(camera.rmx * (tri2.v[1].pos - camera.pos)));
    // auto c2 = unnormalize_coords(&px, perspective_project(camera.rmx * (tri2.v[2].pos - camera.pos)));
    // // if (!(a2.x < 0 || a2.x > W-1 || a2.y < 0 || a2.y > H-1 || b2.x < 0 || b2.x > W-1 || b2.y < 0 || b2.y > H-1 || c2.x < 0 || c2.x > W-1 || c2.y < 0 || c2.y > H-1))
    // rasterize_triangle(&px, zbf, a2, b2, c2, 0x00ff00);

    // auto a3 = unnormalize_coords(&px, perspective_project(camera.rmx * (tri3.v[0].pos - camera.pos)));
    // auto b3 = unnormalize_coords(&px, perspective_project(camera.rmx * (tri3.v[1].pos - camera.pos)));
    // auto c3 = unnormalize_coords(&px, perspective_project(camera.rmx * (tri3.v[2].pos - camera.pos)));
    // // if (!(a3.x < 0 || a3.x > W-1 || a3.y < 0 || a3.y > H-1 || b3.x < 0 || b3.x > W-1 || b3.y < 0 || b3.y > H-1 || c3.x < 0 || c3.x > W-1 || c3.y < 0 || c3.y > H-1))
    // rasterize_triangle(&px, zbf, a3, b3, c3, 0x0000ff);

    SDL_UpdateTexture(texture, nullptr, px.px, 4*W);
    SDL_RenderCopy(renderer, texture, nullptr, nullptr);
    SDL_RenderPresent(renderer);

    SDL_Delay(1000/60);
  }

  return 0;
}