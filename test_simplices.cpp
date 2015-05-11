
#include <iostream>

#include "simplices.hpp"


typedef std::array<double, 2> Point;

int main(int argc, char **argv) {

  Simplices<2, 0> pts;
  Point p = {0.0, 0.0};
  pts.add(p);
  p[0] = 1.0;
  pts.add(p);
  p[1] = 1.0;
  pts.add(p);

  std::cout << "Number of points: " << pts.size() << std::endl;
  for (const auto& p: pts) {
    std::cout << "  " << p[0] << ", " << p[1] << std::endl;
  }
  std::cout << std::endl;


  Simplices<2, 1> lines(pts);
  lines.add(std::array<size_t, 2>{0, 1});
  lines.add(std::array<size_t, 2>{1, 2});
  lines.add(std::array<size_t, 2>{2, 0});

  std::cout << "Number of lines: " << lines.size() << std::endl;
  for (const auto& l: lines) {
    std::cout << "  " << l[0] << ", " << l[1] << std::endl;
  }
  std::cout << std::endl;


  Simplices<2, 2> triangles(lines);
  triangles.add(std::array<size_t, 3>{0, 1, 2});

  std::cout << "Number of triangles: " << triangles.size() << std::endl;
  for (const auto& t: triangles) {
    std::cout << "  " << t[0] << ", " << t[1] << ", " << t[2] << std::endl;
  }


  return 0;
}
