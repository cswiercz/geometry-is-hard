
#include <iostream>

#include <fstream>
#include <sstream>
#include <string>

#include "simplices.hpp"

using namespace std;

Simplices<2, 2> read_triangle_mesh(const std::string& filename) {

  // Read in all the nodes of the triangulation
  Simplices<2, 0> pts;
  size_t nn, _;

  ifstream node_file(filename + ".node");
  node_file >> nn >> _ >> _ >> _;

  std::array<double, 2> x;
  for (size_t i = 0; i < nn; ++i) {
    node_file >> _ >> x[0] >> x[1] >> _;
    pts.add(x);
  }

  node_file.close();


  // Read in all the edges of the triangulation
  Simplices<2, 1> edges(pts);
  size_t ne;

  ifstream edge_file(filename + ".edge");
  edge_file >> ne >> _;

  std::array<size_t, 2> edge;
  for (size_t i = 0; i < ne; ++i) {
    edge_file >> _ >> edge[0] >> edge[1] >> _;
    edge[0] -= 1;
    edge[1] -= 1;
    edges.add(edge);
  }

  edge_file.close();


  // Read in all the elements of the triangulation
  Simplices<2, 2> triangles(edges);
  size_t nt;

  ifstream ele_file(filename + ".ele");
  ele_file >> nt >> _ >> _;

  std::array<size_t, 3> ele;
  for (size_t i = 0; i < nt; ++i) {
    ele_file >> _ >> ele[0] >> ele[1] >> ele[2];
    ele[0] -= 1;  ele[1] -= 1;  ele[2] -= 1;

    // Now what?
  }

  return triangles;
}


