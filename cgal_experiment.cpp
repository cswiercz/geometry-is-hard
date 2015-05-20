
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/IO/File_poly.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> Tr;
typedef CGAL::Delaunay_mesh_size_criteria_2<Tr> Criteria;

typedef K::Point_2 Point;


int main(int argc, char **argv)
{
  Criteria criteria(0.125, 0.5);
  Tr t;

  std::ifstream input(argv[1]);
  if(input) {
    CGAL::read_triangle_poly_file(t, input);
    CGAL::refine_Delaunay_mesh_2(t, criteria);

    std::cout << "Mesh points: " << t.number_of_vertices() << std::endl;
    std::cout << "Mesh triangles: " << t.number_of_faces() << std::endl;
  }

  return 0;
}
