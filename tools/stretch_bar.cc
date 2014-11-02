#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include <vector>

struct Vertex {
  Vertex() { x = y = z = 0.0; }
  Vertex(double xx, double yy, double zz)
    : x(xx), y(yy), z(zz) {}

  static Vertex GetMinVertex(const Vertex& v1, const Vertex& v2) {
    Vertex v;
    v.x = std::min(v1.x, v2.x);
    v.y = std::min(v1.y, v2.y);
    v.z = std::min(v1.z, v2.z);
    return v;
  }

  static Vertex GetMaxVertex(const Vertex& v1, const Vertex& v2) {
    Vertex v;
    v.x = std::max(v1.x, v2.x);
    v.y = std::max(v1.y, v2.y);
    v.z = std::max(v1.z, v2.z);
    return v;
  }

  double x, y, z;
};

// This program reads bar1.off files, then stretch the top and bottom squares.
// Usage: stretch_bar [input_off file] [delta] [output_dmat file].
// => stretch the top and bottom squares by delta.
int main(int argc, char* argv[]) {
  std::ifstream off_file(argv[1]);
  double delta = atof(argv[2]);
  std::ofstream dmat_file(argv[3]);

  // Read OFF header.
  std::string header;
  off_file >> header;
  std::cout << "File header: " << header << std::endl;
  // Read vertex number, face number and edge number.
  int vertex_num, face_num, edge_num;
  off_file >> vertex_num >> face_num >> edge_num;
  std::cout << "Vertex number: " << vertex_num << std::endl
            << "Face number: " << face_num << std::endl
            << "Edge number: " << edge_num << std::endl;

  // Read all the vertices.
  std::vector<Vertex> vertices;
  for (int i = 0; i < vertex_num; ++i) {
    Vertex v;
    off_file >> v.x >> v.y >> v.z;
    vertices.push_back(v);
  }
  // Stretch the top and bottom in y axis.
  const int square_num = 7;
  for (int i = 0; i < square_num * square_num; ++i) {
    Vertex v = vertices[i];
    v.y += delta;
    vertices[i] = v;
  }
  int offset = square_num * square_num;
  for (int i = offset; i < offset + square_num * square_num; ++i) {
    Vertex v = vertices[i];
    v.y -= delta;
    vertices[i] = v;
  }

  // Write down the top and bottom vertices.
  dmat_file << 3 << " " << 2 * square_num * square_num << std::endl;
  for (int i = 0; i < 2 * square_num * square_num; ++i) {
    Vertex v = vertices[i];
    dmat_file << v.x << std::endl;
  }
  for (int i = 0; i < 2 * square_num * square_num; ++i) {
    Vertex v = vertices[i];
    dmat_file << v.y << std::endl;
  }
  for (int i = 0; i < 2 * square_num * square_num; ++i) {
    Vertex v = vertices[i];
    dmat_file << v.z << std::endl;
  }
  dmat_file.close();
  off_file.close();
  return 0;
}
