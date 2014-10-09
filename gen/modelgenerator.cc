#include "modelgenerator.h"

#include <iostream>
#include <fstream>

namespace arap {
namespace gen {

// Generate files for dino model.
void ModelGenerator::GenerateDino(const char* off_file, const char* vertex_tag,
    const char* vertex_file) {
  std::cout << "Dino." << std::endl;
}

// Generate files for spiky plane model.
void ModelGenerator::GenerateSpikyPlane(const char* off_file,
    const char* vertex_tag, const char* vertex_file) {
  std::cout << "Spiky plane." << std::endl;
  // Generate off files.
  std::ofstream off_file_writer;
  off_file_writer.open(off_file);
  off_file_writer << "OFF" << std::endl;
  // [# of vertices] [# of faces] [# of edges].
  int vertex_num_in_row = 21;
  int vertex_num = vertex_num_in_row * vertex_num_in_row;
  int face_num = (vertex_num_in_row - 1) * (vertex_num_in_row - 1) * 2;
  int edge_num = face_num * 3 / 2;
  off_file_writer << vertex_num << " " << face_num << " " << edge_num
    << std::endl;
  // Generate all the vertices.
  double step_size = 1.0;
  for (int i = 0; i < vertex_num_in_row; ++i) {
    double x = i * step_size;
    for (int j = 0; j < vertex_num_in_row; ++j) {
      double y = j * step_size;
      off_file_writer << x << " " << y << " " << "0.0" << std::endl;
    }
  }

  // Generate all the faces.
  for (int i = 0; i < vertex_num_in_row - 1; ++i) {
    for (int j = 0; j < vertex_num_in_row - 1; ++j) {
      // The four vertices around the square (i, j):
      // (i * vertex_num_in_row + j) (i * vertex_num_in_row + j + 1)
      // ((i + 1) * vertex_num_in_row + j) ((i + 1) * vertex_num_in_row + j + 1)
      int upper_left = i * vertex_num_in_row + j;
      int upper_right = upper_left + 1;
      int lower_left = (i + 1) * vertex_num_in_row + j;
      int lower_right = lower_left + 1;
      off_file_writer << 3 << " " << upper_left << " " << lower_right
        << " " << upper_right << std::endl;
      off_file_writer << 3 << " " << upper_left << " " << lower_left
        << " " << lower_right << std::endl;
    }
  }
  off_file_writer.close();
}

// Generate files for twist bar model.
void ModelGenerator::GenerateTwistBar(const char* off_file,
    const char* vertex_tag, const char* vertex_file) {
  std::cout << "Twist bar." << std::endl;
}

// Generate files for armadillo model.
void ModelGenerator::GenerateArmadillo(const char* off_file,
    const char* vertex_tag, const char* vertex_file) {
  std::cout << "Armadillo." << std::endl;
}

// Generate files for Stanford Bunny.
void ModelGenerator::GenerateBunny(const char* off_file,
    const char* vertex_tag, const char* vertex_file) {
  std::cout << "Bunny." << std::endl;
}

}  // namespace gen
}  // namespace arap
