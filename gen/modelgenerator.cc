#include "modelgenerator.h"

#include <iostream>

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
