#include "modelgenerator.h"

#include <iostream>

namespace arap {
namespace gen {

// Generate files for dino model.
void ModelGenerator::GenerateDino(
    const char* off_file_name, const char* dmat_file_name) {
  std::cout << "Generate dino." << std::endl;
}

// Generate files for spiky plane model.
void ModelGenerator::GenerateSpikyPlane(
    const char* off_file_name, const char* dmat_file_name) {
  std::cout << "Generate spiky plane." << std::endl;
}

// Generate files for twist bar model.
void ModelGenerator::GenerateTwistBar(
    const char* off_file_name, const char* dmat_file_name) {
  std::cout << "Generate twist bar." << std::endl;
}

// Generate files for armadillo model.
void ModelGenerator::GenerateArmadillo(
    const char* off_file_name, const char* dmat_file_name) {
  std::cout << "Generate armadillo." << std::endl;
}

// Generate files for Stanford Bunny.
void ModelGenerator::GenerateBunny(
    const char* off_file_name, const char* dmat_file_name) {
  std::cout << "Generate bunny." << std::endl;
}

}  // namespace gen
}  // namespace arap
