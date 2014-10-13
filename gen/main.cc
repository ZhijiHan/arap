#include <iostream>
#include <string>

#include "modelgenerator.h"

// Usage: [model_gen] [model name] [.off filename] [.dmat filename]
// [.dmat filename]
int main(int argc, char *argv[]) {
  if (argc < 5) {
    std::cout << "Error: not enough input parameters.";
    return 0;
  }
  arap::gen::ModelGenerator generator;
  std::string model_name(argv[1]);
  if (model_name == "dino") {
    generator.GenerateDino(argv[2], argv[3], argv[4]);
  } else if (model_name == "spiky-plane") {
    generator.GenerateSpikyPlane(argv[2], argv[3], argv[4]);
  } else if (model_name == "twist-bar") {
    generator.GenerateTwistBar(argv[2], argv[3], argv[4]);
  } else if (model_name == "armadillo") {
    generator.GenerateArmadillo(argv[2], argv[3], argv[4]);
  } else if (model_name == "bunny") {
    generator.GenerateBunny(argv[2], argv[3], argv[4]);
  }
  return 0;
}
