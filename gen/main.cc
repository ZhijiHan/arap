#include <iostream>
#include <string>

#include "modelgenerator.h"

// Usage: [model_gen] [.off filename] [.dmat filename] [model name].
int main(int argc, char *argv[]) {
  if (argc < 4) {
    std::cout << "Error: not enough input parameters.";
    return 0;
  }
  arap::gen::ModelGenerator generator;
  std::string model_name(argv[3]);
  if (model_name == "dino") {
    generator.GenerateDino(argv[1], argv[2]);
  } else if (model_name == "spiky-plane") {
    generator.GenerateSpikyPlane(argv[1], argv[2]);
  } else if (model_name == "twist-bar") {
    generator.GenerateTwistBar(argv[1], argv[2]);
  } else if (model_name == "armadillo") {
    generator.GenerateArmadillo(argv[1], argv[2]);
  } else if (model_name == "bunny") {
    generator.GenerateBunny(argv[1], argv[2]);
  }
  return 0;
}
