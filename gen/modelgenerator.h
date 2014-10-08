#ifndef _ARAP_GEN_MODELGENERATOR_H_
#define _ARAP_GEN_MODELGENERATOR_H_

namespace arap {
namespace gen {

class ModelGenerator {
 public:
  // Generate files for dino model.
  void GenerateDino(const char* off_file_name, const char* dmat_file_name);

  // Generate files for spiky plane model.
  void GenerateSpikyPlane(
      const char* off_file_name, const char* dmat_file_name);

  // Generate files for twist bar model.
  void GenerateTwistBar(const char* off_file_name, const char* dmat_file_name);

  // Generate files for armadillo model.
  void GenerateArmadillo(
      const char* off_file_name, const char* dmat_file_name);

  // Generate files for Stanford Bunny.
  void GenerateBunny(const char* off_file_name, const char* dmat_file_name);

 private:
};

}  // namespace gen
}  // namespace arap

#endif  // _ARAP_GEN_MODELGENERATOR_H_
