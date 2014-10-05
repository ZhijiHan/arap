#include "energy.h"

namespace arap {
namespace demo {

std::ostream& operator<<(std::ostream& output, const Energy& energy) {
  for (EnergyMap::const_iterator it = energy.energy_map_.begin();
      it != energy.energy_map_.end(); ++it) {
    output << it->first << ": " << it->second << "\t";
  }
  output << std::endl;
  return output;
}

std::istream& operator>>(std::istream& input, Energy& energy) {
  std::string type;
  double value;
  input >> type >> value;
  energy.AddEnergyType(type, value);
  return input;
}

bool Energy::AddEnergyType(const std::string type, const double value) {
  return energy_map_.insert(
      std::pair<const std::string, double>(type, value)).second;
}

double Energy::GetEnergyValue(const std::string type) const {
  EnergyMap::const_iterator it = energy_map_.find(type);
  if (it != energy_map_.end()) {
    return it->second;
  }
  std::cout << "Warning: this energy type does not exist!" << std::endl;
  return 0.0;
}

}  // namespace demo
}  // namespace arap
