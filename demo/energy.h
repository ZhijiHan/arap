#ifndef _ARAP_DEMO_ENERGY_H_
#define _ARAP_DEMO_ENERGY_H_

#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace arap {
namespace demo {

typedef std::map<const std::string, double> EnergyMap;

class Energy {
 public:
  Energy() {}

  // Overload output function.
  friend std::ostream& operator<<(std::ostream& output, const Energy& energy);

  // Overload input function.
  friend std::istream& operator>>(std::istream& input, Energy& energy);

  // Insert energy type and value. Returns whether the new type has been added
  // successfully. If there is an energy type with the same value in the map,
  // return false.
  bool AddEnergyType(const std::string type, const double value);

  // Returns the value corresponds to the given type.
  double GetEnergyValue(const std::string type) const;

  // Get all the energy names.
  std::vector<std::string> GetEnergyTypes() const;

 private:
  // Maps energy type name and value.
  EnergyMap energy_map_;
};

}  // namespace demo
}  // namespace arap

#endif  // _ARAP_DEMO_ENERGY_H_
