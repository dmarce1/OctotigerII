/** @file
 * @brief Material partial densities and massless passive tracer definitions.
 */
#pragma once
#include <array>
#include <string>
#include <string_view>
#include <vector>
#include "octotigerII/units/cgs.hpp"

namespace octotigerII::composition {
struct Element {
	std::string_view symbol, name;
	Real atomicMass;
	int atomicNumber;
	bool representativeIsotope = false;
};
std::array<Element, 118> const& periodicTable();
Element const& element(std::string_view name);

struct Constituent {
	int atomicNumber = 0;
	Real massFraction = 0;
	template <typename Archive> void serialize(Archive& a, unsigned) { a & atomicNumber & massFraction; }
};
struct Species {
	std::string name;
	Real atomicMass = 0, atomicNumber = 0, initialFraction = 0;
	std::vector<Constituent> mixture;
	bool tracer() const { return atomicMass == 0 && atomicNumber == 0; }
	void validate() const;
	template <typename Archive> void serialize(Archive& a, unsigned) {
		a & name & atomicMass & atomicNumber & initialFraction & mixture;
	}
};
struct Options {
	bool enabled = false;
	std::vector<Species> species;
	void validate() const;
	template <typename Archive> void serialize(Archive& a, unsigned) { a & enabled & species; }
};
/// name:initialFraction:element OR A=mass,Z=number OR He=70%,O=30%.
/// Separate species with semicolons. Percentages are by mass.
std::vector<Species> parseSpecies(std::string const&);
std::vector<units::Density> initialDensities(Options const&, units::Density);
units::Density totalDensity(Options const&, std::vector<units::Density> const&);
} // namespace octotigerII::composition
