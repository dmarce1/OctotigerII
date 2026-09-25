// Element names and abridged weights: CIAAW 2024.
// Elements without a standard weight use a representative isotope mass number
// (IUPAC table, 4 May 2022). These are defaults, not an isotope network.
#include "octotigerII/composition/species.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <set>
#include <stdexcept>

namespace octotigerII::composition {
namespace {
std::string trim(std::string s) {
	auto begin = s.find_first_not_of(" \t\r\n");
	if (begin == std::string::npos) return {};
	return s.substr(begin, s.find_last_not_of(" \t\r\n") - begin + 1);
}
std::string lower(std::string_view s) {
	std::string out = trim(std::string(s));
	for (auto& c : out) c = char(std::tolower(static_cast<unsigned char>(c)));
	return out;
}
std::vector<std::string> split(std::string const& s, char delimiter) {
	std::vector<std::string> result;
	std::size_t begin = 0;
	for (;;) {
		auto end = s.find(delimiter, begin);
		result.push_back(trim(s.substr(begin, end == std::string::npos ? end : end - begin)));
		if (result.back().empty()) throw std::invalid_argument("Empty species/composition entry");
		if (end == std::string::npos) return result;
		begin = end + 1;
	}
}
Real number(std::string text, bool percentage = false) {
	text = trim(text);
	bool percent = !text.empty() && text.back() == '%';
	if (percent && !percentage) throw std::invalid_argument("Percentage is not an atomic mass or number");
	if (percent) text.pop_back();
	std::size_t n = 0;
	Real const x = std::stod(text, &n) / (percent ? 100 : 1);
	if (n != text.size() || !std::isfinite(x)) throw std::invalid_argument("Invalid composition number: " + text);
	return x;
}
}
std::array<Element, 118> const& periodicTable() {
	static const std::array<Element, 118> table{{
		{"H", "hydrogen", 1.0080, 1, false},
		{"He", "helium", 4.0026, 2, false},
		{"Li", "lithium", 6.94, 3, false},
		{"Be", "beryllium", 9.0122, 4, false},
		{"B", "boron", 10.81, 5, false},
		{"C", "carbon", 12.011, 6, false},
		{"N", "nitrogen", 14.007, 7, false},
		{"O", "oxygen", 15.999, 8, false},
		{"F", "fluorine", 18.998, 9, false},
		{"Ne", "neon", 20.180, 10, false},
		{"Na", "sodium", 22.990, 11, false},
		{"Mg", "magnesium", 24.305, 12, false},
		{"Al", "aluminium", 26.982, 13, false},
		{"Si", "silicon", 28.085, 14, false},
		{"P", "phosphorus", 30.974, 15, false},
		{"S", "sulfur", 32.06, 16, false},
		{"Cl", "chlorine", 35.45, 17, false},
		{"Ar", "argon", 39.95, 18, false},
		{"K", "potassium", 39.098, 19, false},
		{"Ca", "calcium", 40.078, 20, false},
		{"Sc", "scandium", 44.956, 21, false},
		{"Ti", "titanium", 47.867, 22, false},
		{"V", "vanadium", 50.942, 23, false},
		{"Cr", "chromium", 51.996, 24, false},
		{"Mn", "manganese", 54.938, 25, false},
		{"Fe", "iron", 55.845, 26, false},
		{"Co", "cobalt", 58.933, 27, false},
		{"Ni", "nickel", 58.693, 28, false},
		{"Cu", "copper", 63.546, 29, false},
		{"Zn", "zinc", 65.38, 30, false},
		{"Ga", "gallium", 69.723, 31, false},
		{"Ge", "germanium", 72.630, 32, false},
		{"As", "arsenic", 74.922, 33, false},
		{"Se", "selenium", 78.971, 34, false},
		{"Br", "bromine", 79.904, 35, false},
		{"Kr", "krypton", 83.798, 36, false},
		{"Rb", "rubidium", 85.468, 37, false},
		{"Sr", "strontium", 87.62, 38, false},
		{"Y", "yttrium", 88.906, 39, false},
		{"Zr", "zirconium", 91.222, 40, false},
		{"Nb", "niobium", 92.906, 41, false},
		{"Mo", "molybdenum", 95.95, 42, false},
		{"Tc", "technetium", 97, 43, true},
		{"Ru", "ruthenium", 101.07, 44, false},
		{"Rh", "rhodium", 102.91, 45, false},
		{"Pd", "palladium", 106.42, 46, false},
		{"Ag", "silver", 107.87, 47, false},
		{"Cd", "cadmium", 112.41, 48, false},
		{"In", "indium", 114.82, 49, false},
		{"Sn", "tin", 118.71, 50, false},
		{"Sb", "antimony", 121.76, 51, false},
		{"Te", "tellurium", 127.60, 52, false},
		{"I", "iodine", 126.90, 53, false},
		{"Xe", "xenon", 131.29, 54, false},
		{"Cs", "caesium", 132.91, 55, false},
		{"Ba", "barium", 137.33, 56, false},
		{"La", "lanthanum", 138.91, 57, false},
		{"Ce", "cerium", 140.12, 58, false},
		{"Pr", "praseodymium", 140.91, 59, false},
		{"Nd", "neodymium", 144.24, 60, false},
		{"Pm", "promethium", 145, 61, true},
		{"Sm", "samarium", 150.36, 62, false},
		{"Eu", "europium", 151.96, 63, false},
		{"Gd", "gadolinium", 157.25, 64, false},
		{"Tb", "terbium", 158.93, 65, false},
		{"Dy", "dysprosium", 162.50, 66, false},
		{"Ho", "holmium", 164.93, 67, false},
		{"Er", "erbium", 167.26, 68, false},
		{"Tm", "thulium", 168.93, 69, false},
		{"Yb", "ytterbium", 173.05, 70, false},
		{"Lu", "lutetium", 174.97, 71, false},
		{"Hf", "hafnium", 178.49, 72, false},
		{"Ta", "tantalum", 180.95, 73, false},
		{"W", "tungsten", 183.84, 74, false},
		{"Re", "rhenium", 186.21, 75, false},
		{"Os", "osmium", 190.23, 76, false},
		{"Ir", "iridium", 192.22, 77, false},
		{"Pt", "platinum", 195.08, 78, false},
		{"Au", "gold", 196.97, 79, false},
		{"Hg", "mercury", 200.59, 80, false},
		{"Tl", "thallium", 204.38, 81, false},
		{"Pb", "lead", 207.2, 82, false},
		{"Bi", "bismuth", 208.98, 83, false},
		{"Po", "polonium", 209, 84, true},
		{"At", "astatine", 210, 85, true},
		{"Rn", "radon", 222, 86, true},
		{"Fr", "francium", 223, 87, true},
		{"Ra", "radium", 226, 88, true},
		{"Ac", "actinium", 227, 89, true},
		{"Th", "thorium", 232.04, 90, false},
		{"Pa", "protactinium", 231.04, 91, false},
		{"U", "uranium", 238.03, 92, false},
		{"Np", "neptunium", 237, 93, true},
		{"Pu", "plutonium", 244, 94, true},
		{"Am", "americium", 243, 95, true},
		{"Cm", "curium", 247, 96, true},
		{"Bk", "berkelium", 247, 97, true},
		{"Cf", "californium", 251, 98, true},
		{"Es", "einsteinium", 252, 99, true},
		{"Fm", "fermium", 257, 100, true},
		{"Md", "mendelevium", 258, 101, true},
		{"No", "nobelium", 259, 102, true},
		{"Lr", "lawrencium", 266, 103, true},
		{"Rf", "rutherfordium", 267, 104, true},
		{"Db", "dubnium", 268, 105, true},
		{"Sg", "seaborgium", 269, 106, true},
		{"Bh", "bohrium", 270, 107, true},
		{"Hs", "hassium", 269, 108, true},
		{"Mt", "meitnerium", 278, 109, true},
		{"Ds", "darmstadtium", 281, 110, true},
		{"Rg", "roentgenium", 282, 111, true},
		{"Cn", "copernicium", 285, 112, true},
		{"Nh", "nihonium", 286, 113, true},
		{"Fl", "flerovium", 289, 114, true},
		{"Mc", "moscovium", 290, 115, true},
		{"Lv", "livermorium", 293, 116, true},
		{"Ts", "tennessine", 294, 117, true},
		{"Og", "oganesson", 294, 118, true},
	}};
	return table;
}
Element const& element(std::string_view text) {
	auto name = lower(text);
	if (name == "aluminum") name = "aluminium";
	if (name == "cesium") name = "caesium";
	if (name == "sulphur") name = "sulfur";
	for (auto const& e : periodicTable())
		if (name == lower(e.symbol) || name == e.name) return e;
	throw std::invalid_argument("Unknown element: " + std::string(text));
}
void Species::validate() const {
	if (name.empty() || !(std::isalpha(static_cast<unsigned char>(name[0])) || name[0] == '_') ||
		!std::all_of(name.begin(), name.end(), [](unsigned char c) { return std::isalnum(c) || c == '_'; }))
		throw std::invalid_argument("Species name must be an identifier: " + name);
	if (!std::isfinite(atomicMass) || !std::isfinite(atomicNumber) || !std::isfinite(initialFraction) || initialFraction < 0 ||
		(!tracer() && (!(atomicMass > 0) || !(atomicNumber > 0) || atomicNumber > atomicMass)))
		throw std::invalid_argument("Species needs A >= Z > 0, or A=Z=0 for a tracer: " + name);
	if (!mixture.empty()) {
		Real sum = 0, inverseMass = 0, charge = 0;
		std::set<int> seen;
		for (auto const& c : mixture) {
			if (c.atomicNumber < 1 || c.atomicNumber > 118 || !seen.insert(c.atomicNumber).second ||
				!std::isfinite(c.massFraction) || c.massFraction < 0) throw std::invalid_argument("Invalid elemental mixture: " + name);
			auto const& e = periodicTable()[c.atomicNumber - 1];
			sum += c.massFraction;
			inverseMass += c.massFraction / e.atomicMass;
			charge += c.massFraction * e.atomicNumber / e.atomicMass;
		}
		if (std::abs(sum - 1) > 1e-12 || !(inverseMass > 0) ||
			std::abs(atomicMass * inverseMass - 1) > 1e-12 || std::abs(atomicNumber * inverseMass - charge) > 1e-12)
			throw std::invalid_argument("Mixture fractions must sum to one and agree with A,Z: " + name);
	}
}
void Options::validate() const {
	std::set<std::string> names;
	Real sum = 0;
	std::size_t material = 0;
	for (auto const& s : species) {
		s.validate();
		if (!names.insert(lower(s.name)).second) throw std::invalid_argument("Duplicate species name: " + s.name);
		if (!s.tracer()) { sum += s.initialFraction; ++material; }
	}
	if (enabled && (!material || std::abs(sum - 1) > 1e-12))
		throw std::invalid_argument("Material initial mass fractions must sum to one; tracers do not contribute");
}
std::vector<Species> parseSpecies(std::string const& text) {
	std::vector<Species> result;
	for (auto const& definition : split(text, ';')) {
		auto parts = split(definition, ':');
		if (parts.size() != 3) throw std::invalid_argument("Expected name:initialFraction:composition");
		Species s;
		s.name = parts[0];
		s.initialFraction = number(parts[1], true);
		if (parts[2].find('=') == std::string::npos) {
			auto const& e = element(parts[2]);
			s.atomicMass = e.atomicMass; s.atomicNumber = e.atomicNumber;
			s.mixture.push_back({e.atomicNumber, 1});
		} else {
			bool aSeen = false, zSeen = false;
			Real inverseMass = 0, charge = 0, sum = 0;
			for (auto const& entry : split(parts[2], ',')) {
				auto pair = split(entry, '=');
				if (pair.size() != 2) throw std::invalid_argument("Expected element=massFraction or A=mass,Z=number");
				auto key = lower(pair[0]);
				if (key == "a" || key == "atomicmass") {
					if (aSeen) throw std::invalid_argument("Duplicate atomic mass");
					aSeen = true; s.atomicMass = number(pair[1]);
				} else if (key == "z" || key == "atomicnumber") {
					if (zSeen) throw std::invalid_argument("Duplicate atomic number");
					zSeen = true; s.atomicNumber = number(pair[1]);
				} else {
					auto const& e = element(key);
					Real const fraction = number(pair[1], true);
					s.mixture.push_back({e.atomicNumber, fraction});
					sum += fraction;
					inverseMass += fraction / e.atomicMass;
					charge += fraction * e.atomicNumber / e.atomicMass;
				}
			}
			if (aSeen || zSeen) {
				if (!aSeen || !zSeen || !s.mixture.empty()) throw std::invalid_argument("Specify both A and Z, or an elemental mixture");
			} else {
				if (!(inverseMass > 0) || std::abs(sum - 1) > 1e-12) throw std::invalid_argument("Element mass fractions must sum to one (100%)");
				s.atomicMass = 1 / inverseMass;
				s.atomicNumber = charge / inverseMass;
			}
		}
		s.validate();
		result.push_back(std::move(s));
	}
	Options{true, result}.validate();
	return result;
}
std::vector<units::Density> initialDensities(Options const& options, units::Density rho) {
	std::vector<units::Density> result;
	for (auto const& s : options.species) result.push_back(s.initialFraction * rho);
	return result;
}
units::Density totalDensity(Options const& options, std::vector<units::Density> const& values) {
	if (values.size() != options.species.size()) throw std::invalid_argument("Species field count mismatch");
	units::Density rho{};
	for (std::size_t s = 0; s < values.size(); ++s) if (!options.species[s].tracer()) rho += values[s];
	return rho;
}
} // namespace octotigerII::composition
