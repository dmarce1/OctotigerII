// Copyright (c) 2026 AUTHORS
// Distributed under the Boost Software License, Version 1.0.

#include "octotigerII/finiteVolume/solver.hpp"


namespace octotigerII::finiteVolume {

std::string_view finiteVolumeSchemeName() {
	return "unsplit MUSCL-Hancock with face-centered fluxes";
}

}	 // namespace octotigerII::finiteVolume
