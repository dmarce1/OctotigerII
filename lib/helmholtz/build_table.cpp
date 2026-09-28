#include <iostream>
#include <octotigerII/helmholtz/helmholtz.hpp>
#include <stdexcept>
int main(int argc, char **argv) {
    try {
        if (argc != 2 && argc != 8 && argc != 9) {
            std::cerr << "Usage: helmholtz_build_table OUTPUT [ND NT LOG_D_MIN LOG_D_MAX LOG_T_MIN LOG_T_MAX "
                         "[DERIVATIVE_STEP]]\n";
            return 2;
        }
        octotigerII::helmholtz::Grid g;
        if (argc >= 8)
            g = {std::stoi(argv[2]), std::stoi(argv[3]), std::stod(argv[4]),
                 std::stod(argv[5]), std::stod(argv[6]), std::stod(argv[7])};
        const auto report =
            octotigerII::helmholtz::generateTable(argv[1], g, argc == 9 ? std::stod(argv[8]) : 1e-3);
        std::cout << "Generated " << g.densityPoints << " x " << g.temperaturePoints
                  << " Timmes-format table at " << argv[1] << '\n';
        std::cout << "Sign diagnostics: " << report.nonpositiveElectronHeatCapacityNodes
                  << " nonpositive electron heat capacities; " << report.negativeElectronCompressibilityNodes
                  << " negative electron compressibilities. Metadata: " << argv[1] << ".meta.json\n";
        if (report.nonpositiveElectronHeatCapacityNodes || report.negativeElectronCompressibilityNodes)
            std::cerr << "The original direct EOS loses accuracy in parts of this grid; "
                         "this generated table is not validated for production use over its full range.\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
