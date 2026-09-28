#include "direct.hpp"
#include <iomanip>
#include <iostream>
#include <octotigerII/helmholtz/helmholtz.hpp>
#include <string>
int main(int argc, char **argv) {
    try {
        std::cout << std::scientific << std::setprecision(17);
        double d, t, a, z;
        if (argc == 2 && std::string(argv[1]) == "--direct") {
            while (std::cin >> d >> t >> a >> z) {
                auto q = octotigerII::helmholtz::detail::direct(d, t);
                std::cout << q.p << ' ' << q.e << ' ' << q.s << ' ' << q.pd << ' ' << q.pt << ' ' << q.pdd
                          << ' ' << q.pdt << ' ' << q.ptt << ' ' << q.st << ' ' << q.sdd << ' ' << q.eta
                          << ' ' << q.etad << ' ' << q.etat << ' ' << q.etadt << ' ' << q.n << ' ' << q.nd
                          << ' ' << q.nt << ' ' << q.ndt << ' ' << '\n';
            }
        } else {
            if (argc != 2 && argc != 8)
                return 2;
            octotigerII::helmholtz::Grid grid;
            if (argc == 8)
                grid = {std::stoi(argv[2]), std::stoi(argv[3]), std::stod(argv[4]),
                        std::stod(argv[5]), std::stod(argv[6]), std::stod(argv[7])};
            octotigerII::helmholtz::Eos eos(argv[1], grid);
            while (std::cin >> d >> t >> a >> z) {
                const auto q = eos.evaluate(d, t, a, z);
                std::cout << q.ptot << ' ' << q.dpt << ' ' << q.dpd << ' ' << q.dpa << ' ' << q.dpz << ' '
                          << q.etot << ' ' << q.det << ' ' << q.ded << ' ' << q.dea << ' ' << q.dez << ' '
                          << q.stot << ' ' << q.dst << ' ' << q.dsd << ' ' << q.dsa << ' ' << q.dsz << ' '
                          << q.pgas << ' ' << q.dpgast << ' ' << q.dpgasd << ' ' << q.dpgasa << ' '
                          << q.dpgasz << ' ' << q.egas << ' ' << q.degast << ' ' << q.degasd << ' '
                          << q.degasa << ' ' << q.degasz << ' ' << q.sgas << ' ' << q.dsgast << ' '
                          << q.dsgasd << ' ' << q.dsgasa << ' ' << q.dsgasz << ' ' << q.prad << ' '
                          << q.dpradt << ' ' << q.dpradd << ' ' << q.dprada << ' ' << q.dpradz << ' '
                          << q.erad << ' ' << q.deradt << ' ' << q.deradd << ' ' << q.derada << ' '
                          << q.deradz << ' ' << q.srad << ' ' << q.dsradt << ' ' << q.dsradd << ' '
                          << q.dsrada << ' ' << q.dsradz << ' ' << q.pion << ' ' << q.dpiont << ' '
                          << q.dpiond << ' ' << q.dpiona << ' ' << q.dpionz << ' ' << q.eion << ' '
                          << q.deiont << ' ' << q.deiond << ' ' << q.deiona << ' ' << q.deionz << ' '
                          << q.sion << ' ' << q.dsiont << ' ' << q.dsiond << ' ' << q.dsiona << ' '
                          << q.dsionz << ' ' << q.xni << ' ' << q.pele << ' ' << q.ppos << ' ' << q.dpept
                          << ' ' << q.dpepd << ' ' << q.dpepa << ' ' << q.dpepz << ' ' << q.eele << ' '
                          << q.epos << ' ' << q.deept << ' ' << q.deepd << ' ' << q.deepa << ' ' << q.deepz
                          << ' ' << q.sele << ' ' << q.spos << ' ' << q.dsept << ' ' << q.dsepd << ' '
                          << q.dsepa << ' ' << q.dsepz << ' ' << q.xnem << ' ' << q.xne << ' ' << q.dxnet
                          << ' ' << q.dxned << ' ' << q.dxnea << ' ' << q.dxnez << ' ' << q.xnp << ' '
                          << q.zeff << ' ' << q.etaele << ' ' << q.detat << ' ' << q.detad << ' ' << q.detaa
                          << ' ' << q.detaz << ' ' << q.etapos << ' ' << q.pcou << ' ' << q.dpcout << ' '
                          << q.dpcoud << ' ' << q.dpcoua << ' ' << q.dpcouz << ' ' << q.ecou << ' '
                          << q.decout << ' ' << q.decoud << ' ' << q.decoua << ' ' << q.decouz << ' '
                          << q.scou << ' ' << q.dscout << ' ' << q.dscoud << ' ' << q.dscoua << ' '
                          << q.dscouz << ' ' << q.plasg << ' ' << q.dse << ' ' << q.dpe << ' ' << q.dsp << ' '
                          << q.cv_gas << ' ' << q.cp_gas << ' ' << q.gam1_gas << ' ' << q.gam2_gas << ' '
                          << q.gam3_gas << ' ' << q.nabad_gas << ' ' << q.cs_gas << ' ' << q.cv << ' ' << q.cp
                          << ' ' << q.gam1 << ' ' << q.gam2 << ' ' << q.gam3 << ' ' << q.nabad << ' ' << q.cs
                          << ' ' << '\n';
            }
        }
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
