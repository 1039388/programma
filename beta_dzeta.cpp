#include <filesystem>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>
#include <vector>
#include <tuple>
#include <algorithm>

namespace fs = std::filesystem;

int main()
{
    fs::path rw_dir = fs::current_path();

    std::ofstream out("beta_dzeta.txt");
    out << "# beta dzeta_re dzeta_im\n";

    std::regex beta_regex(R"(_beta=([0-9.eE+-]+))");

    // парсим строку вида:
    // omega:(0.0471089,0.0819412) dzeta:(2.05829e-07,-5.49498e-08) bc:9.30872e-06
    std::regex dzeta_regex(
        R"(dzeta:\(\s*([0-9.eE+-]+)\s*,\s*([0-9.eE+-]+)\s*\))"
    );

    std::vector<std::tuple<double,double,double>> data;

    for (const auto& entry : fs::directory_iterator(rw_dir))
    {
        if (!entry.is_directory())
            continue;

        std::string dirname = entry.path().filename().string();
        std::smatch beta_match;

        if (!std::regex_search(dirname, beta_match, beta_regex))
            continue;

        double beta = std::stod(beta_match[1]);

        fs::path dzeta_file = entry.path() / "omega_dzeta.txt";
        if (!fs::exists(dzeta_file))
            continue;

        std::ifstream in(dzeta_file);
        if (!in)
            continue;

        std::string line;
        std::getline(in, line);

        std::smatch dzeta_match;
        if (!std::regex_search(line, dzeta_match, dzeta_regex))
            continue;

        double dzeta_re = std::stod(dzeta_match[1]);
        double dzeta_im = std::stod(dzeta_match[2]);

        data.emplace_back(beta, dzeta_re, dzeta_im);
    }

    std::sort(data.begin(), data.end(),
              [](auto& a, auto& b){ return std::get<0>(a) < std::get<0>(b); });

    for (auto& [beta, re, im] : data)
        out << beta << " " << re << " " << im << "\n";

    std::cout << "Готово: beta_dzeta.txt в " << rw_dir << "\n";
    return 0;
}
