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

    std::ofstream out("beta_omega.txt");
    out << "# beta omega_re omega_im\n";

    std::regex beta_regex(R"(_beta=([0-9.eE+-]+))");
    std::regex omega_regex(
        R"(omega:\(\s*([0-9.eE+-]+)\s*,\s*([0-9.eE+-]+)\s*\))"
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

        fs::path omega_file = entry.path() / "omega_dzeta.txt";
        if (!fs::exists(omega_file))
            continue;

        std::ifstream in(omega_file);
        if (!in)
            continue;

        std::string line;
        std::getline(in, line);

        std::smatch omega_match;
        if (!std::regex_search(line, omega_match, omega_regex))
            continue;

        double omega_re = std::stod(omega_match[1]);
        double omega_im = std::stod(omega_match[2]);

        data.emplace_back(beta, omega_re, omega_im);
    }

    std::sort(data.begin(), data.end(),
              [](auto& a, auto& b){ return std::get<0>(a) < std::get<0>(b); });

    for (auto& [beta, re, im] : data)
        out << beta << " " << re << " " << im << "\n";

    std::cout << "Готово: beta_omega.txt в " << rw_dir << "\n";
    return 0;
}
