#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <vector>
#include <algorithm>
#include <cmath>

namespace fs = std::filesystem;

double extract_rw(const std::string& dirname)
{
    std::stringstream ss(dirname);
    std::string token;

    std::getline(ss, token, '_'); // r
    std::getline(ss, token, '_'); // w
    std::getline(ss, token, '_'); // значение r_w

    return std::stod(token);
}

int main()
{
    double target_beta;
    std::cout << "Enter beta: ";
    std::cin >> target_beta;

    std::vector<std::pair<double,double>> data; 
    // pair<r_w, omega_im>

    for (const auto& entry : fs::directory_iterator("."))
    {
        if (!entry.is_directory()) continue;

        std::string dirname = entry.path().filename().string();

        if (dirname.rfind("r_w_", 0) != 0) continue;

        double rw = extract_rw(dirname);

        std::string filepath = entry.path().string() + "/beta_omega.txt";
        std::ifstream file(filepath);

        if (!file) continue;

        std::string line;
        std::getline(file, line); // пропустить заголовок

        double beta, omega_re, omega_im;

        while (file >> beta >> omega_re >> omega_im)
        {
            if (std::abs(beta - target_beta) < 1e-6)
            {
                data.push_back({rw, omega_im});
                break;
            }
        }
    }

    // сортировка по r_w
    std::sort(data.begin(), data.end(),
        [](const auto& a, const auto& b)
        {
            return a.first < b.first;
        });

    std::ostringstream filename;
    filename << "omega_im_vs_rw_beta_" << target_beta << ".txt";
    std::ofstream out(filename.str());

    for (auto& p : data)
    {
        out << p.second << " " << p.first << std::endl; 
        // omega_im r_w
    }

    out.close();
}