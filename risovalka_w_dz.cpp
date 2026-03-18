#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

int main() {
    std::string folder1 = "r_w_1.1_a(z)_ideal";
    std::string folder2 = "r_w_1.1_a(z)_resist_+_";
    std::string folder3 = "r_w_1.1_a(z)_resist_-_";   // <<< ДОБАВЛЕНО

    std::string fileOmega = "beta_omega.txt";
    std::string fileDzeta = "beta_dzeta.txt";

    std::string scriptName = "plot_script0.gp";
    std::ofstream script(scriptName);

    if (script.is_open()) {

        script << "set grid\n";
        script << "set key outside\n";
        script << "set title 'w(β) и dzeta(β) '\n";
        script << "set xlabel 'beta'\n";
        script << "set ylabel 'values'\n";

        // ============================================================
        // Макрос для объединённого графика
        // ============================================================
        script <<
        "plot_cmd = \"plot \\\n"
        // ---- omega ----
        "'"<< folder1 << "/" << fileOmega << "' using 1:2 with linespoints title 'Ideal Re w' lw 2, \\\n"
        "'"<< folder1 << "/" << fileOmega << "' using 1:3 with linespoints title 'Ideal Im w' lw 2, \\\n"

        "'"<< folder2 << "/" << fileOmega << "' using 1:2 with linespoints title 'Resist+ Re w' lw 2, \\\n"
        "'"<< folder2 << "/" << fileOmega << "' using 1:3 with linespoints title 'Resist+ Im w' lw 2, \\\n"

        "'"<< folder3 << "/" << fileOmega << "' using 1:2 with linespoints title 'Resist- Re w' lw 2, \\\n"
        "'"<< folder3 << "/" << fileOmega << "' using 1:3 with linespoints title 'Resist- Im w' lw 2, \\\n"

        // ---- dzeta ----
        "'"<< folder1 << "/" << fileDzeta << "' using 1:2 with linespoints title 'Ideal Re dzeta' lw 2, \\\n"
        "'"<< folder1 << "/" << fileDzeta << "' using 1:3 with linespoints title 'Ideal Im dzeta' lw 2, \\\n"

        "'"<< folder2 << "/" << fileDzeta << "' using 1:2 with linespoints title 'Resist+ Re dzeta' lw 2, \\\n"
        "'"<< folder2 << "/" << fileDzeta << "' using 1:3 with linespoints title 'Resist+ Im dzeta' lw 2, \\\n"

        "'"<< folder3 << "/" << fileDzeta << "' using 1:2 with linespoints title 'Resist- Re dzeta' lw 2, \\\n"
        "'"<< folder3 << "/" << fileDzeta << "' using 1:3 with linespoints title 'Resist- Im dzeta' lw 2\"\n";

        // ============================================================
        // PNG-файл
        // ============================================================
        script << "set terminal pngcairo size 900,700 font 'Verdana,10'\n";
        script << "set output 'combined_plot_omega_zeta.png'\n";
        script << "@plot_cmd\n";
        script << "unset output\n";

        // ============================================================
        // Окно gnuplot
        // ============================================================
        script << "set terminal qt size 900,700 font 'Verdana,10'\n";
        script << "@plot_cmd\n";
        script << "pause mouse close\n";

        script.close();
        std::cout << "Скрипт gnuplot создан: " << scriptName << std::endl;

        int result = system("gnuplot plot_script0.gp");
        if (result == 0)
            std::cout << "График успешно сохранён и отображён." << std::endl;
        else
            std::cerr << "Ошибка при выполнении gnuplot." << std::endl;

    } else {
        std::cerr << "Ошибка создания файла скрипта!" << std::endl;
        return 1;
    }

    return 0;
}
