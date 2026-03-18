#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib> // Для system

int main() {
    std::string folder1 = "r_w_1.1_a(z)_ideal";
    std::string folder2 = "r_w_1.1_a(z)_resist_+_";
    std::string folder3 = "r_w_1.1_a(z)_resist_-_";   // <<< ДОБАВЛЕНО
    std::string fileName = "beta_omega.txt";

    std::string scriptName = "plot_script0.gp";
    std::ofstream script(scriptName);

    if (script.is_open()) {

        script << "set grid\n";
        script << "set title 'ω(β)'\n";
        script << "set xlabel 'beta'\n";
        script << "set ylabel 'ω'\n";
        script << "set key outside\n";

        // --- Расширенная команда plot, содержащая три папки ---
        script << "my_plot_cmd = \"plot '" << folder1 << "/" << fileName
               << "' using 1:($2) with linespoints title 'Ideal Re' lw 2, \\\n"
               << "                      '" << folder1 << "/" << fileName
               << "' using 1:($3) with linespoints title 'Ideal Im' lw 2, \\\n"
               << "                      '" << folder2 << "/" << fileName
               << "' using 1:($2) with linespoints title 'Resist+ Re' lw 2, \\\n"
               << "                      '" << folder2 << "/" << fileName
               << "' using 1:($3) with linespoints title 'Resist+ Im' lw 2, \\\n"
               << "                      '" << folder3 << "/" << fileName
               << "' using 1:($2) with linespoints title 'Resist- Re' lw 2, \\\n"
               << "                      '" << folder3 << "/" << fileName
               << "' using 1:($3) with linespoints title 'Resist- Im' lw 2\" \n";

        // --- PNG вывод ---
        script << "set terminal pngcairo size 800,600 enhanced font 'Verdana,10'\n";
        script << "set output 'comparison_plot_w_1_01_a.png'\n";
        script << "@my_plot_cmd\n";
        script << "unset output\n";

        // --- Qt окно ---
        script << "set terminal qt size 800,600 font 'Verdana,10' title 'Gnuplot Simulation Result'\n";
        script << "@my_plot_cmd\n";

        script << "pause mouse close\n";

        script.close();
        std::cout << "Скрипт gnuplot создан: " << scriptName << std::endl;

        std::cout << "Запуск Gnuplot..." << std::endl;
        int result = system("gnuplot plot_script0.gp");

        if (result == 0) {
            std::cout << "График успешно сохранён и отображён." << std::endl;
        } else {
            std::cerr << "Ошибка при выполнении gnuplot." << std::endl;
        }

    } else {
        std::cerr << "Ошибка создания файла скрипта!" << std::endl;
        return 1;
    }

    return 0;
}
