#include <iostream>
#include <fstream>
#include <string>
#include <vector>

int main() {
    // Названия папок из вашего примера
    std::string folder1 = "r_w_1.01_a(z)_ideal";
    std::string folder2 = "r_w_1.01_a(z)_resist_+_";
    std::string fileName = "beta_omega.txt";

    // Имя файла скрипта gnuplot
    std::string scriptName = "plot_script.gp";

    std::ofstream script(scriptName);

    if (script.is_open()) {
        script << "set terminal pngcairo size 800,600 enhanced font 'Verdana,10'\n";
        script << "set output 'modul_w2.png'\n";
        script << "set title '|ω|²(β)'\n";
        script << "set xlabel 'beta'\n";
        script << "set ylabel '|ω|²'\n";
        script << "set grid\n";
        script << "set key outside\n";

        // Формула: берем 2-ю колонку (re) и 3-ю (im), возводим в квадрат и суммируем
        // ($2**2 + $3**2)
        script << "plot '" << folder1 << "/" << fileName << "' using 1:($2**2+$3**2) "
               << "with linespoints title 'Ideal' lw 2, \\\n";
        script << "     '" << folder2 << "/" << fileName << "' using 1:($2**2+$3**2) "
               << "with linespoints title 'Resist ' lw 2\n";
        script.close();
        std::cout << "Скрипт gnuplot создан: " << scriptName << std::endl;

        // Попытка запустить gnuplot (должен быть установлен в системе)
        system("gnuplot plot_script.gp");
        std::cout << "График сохранен в файл comparison_plot.png" << std::endl;
    } else {
        std::cerr << "Ошибка создания файла скрипта!" << std::endl;
        return 1;
    }

    return 0;
}