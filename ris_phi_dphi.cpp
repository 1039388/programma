#include <iostream>
#include <fstream>
#include <string>

void createPlot(const std::string& scriptName, 
                const std::string& outputImage, 
                const std::string& title, 
                const std::string& reFile, 
                const std::string& imFile, 
                const std::string& ylabel) {
    
    std::ofstream script(scriptName);
    if (script.is_open()) {
        script << "set terminal pngcairo size 800,600 enhanced font 'Verdana,10'\n";
        script << "set output '" << outputImage << "'\n";
        script << "set title '" << title << "'\n";
        script << "set xlabel 'z'\n";
        script << "set ylabel '" << ylabel << "'\n";
        script << "set grid\n";
        script << "set key right top\n";

        // Рисуем две колонки данных из разных файлов
        script << "plot '" << reFile << "' using 1:2 with lines title 'Re' lw 2 lc rgb 'blue', \\\n";
        script << "     '" << imFile << "' using 1:2 with lines title 'Im' lw 2 lc rgb 'red'\n";

        script.close();
        std::string command = "gnuplot " + scriptName;
        system(command.c_str());
        std::cout << "Создан график: " << outputImage << std::endl;
    }
}

int main() {
    // График для phi (Re_phi.txt и Im_phi.txt)
    createPlot("plot_phi.gp", "phi_plot.png", "Комплексный потенциал {/Symbol f}", 
               "Re_phi.txt", "Im_phi.txt", "phi");

    // График для dphi (Re_dphi.txt и Im_dphi.txt)
    createPlot("plot_dphi.gp", "dphi_plot.png", "Производная потенциала d{/Symbol f}/dz", 
               "Re_dphi.txt", "Im_dphi.txt", "dphi/dz");

    return 0;
}