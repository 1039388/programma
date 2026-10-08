# programma
 ω=r*e^(i phase)
# Untitled-1.cpp - программа
reshatel(Model& model, Wall& wall, int resuis, double dx,double r_init, double phase_init,int bc_reg, int mod_switch,std::string imya) - функция, решающая уравнение Лодестро.
В качестве результата выдаёт std::pair<double,double> - r и phase частоты, а также сразу создаёт папку с некоторым именем и сохраняет результат.
model - структура модели плазмы (struct PlasmaModel_A1, PlasmaModel_A3; struct PlasmaModel_A1_int - с интерполированным вакуумным магнитным полем). Туда же будет входить проводимость стенки (Dzeta)
wall - структура стенки (struct Wall_St, struct Wall_Pr(model) ...)
resuis - был сделан, чтобы переключать название папок.
 double dx - шаг интегрирования (чаще всего хватало 0.001, но надо проверять)
 double r_init, double phase_init - начальное приближение
 int bc_reg - переключатель начального условия bc( bc_reg=0 => bc=phi'(L)=0; bc_reg!=0 => bc=phi(L)=0)
если mode_switch=1 - пристрелка только по радиусу c фиксированным phase
std::string imya - название папки
reshatel_0 - решение уравнения при конкретной ω, входит в reshatel
уравнение решается с помощью runge_kutta_fehlberg78 - в модельных случаях на скорость вычислений не влияет;
запуск функции будет выглядеть примерно так:
para1=reshatel<LoDestroEquation>(plasma1, walstr,0,0.01, para1.first,para1.second,0,0,"probaA3");
struct LoDestroEquation - структура с уравнением Лодестро (LoDestroEquation_int2 - для интерполированного поля, совместима только с PlasmaModel_A1_int)
В решателе сидит костыль, выделяющий только неустойчивые ветви.
Помимо этого, внутри reshatel() есть отдельные параметры: r_step_min, phase_step_min - ограничивают бесконечные уменьшения шага, delta - задаёт точность зануления на правой границе( на практике чаще достигался минимальный шаг пристрелки, а на правой границе с достаточной точностью был 0)

Для запуска потребуется библиотека Boost math.
Исчерпывающий пример находится в main файла Untitled-1.cpp.

Дополнительно есть beta_omega.exe - если запустить его в папке с папками beta, то создаст txt файл со списком Re(w) Im(w) beta.
