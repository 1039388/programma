# programma
 ω=r*e^(i phase)
# Untitled-1.cpp - программа
reshatel(Model& model, Wall& wall, int resuis, double dx,double r_init, double phase_init,int bc_reg, int mod_switch,std::string imya) - функция, решающая уравнение Лодестро.
В качестве результата выдаёт std::pair<double,double> - r и phase частоты, а также сразу создаёт папку с некоторым именем и сохраняет результат.<br>
model - структура модели плазмы (struct PlasmaModel_A1, PlasmaModel_A3; struct PlasmaModel_A1_int - с интерполированным вакуумным магнитным полем). Туда же будет входить проводимость стенки (Dzeta)<br>
wall - структура стенки (struct Wall_St, struct Wall_Pr(model) ...)<br>
resuis - был сделан, чтобы переключать название папок.<br>
 double dx - шаг интегрирования (чаще всего хватало 0.001, но надо проверять)<br>
 double r_init, double phase_init - начальное приближение<br>
 int bc_reg - переключатель начального условия bc( bc_reg=0 => bc=phi'(L)=0; bc_reg!=0 => bc=phi(L)=0)<br>
если mode_switch=1 - пристрелка только по радиусу c фиксированным phase<br>
std::string imya - название папки<br>
reshatel_0 - решение уравнения при конкретной ω, входит в reshatel<br>
уравнение решается с помощью runge_kutta_fehlberg78 - в модельных случаях на скорость вычислений не влияет;<br>
запуск функции будет выглядеть примерно так:
para1=reshatel<LoDestroEquation>(plasma1, walstr,0,0.01, para1.first,para1.second,0,0,"probaA3");<br>
struct LoDestroEquation - структура с уравнением Лодестро (LoDestroEquation_int2 - для интерполированного поля, совместима только с PlasmaModel_A1_int)<br>
В решателе сидит костыль, выделяющий только неустойчивые ветви.(изредка отламывается)<br>
Помимо этого, внутри reshatel() есть отдельные параметры: r_step_min, phase_step_min - ограничивают бесконечные уменьшения шага, delta - задаёт точность зануления на правой границе( на практике чаще достигался минимальный шаг пристрелки, а на правой границе с достаточной точностью был 0)<br>

Для запуска потребуется библиотека Boost math.<br>
Исчерпывающий пример находится в main файла Untitled-1.cpp.
<br>
Дополнительно есть beta_omega.exe - если запустить его в папке с папками \beta, то создастся txt файл со списком Re(w) Im(w) beta.
