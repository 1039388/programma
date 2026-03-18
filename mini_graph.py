import matplotlib.pyplot as plt
import numpy as np
# Ваши массивы
x = [1, 2, 3, 4]
y = [np.log(0.425-0.41875),np.log(0.41875-0.4186401367),np.log( 0.4186401367-0.4186309814),np.log( abs(0.4186309814-0.4186332703))]

# Построение графика
plt.plot(x, y, 'o-', color='purple', linewidth=2, markersize=8, 
         markerfacecolor='yellow', markeredgecolor='black', markeredgewidth=1)
plt.xlabel('k')
plt.ylabel('log(delta(w))')
plt.title('log(w^k-w^{k-1})(k)')
plt.grid(True)
plt.show()