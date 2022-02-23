import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import seaborn as sns

dados = np.loadtxt('dados_para_ajustar.txt')
dist = dados[:, 0]
evdW = dados[:, 1]

def func(dist, c, n):
    return - (c / dist) ** 6 + n

const , pcov = curve_fit(func, dist, evdW, [1, 5], maxfev=1000)

print(const[0])
print(const[1])
#$print(const[2])
#print(const[1])
print(pcov)
print(dist)

plt.close()
sns.set(style="ticks")
fig, ax = plt.subplots(figsize=(8, 6))

# Ajusta subplots.
fig.subplots_adjust(
                      left = 0.12,
                      right = 0.90,    # {Define as distâncias entre os extremos}
                      bottom = 0.12,
                      top = 0.93,
                      hspace = 0.24,   # Organiza espaçoes entre os subplots
                      wspace = 0.23    # Organiza espaçoes entre os subplots
                   )

fig11 = ax.plot(dist, evdW, 'o', linewidth=2.6, label="EvdW")
fig12 = ax.plot(dist, func(dist, *const), '-g', linewidth=2.6, label="Ajuste")
fig.legend(loc='upper right', shadow=False, fontsize='large',
           bbox_to_anchor=(0.88, 0.78), frameon=False)
fig.suptitle('Interação F2---F2', fontsize=16)
plt.grid(False)
# Descritores dos eixos
fig.text(
          0.04,                      # Ordena posição x
          0.5,                       # Ordena posição y
          r'U (kcal/mol)',
          ha = 'center',
          va = 'center',
          fontsize = 'xx-large',
          rotation = 'vertical')

fig.text(
          0.5,                      # Ordena posição x
          0.03,                     # Ordena posição y
          r'Distância $(\AA)$',
          ha = 'center',
          va = 'center',
          fontsize = 'xx-large')

plt.savefig(
             'scan_ajuste_new_1.pdf',
             dpi=300,
             orientation = 'portrait',
             transparent = True,)
print(f'O valor de C_vdW Ótimo é {const[0]}')
#plt.savefig('ajuste.pdf', dpi=300, orientation = 'portrait', transparent = True,)
