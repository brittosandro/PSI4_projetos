import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import seaborn as sns
import psi4


# distâncias de interesse em angstrons
dist = np.arange(3.159, 9.3, 0.2)

# listas de energias com zeros
eelst = np.zeros((len(dist)))
eind = np.zeros((len(dist)))
edisp = np.zeros((len(dist)))
eexch = np.zeros((len(dist)))
esapt = np.zeros((len(dist)))

for i in range(len(dist)):
    dimero_F2 = """
    0 1
    F                 0.000000    0.000000    0.698384
    F                 0.000000    0.000000   -0.698384
    --
    0 1
    F """ + str(dist[i]) + """  0.00000000    0.698384
    F """ + str(dist[i]) + """  0.00000000   -0.698384
    units angstrom
    symmetry c1
    """
    # constroi a molécula
    psi4.geometry(dimero_F2)

    # calcula a energia
    psi4.energy('sapt0/aug-cc-pvtz')

    # adionar as energias decompostas em suas respectivas listas
    # As energias são convertidas de Hartree para kcal/mol
    eelst[i] = psi4.variable('SAPT ELST ENERGY') * 627.509608
    eind[i] = psi4.variable('SAPT IND ENERGY') * 627.509608
    edisp[i] = psi4.variable('SAPT DISP ENERGY') * 627.509608
    eexch[i] = psi4.variable('SAPT EXCH ENERGY') * 627.509608
    esapt[i] = psi4.variable('SAPT TOTAL ENERGY') * 627.509608
    psi4.core.clean()

eelst1 = eelst
eind1 = eind
edisp1 = edisp
eexch1 = eexch
esapt1 = esapt

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

evdW = []
for i, j, k in zip(eelst1, eind1, edisp1):
    evdW.append(i+j+k)
#print('--------')
#print(evdW)

with open('dados_para_ajustar.txt', 'w') as f:
    for i, j in zip(dist, evdW):
        print(f'{i}            {j}', end='\n', file=f)

'''
fig1 = ax.plot(dist, eelst1, 'r', linestyle='-', linewidth=2.6, label='Eletrostática')
fig2 = ax.plot(dist, eind1, 'g', linestyle='-', linewidth=2.6, label='Indução')
fig3 = ax.plot(dist, edisp1, 'm', linestyle='-', linewidth=2.6, label='Disperção')
#fig4 = ax.plot(dist, eexch, 'bo', linestyle='-', linewidth=2.6, label='SAPT0 exch')
fig4 = ax.plot(dist, evdW, 'b', linestyle='-', linewidth=2.6, label='Uvdw')
fig5 = ax.plot(dist, esapt1, 'k', linestyle='-', linewidth=2.6, label='SAPT0 total')
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
             'scan_total_kcalmol.pdf',
             dpi=300,
             orientation = 'portrait',
             transparent = True,)
#plt.show()
'''

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
             'scan_ajuste_psiAPI.pdf',
             dpi=300,
             orientation = 'portrait',
             transparent = True,)
print(f'O valor de C_vdW Ótimo é {const[0]}')
#plt.savefig('ajuste.pdf', dpi=300, orientation = 'portrait', transparent = True,)
