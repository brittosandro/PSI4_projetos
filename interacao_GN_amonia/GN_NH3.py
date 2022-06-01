from glob import glob
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt
import numpy as np
import psi4
import seaborn as sns
import subprocess
import time


geometrias_amonia = glob('*_sapt.xyz')
#print(geometrias_amonia)

#metodos = ['sapt0', 'sapt2', 'sapt2+']
metodos = ['sapt0']
bases = ['jun-cc-pvdz']
gases_nobres = ['He', 'Ne', 'Ar', 'Kr']
#gases_nobres = ['He']

for geo in geometrias_amonia:
    with open(geo, 'r') as f:
        str_geo = f.read()

    for gas_nobre in gases_nobres:
        for metodo in metodos:
            for base in bases:
                proc1 = subprocess.run(['mkdir', f'{gas_nobre}_{metodo}_{base}'])

                distancias = np.arange(1.4, 10.1, 0.1)
                # listas de energias com zeros
                eelst = np.zeros((len(distancias)))
                eind = np.zeros((len(distancias)))
                edisp = np.zeros((len(distancias)))
                eexch = np.zeros((len(distancias)))
                esapt = np.zeros((len(distancias)))
                #ect = np.zeros((len(distancias)))

                for i, dist in enumerate(distancias):
                    # Construindo a geometria do Dimero
                    dimero = """
                    0 1 """ + """\n""" + str_geo + """--
                    0 1 """ + """\n""" + gas_nobre + """  """ + str(distancias[i]) + """     0.0000000      0.0000000000
                    units angstrom
                    symmetry c1
                    """

                    # constroi a molecula
                    psi4.geometry(dimero)

                    # calcula a energia
                    nivel = f'{metodo}/{base}'
                    psi4.energy(nivel)

                    # adiona as energias decompostas em suas respectivas listas
                    # As energias sao convertidas de Hartree para meV.
                    #ect[i] = psi4.variable('SAPT CT ENERGY') * 27211.4
                    eelst[i] = psi4.variable('SAPT ELST ENERGY') * 27211.4
                    eind[i] = psi4.variable('SAPT IND ENERGY') * 27211.4
                    edisp[i] = psi4.variable('SAPT DISP ENERGY') * 27211.4
                    eexch[i] = psi4.variable('SAPT EXCH ENERGY') * 27211.4
                    esapt[i] = psi4.variable('SAPT TOTAL ENERGY') * 27211.4

                    psi4.core.clean()

                sitio_inte = geo.replace('sapt.xyz', '').replace('amonia', '')
                nome_arq_out = metodo +  sitio_inte + base + '.dat'


                time.sleep(2)
                with open(nome_arq_out, 'w') as f:
                    print(f'# Dist    Energia Eletrostatica       Energia Inducao        Energia Dispercao      Energia EXCH       Energia Sapt        Energia CT', end='\n', file=f)
                    for d, el, ei, ed, eex, es in zip(distancias, np.around(eelst, 7), np.around(eind, 7), np.around(edisp, 7),
                                                   np.around(eexch, 7), np.around(esapt, 7)):
                        print(f'{d:5.2f}  {el:16.7f}   {ei:25.7f}   {ed:20.7f}   {eex:20.7f}   {es:15.7f}', end='\n', file=f)


                proc2 = subprocess.run(['mv', nome_arq_out, f'{gas_nobre}_{metodo}_{base}'])



'''
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

        fig1 = ax.plot(distancias, eelst, 'r+', linestyle='-', linewidth=2.6, label='Eletrostática')
        fig2 = ax.plot(distancias, eind, 'g^', linestyle='-', linewidth=2.6, label='Indução')
        fig3 = ax.plot(distancias, edisp, 'mx', linestyle='-', linewidth=2.6, label='Disperção')
        fig4 = ax.plot(distancias, eexch, 'bo', linestyle='-', linewidth=2.6, label='SAPT0 exch')
        fig5 = ax.plot(distancias, esapt, 'k*', linestyle='-', linewidth=2.6, label='SAPT0 total')
        fig.legend(loc='upper right', shadow=False, fontsize='large',
                   bbox_to_anchor=(0.88, 0.78), frameon=False)
        fig.suptitle('Interação He---NH3', fontsize=16)
        plt.grid(False)
        # Descritores dos eixos
        fig.text(
                 0.04,                      # Ordena posição x
                 0.5,                       # Ordena posição y
                 r'U (meV)',
                 ha = 'center',
                 va = 'center',
                 fontsize = 'xx-large',
                 rotation = 'vertical'
        )
        fig.text(
                 0.5,                      # Ordena posição x
                 0.03,                     # Ordena posição y
                 r'Distância $(\AA)$',
                 ha = 'center',
                 va = 'center',
                 fontsize = 'xx-large')

        nome = arq.replace('_sapt.xyz', '')
        plt.savefig(f'{nome}_{metodo}_{base}', dpi=300, orientation = 'portrait', transparent = True,)
        plt.savefig(f'{nome}_{metodo}_{base}', dpi=300, orientation = 'portrait', transparent = True, format='pdf')
'''

'''
# distâncias de interesse em angstrons
dist = np.arange(1.62781890, 2.3, 0.2)
dist_F = np.arange(2.50721912, )

# listas de energias com zeros
eelst = np.zeros((len(dist)))
eind = np.zeros((len(dist)))
edisp = np.zeros((len(dist)))
eexch = np.zeros((len(dist)))
esapt = np.zeros((len(dist)))


for i in range(len(dist)):
    dimero_HF = """
    0 1
    N       0.0000000000     0.0000000000     0.0000000000
    H       0.9999999978    -0.0000000000    -0.0000000000
    H      -0.3333333348     0.9428090383    -0.0000000000
    H      -0.3333333355    -0.4714045245    -0.8164965751
    --
    0 1
    He """+ str(dist[i]) + """ 0.0000000      0.0000000000
    units angstrom
    symmetry c1
    """
    # constroi a molécula
    psi4.geometry(dimero_HF)

    # calcula a energia
    psi4.energy('sapt0/jun-cc-pvdz')

    # adionar as energias decompostas em suas respectivas listas
    # As energias são convertidas de Hartree para kcal/mol.
    eelst[i] = psi4.variable('SAPT ELST ENERGY') * 27211.4
    eind[i] = psi4.variable('SAPT IND ENERGY') * 27211.4
    edisp[i] = psi4.variable('SAPT DISP ENERGY') * 27211.4
    eexch[i] = psi4.variable('SAPT EXCH ENERGY') * 27211.4
    esapt[i] = psi4.variable('SAPT TOTAL ENERGY') * 27211.4
    psi4.core.clean()


print(eelst)
print('----')
print(eind)
print('----')
print(edisp)


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

fig1 = ax.plot(dist, eelst, 'r', linestyle='-', linewidth=2.6, label='Eletrostática')
fig2 = ax.plot(dist, eind, 'g', linestyle='-', linewidth=2.6, label='Indução')
fig3 = ax.plot(dist, edisp, 'm', linestyle='-', linewidth=2.6, label='Disperção')
#fig4 = ax.plot(dist, eexch, 'bo', linestyle='-', linewidth=2.6, label='SAPT0 exch')
fig5 = ax.plot(dist, esapt, 'k', linestyle='-', linewidth=2.6, label='SAPT0 total')
fig.legend(loc='upper right', shadow=False, fontsize='large',
           bbox_to_anchor=(0.88, 0.78), frameon=False)
fig.suptitle('Interação He---NH3 (bond)', fontsize=16)
plt.grid(False)
# Descritores dos eixos
fig.text(
          0.04,                      # Ordena posição x
          0.5,                       # Ordena posição y
          r'U (meV)',
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

plt.savefig('scan_psiAPI.pdf', dpi=300, orientation = 'portrait', transparent = True,)

evdW = []
for i, j, k in zip(eelst, eind, edisp):
    evdW.append(i+j+k)
print('--------')
print(evdW)

def func(dist, c):
    return - (c / dist) ** 6

const , pcov = curve_fit(func, dist, evdW, [1], maxfev=1000)

print(const[0])
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
fig12 = ax.plot(dist, func(dist, const), '-g', linewidth=2.6, label="aJUSTE")
fig.legend(loc='upper right', shadow=False, fontsize='large',
           bbox_to_anchor=(0.88, 0.78), frameon=False)
fig.suptitle('Interação HF---HF', fontsize=16)
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
             'scan_testando_psiAPI.pdf',
             dpi=300,
             orientation = 'portrait',
             transparent = True,)
print(f'O valor de C_vdW Ótimo é {const[0]}')
plt.savefig('ajuste.pdf', dpi=300, orientation = 'portrait', transparent = True,)

#plt.show()
'''
