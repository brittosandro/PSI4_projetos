import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from distutils.spawn import find_executable
import seaborn as sns
import os
import numpy as np


# Esse script cálcula as constante de van der Waals
# a partir de potenciais analíticos derivados de interações
# entre moléculas polares e apolares.
# O script também calcula o potencial considerando o modelo de
# esferas rígidas.


def C_keeson(T):
    '''
    Essa função calcula a constante de Keeson.
    A função retorna o valor da constante de Keesnon (ck) em
    J*m**6 (joule x metros a sexta).
    
    Unidades dos Parâmetros
    ----------------------
    
    - Momento de dipolo (mu): C*m (coulomb x metro).
    - epsilon_0: C**2/N*m**2.
    - constante de Boltzmann (kb): J/K.
    - Temperatura (T): kelvin (K).
    - Constante dielétrica (k): adimencional.
    '''
    
    pi = 3.1416
    epsilon_0 = 8.85e-12
    k = 1
    kb = 1.38e-23
    dipolo1 = 1.8526
    dipolo2 = 1.8526
    mu1 = dipolo1 * 3.33e-30
    mu2 = dipolo2 * 3.33e-30
    ck = (1 / 3) * ((1 / (4 * pi * epsilon_0 * k))**2) * (1 / (kb * T)) * (mu1 * mu2)**2
    return ck

def C_debye():
    '''
    Essa função calcula a constante de Debye.
    A função retorna o valor da constante de Debye (cd)
    em J*m**6 (joule x metros a sexta).
    
    Unidades dos Parâmetros
    ----------------------
    
    - Momento de dipolo (mu): C*m (coulomb x metro).
    - epsilon_0: C**2/N*m**2.
    - Constante dielétrica (k): adimencional.
    - A polarizabilidade volumetrica (alfa_0): m**3.
    '''
    
    pi = 3.1416
    epsilon_0 = 8.85e-12
    k = 1
    dipolo1 = 1.8526
    # dipolo da molécula 1
    mu1 = dipolo1 * 3.33e-30
    # polarizabilidade da molécula 2
    alpha_0 = 0.848052e-30
    cb = 2 * (1 / (4 * pi * epsilon_0)) * (1 / k**2) * (mu1**2 * alpha_0)
    return cb

def C_london():
    '''
    Essa função calcula a constante de London.
    A função retorna o valor da constante de London (Cl)
    em J*m**6 (joule x metros a sexta).
    
    Unidades dos Parâmetros
    ----------------------
    
    A Polarizabilidade: m**3.
    A energia de inonização: eV.
    '''
    
    # Energia de ionização da molécula 1
    I1 = 11.4759
    # Energia de ionização da molécula 1
    I2 = 11.4759
    I_efe = (I1 * I2) / (I1 + I2)
    # Polarizabilidade da molécula 1
    alfa_01 = 0.848052e-30
    # Polarizabilidade da moléula 2
    alfa_02 = 0.848052e-30
    k = 1
    Cl = (3 / 2) * (1 / k**2) * (alfa_01 * alfa_02) * I_efe
    return Cl * 1.60218e-19

def C_vdW(T):
    '''
    Essa função calcula a constante de van der Waals.
    A função retorna o valor da constante em J*m**6 (joule x metros a sexta).
    '''
    TT = T
    cvdW = C_keeson(TT) + C_debye() + C_london()
    return cvdW

def U(T, r, r0):
    '''
    Potencial modelado como esferas rígidas.
    '''
    if r > r0:
        return - C_vdW(T) / (r**6)
    else:
        return 1.0e-20


if __name__ == '__main__':
    print(f'Constante de Keeson C_kesson = {C_keeson(300)}')
    print(f'Constante de Debye  C_debye = {C_debye()}')
    print(f'Constante de London C_london = {C_london()}')
    print("--------------------------------------------------")
    print(f'Constante de van der Waals C_vdW = {C_vdW(300)}')

    r0 = 2.4e-10                          # Disntância em metros (m)
    r = np.linspace(2.4e-10, 7e-10, 90)   # Disntâncias em metros (m)
    U_tot = [U(300, i, r0) for i in r]

    with open('dados_constantes.txt', 'w') as f:
        print(f'Esses são valores das constantes calculados com resultados,', end='\n', file=f)
        print(f'provenientes de cálculos teóricos.', end='\n', file=f)
        print("--------------------------------------------------------", end='\n', file=f)
        print(f'Constante de Keeson C_kesson = {C_keeson(300)}', end='\n', file=f)
        print(f'Constante de Debye  C_debye = {C_debye()}', end='\n', file=f)
        print(f'Constante de London C_london = {C_london()}', end='\n', file=f)
        print("--------------------------------------------------------", end='\n', file=f)
        print(f'Constante de van der Waals C_vdW = {C_vdW(300)}', end='\n', file=f)

    with open('calc_base_aug_cc_pvtz.txt', 'w') as f:
        print(f'# ------------------------------------------------------------- #', end='\n', file=f)
        print(f'#      Distância (r[m])           Energia Potencial (U[J*m**6]) #', end='\n', file=f)
        print(f'# ------------------------------------------------------------- #', end='\n', file=f)
        for i, j in zip(r, U_tot):
            print(f'{i:25.17}           {j:3.17}', end='\n', file=f)


    #reta = [0 for i in range(len(r))]

    fig, ax = plt.subplots()
    if find_executable('latex') and find_executable('dvipng'):
        mpl.rcParams.update({'font.size':18, 'text.usetex':True, 'font.family':
                           'serif', 'ytick.major.pad':4})
    else:
        mpl.rcParams.update({'font.size':18, 'font.family':'serif',
                           'ytick.major.pad':4})

    sns.set(style="ticks")

    fig1 = ax.plot(r, U_tot, lw=2.5, label=r"$U_{Total}$")
    #fig2 = ax.plot(r, reta, lw=2.5)

    fig.legend(loc='upper right', shadow=False, fontsize='large', bbox_to_anchor=(0.90, 0.82), frameon=False)
    fig.text(
          0.50,                      # Ordena posição x
          0.02,                      # Ordena posição y
          r"Distância ($m$)",
          ha = 'center',
          va = 'center',
          fontsize = 'xx-large')

    fig.text(
          0.03,
          0.5,
          r"$U$ ($J$)",
          ha = 'center',
          va = 'center',
          fontsize = 'xx-large',
          rotation = 'vertical')

    plt.savefig(os.path.splitext('fig_utotal_teorico1')[0] + ".pdf", bbox_inches='tight')
    plt.savefig(os.path.splitext('fig_utotal_teorico')[0] + ".png", dpi=300, orientation='portrait',
                transparent = True, format='png', bbox_inches = 'tight', pad_inches = .1)
