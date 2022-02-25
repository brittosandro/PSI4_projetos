import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
#mpl.use('Agg')
#from distutils.spawn import find_executable
import seaborn as sns
import os

def U_fixo(z, teta):
    Q = 1
    mu = 1.8 * (3.336 * 10 ** (-20))
    k = 1
    k_0 = 14.339 / (1.6 * 10 ** (-19))
    return (k_0 * mu * np.cos(teta)) * (1 / (k * z ** 2)) * 96.485307499

def U_medio(z, T):
    Q = 1
    mu = 1.8 * (3.336 * 10 ** (-20))
    k = 1
    k_0 = 14.339 / (1.6 * 10 ** (-19))
    kb = 8.617 * 10 ** (-5)
    return (- ((1 / (3 * kb * T)) * ((k_0 * mu * Q) / k) ** 2) * (1 / (z ** 4))) * 96.485307499

def U_induzido(z):
    '''
    As unidade que estamos considerando para calcular o potencial são:

    - carga Q em coulombs -> Q = n * e, em que e = 1.602176 * 10 ** (-19) C
    - kapa k não apresenta unidades
    - constante de coulomb k_a = 14.339 eV A(angstron) / (1.602176 * 10 ** (-19) C) ** 2
    - alpha zero alpha_0 = A^3 (angstron)^3

    O resultado do potencial é em kJ / mol, pois estamos multiplicando por 96.485307499
    '''

    Q = 1 * (1.602176 * 10 ** (-19))
    k = 2.3
    k_a = 14.339 / (1.602176 * 10 ** (-19)) ** 2
    alpha_0 = 4.8
    return - (k_a * (Q ** 2) * alpha_0) * (1 / k ** 2) * (1 / z ** 4) * 96.485307499

def U_keeson(z, T):
    '''
    Essa função calcula o potencial da energia eletrostática média de Keeson.
    Essa é uma interação do tipo dipolo - dipolo.
    A função retorna o valor do potencial em eV.
    OBSERVAÇÃO: Caso queira o valor do potencial em kJ/mol, multiplique Uk por
    96.485307499.
    As unidades consideradas são:

    - energia potencial eletrostática média (Uk) -> eV
    - Constante de Boltzmann (kb) -> É escrita em eV/K (eletron-volts por kelvin)
    - Temperatura (T) -> K (kelvin)
    - Constante de coulomb (ka) -> É escrita em eV*A/e**2 (e=1.60217662×10**-19C)
    - Momento de dipolo (mu)-> É escrito em C*A (coulomb * angstron)
    - Distância (z) -> É escrita em A (anstron)
    - Constante dielétrica (k) -> adimencional [Esse valore você deve alterar
    segundo seu sistema]
    '''
    kb = 8.617333262 * 10**(-5)
    ka = 14.339 / (1.60217662 * 10**(-19))**2
    dipolo1 = 1.7
    dipolo2 = 1.7
    k = 24
    mu1 = dipolo1 * (3.336 * 10**(-20))
    mu2 = dipolo2 * (3.336 * 10**(-20))
    Uk = - (1 / (3 * kb * T)) * (ka**2) * (1 / k**2) * ((mu1 * mu2)**2) * (1 / (z**6))
    return Uk

def U_debye(z):
    '''
    Essa função calcula o potencial de interação de Debye.
    Essa interação é do tipo dipolo - dipolo induzido.
    A função retorna o valor do potencial em eV.
    OBSERVAÇÃO: Caso queira o valor do potencial em kJ/mol, multiplique Ud por
    96.485307499.
    As unidades consideradas são:

    - Energia potencial de indução (Ud) -> eV.
    - Constante de coulomb (ka) -> É escrita em eV*A/e**2 (e=1.60217662×10**-19C)
    - Momento de dipolo (mu)-> É escrito em C*A (coulomb * angstron).
    - Polarizabilidade volumetrica (alfa_0) -> A**3 (angstrons cúbicos). Lembre-se
    que alfa = 4*pi*epsilon_0*alfa_0.
    - Distância (z) -> É escrita em A (anstron).
    - Constante dielétrica (k) -> adimencional [Esse valore você deve alterar
    segundo seu sistema].
    '''
    ka = 14.339 / (1.60217662 * 10**(-19))**2
    dipolo1 = 1.7
    alfa_0 = 5.1
    mu1 = dipolo1 * (3.336 * 10**(-20))
    k = 24
    Ub = -2 * ka * (1 / k**2) * (mu1**2 * alfa_0) * (1 / (z**6))
    return Ub

def U_london(z):
    '''
    Essa função calcula o potencial de interação de London.
    Essa interação é do tipo dipolo induzido - dipolo induzido.
    A função retorna o valor do potencial em eV.
    OBSERVAÇÃO: Caso queira o valor do potencial em kJ/mol, multiplique Ud por
    96.485307499.
    As unidades consideradas são:

    - Energia potencial de indução (Ul) -> eV.
    - Polarizabilidade volumetrica (alfa_0) -> A**3 (angstrons cúbicos). Lembre-se
    que alfa1 = 4*pi*epsilon_0*alfa_01 e alfa2 = 4*pi*epsilon_0*alfa_02.
    - Energia de ionização (I_efe)-> eV.
    - Distância (z) -> É escrita em A (anstron).
    - Constante dielétrica (k) -> adimencional [Esse valore você deve alterar
    segundo seu sistema].
    '''
    I1 = 10.2
    I2 = 10.2
    I_efe = (I1 * I2) / (I1 + I2)
    alfa_01 = 2.30
    alfa_02 = 2.30
    k = 1
    Ul = -(3 / 2) * (1 / k**2) * (alfa_01 * alfa_02) * I_efe * (1 / (z**6))
    return Ul

def C_keeson(T):
    '''
    Essa função calcula a constante de Keeson.
    A função retorna o valor da constante ck em J*m**6 (joule x metros a sexta).
    - Momento de dipolo (mu)-> É escrito em C*m (coulomb x metro).
    - epsilon_0 -> C**2/N*m**2.
    - constante de Boltzmann (kb) -> J/K.
    - Temperatura (T) -> kelvin.
    - Constante dielétrica (k) -> adimencional.
    '''
    pi = 3.1416
    epsilon_0 = 8.85e-12
    k = 1
    kb = 1.38e-23
    dipolo1 = 1.46
    dipolo2 = 1.46
    mu1 = dipolo1 * 3.33e-30
    mu2 = dipolo2 * 3.33e-30
    ck = (1 / 3) * ((1 / (4 * pi * epsilon_0 * k))**2) * (1 / (kb * T)) * (mu1 * mu2)**2
    return ck

def C_debye():
    '''
    Essa função calcula a constante de Debye.
    A função retorna o valor da constante cd em J*m**6 (joule x metros a sexta).
    - Momento de dipolo (mu)-> É escrito em C*m (coulomb x metro).
    - epsilon_0 -> C**2/N*m**2.
    - Constante dielétrica (k) -> adimencional.
    - A polarizabilidade volumetrica (alfa_0) -> m**3.
    '''
    pi = 3.1416
    epsilon_0 = 8.85e-12
    k = 1
    dipolo1 = 1.46
    mu1 = dipolo1 * 3.33e-30
    alpha_0 = 2.3e-30
    cb = 2 * (1 / (4 * pi * epsilon_0)) * (1 / k**2) * (mu1**2 * alpha_0)
    return cb

def C_london():
    '''
    Essa função calcula a constante de London.
    A função retorna o valor da constante em J*m**6 (joule x metros a sexta).
    A Polarizabilidade é descrita em m**3.
    A energia de inonização é descrita em eV.
    '''
    I1 = 10.2
    I2 = 10.2
    I_efe = (I1 * I2) / (I1 + I2)
    alfa_01 = 2.3e-30
    alfa_02 = 2.3e-30
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


#print(U_keeson(20, 300))
#print(U_debye(20))
#print(C_keeson(300))
#print(C_debye())
#print(C_london())
print(C_vdW(300))
#z = np.linspace(2, 20, 90)
#y = [U_keeson(i, 30000) for i in z]

#plt.plot(z, y)
#plt.show()



'''
z = np.linspace(2.5, 10, 90)
y0 = [U_induzido(i) for i in z]

print(U_induzido(20.0))


# Essa interação é do tipo 1: y - HF ----- Li
#interacao1 = np.loadtxt("interacao1.txt")
#z_int1 = interacao1[:, 0]
#u_int1 = interacao1[:, 1]

fig, ax = plt.subplots()

if find_executable('latex') and find_executable('dvipng'):
    mpl.rcParams.update({'font.size':18, 'text.usetex':True, 'font.family':
                           'serif', 'ytick.major.pad':4})
else:
    mpl.rcParams.update({'font.size':18, 'font.family':'serif',
                           'ytick.major.pad':4})

sns.set(style="ticks")

#fig1 = ax.plot(z, y0, lw=2.5, label=r"$U_{inst}(180^{\circ})$")
#fig1 = ax.plot(z, y000, lw=2.5, label=r"$U_{inst}(90^{\circ})$")
#fig11 = ax.plot(z, y00, lw=2.5, label=r"$U_{inst}(0^{\circ})$")
#fig2 = ax.plot(z, y1, lw=2.5, label=r"$U_{medio}(300K)$")
#fig3 = ax.plot(z, y2, lw=2.5, label=r"$U_{medio}(z, 30000K)$")
#fig4 = ax.plot(z_int1, u_int1, lw=2.5, label=r"$U_{sitio 1}$")
#fig5 = ax.plot(z_int2, u_int2, lw=2.5, label=r"$U_{sitio 2}$")
#fig6 = ax.plot(z_int3, u_int3, lw=2.5, label=r"$U_{sitio 3}$")
#plt.yticks([-80, -60, -40, -20, 0])


fig.legend(loc='upper right', shadow=False, fontsize='large', bbox_to_anchor=(0.90, 0.82), frameon=False)


fig.text(
          0.50,                      # Ordena posição x
          0.02,                      # Ordena posição y
          r"Distância ($\textrm{\AA}$)",
          ha = 'center',
          va = 'center',
          fontsize = 'xx-large')

fig.text(
          0.03,
          0.5,
          r"$U (kJ/mol)$",
          ha = 'center',
          va = 'center',
          fontsize = 'xx-large',
          rotation = 'vertical')

#plt.show()

plt.savefig(os.path.splitext('fig_umedio_uinst1')[0] + ".pdf", bbox_inches='tight')
plt.savefig(os.path.splitext('fig_umedio_uinst2')[0] + ".png", dpi=300, orientation='portrait',
            transparent = True, format='png', bbox_inches = 'tight', pad_inches = .1)
'''
