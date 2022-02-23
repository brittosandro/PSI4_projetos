from constantes_de_van_der_Waals import C_vdW
import matplotlib.pyplot as plt
import numpy as np

# Esse módulo calcula as constantes teóricas da equação de estado de
# van der Waal a e b.

def b(R0):
    '''
    Essa função calcula a constante b da equação de van der Waal de maneira
    teórica, a partir de cálculos ab initio do volume molar para encontrarmos
    a distância característica.
    A função retorna esse parâmetro em m**3/mol (metros cúbicos por mol).

    - A distância característica (R0) é dada em metros.
    '''
    # Constante de avogadro
    Na = 6.02214076e23
    pi = 3.1416
    return (2 / 3) * (pi * Na) * (R0**3)

def a(T, R0):
    '''
    Essa função calcula a constante a da equação de van der Waals de maneira
    teórica, a partir de cálculos ab initio do volume molar para encontrarmos
    a distância característica.
    A função retorna esse parâmetro em J*m**3/mol**2.

    - A distância característica (R0) é dada em metros.
    '''
    # Constante de avogadro
    Na = 6.02214076e23
    pi = 3.1416
    return (2 / 3) * (pi * Na**2) * C_vdW(T) * (1 / (R0**3))

def a_sapt(C_vdW, R0):
    '''
    Essa função calcula a constante 'a' da equação de van der Waals de maneira
    teórica, a partir de cálculos ab initio SAPT. Ajustamos o valor de C_vdW
    e substituimos na fórmula para encontrarmos o valor de a.
    A função retorna esse parâmetro em J*m**3/mol**2.

    - A distância característica (R0) é dada em metros.
    '''
    Na = 6.02214076e23
    pi = 3.1416
    return (2 / 3) * (pi * Na**2) * C_vdW * (1 / (R0**3))

def segundo_coef_virial(a, b, T):
    '''
    Essa função calcula o segundo coeficiente do virial. Seus parâmetros
    são: a constante de van der Waals 'a', a constante de van der Waals 'b' bem
    como a temperatura T. A função retorna o valor do segundo coefiente do
    virial em m**3/mol.

    As unidades de medida de cada um dos parâmentros são:
    - a : J*m**3/mol**2
    - b : m**3/mol
    - kb: J/K
    - Na: mol**-1
    '''
    Na = 6.02214076e23
    kb = 1.380649e-23
    return b - (a / (Na * kb * T))

def segundo_coef_virial_sapt(a_sapt, b, T):
    '''
    Essa função calcula o segundo coeficiente do virial. Seus parâmetros
    são: a constante de van der Waals 'a_sapt' (essa constante vem de um
    cálculo SAPT a partir de um ajuste), a constante de van der Waals 'b' bem
    como a temperatura T. A função retorna o valor do segundo coefiente do
    virial em m**3/mol.

    As unidades de medida de cada um dos parâmentros são:
    - a : J*m**3/mol**2
    - b : m**3/mol
    - kb: J/K
    - Na: mol**-1
    '''
    Na = 6.02214076e23
    kb = 1.380649e-23
    return b - (a_sapt / (Na * kb * T))


#def plot_seg_coef_virial():
#    '''
#    Essa função plota os valores do coeficiente do virial para um conjunto de
#    temperaturas. A função recebe os valores 'a', 'b' e o intervalo de tempera_
#    ras que deverão ser descritos no gráfico.
#    '''


if __name__ == '__main__':
    #a = a(300, 3.4e-10)
    #b = b(3.4e-10)
    #scv = segundo_coef_virial(a, b, 300)
    #print(f'a = {a}')
    #print(f'b = {b}')
    #print(f'ScV = {scv}')
    R0 = 2.4e-10
    CvdW_sapt = 0.013e12
    Temperaturas = np.arange(300.0, 582.0, 12)
    b = b(R0)

    scv = [segundo_coef_virial(a(T, R0), b, T) for T in Temperaturas]
    with open('coef_virial.txt', 'w') as f:
        print(f'# ------------------------------------------------------------- #', end='\n', file=f)
        print(f'#      Temperatura([K])           Segundo Coef Virial([J*m**6]) #', end='\n', file=f)
        print(f'# ------------------------------------------------------------- #', end='\n', file=f)
        for i, j in zip(Temperaturas, scv):
            print(f'{i:15.17}                    {j:19.17}', end='\n', file=f)

    scv_sapt = [segundo_coef_virial(a(CvdW_sapt, R0), b, T) for T in Temperaturas]
    with open('coef_virial_SAPT.txt', 'w') as f:
        print(f'# ------------------------------------------------------------- #', end='\n', file=f)
        print(f'#      Temperatura([K])           Segundo Coef Virial([J*m**6]) #', end='\n', file=f)
        print(f'# ------------------------------------------------------------- #', end='\n', file=f)
        for i, j in zip(Temperaturas, scv_sapt):
            print(f'{i:15.17}                    {j:19.17}', end='\n', file=f)

    plt.plot(Temperaturas, scv)
    plt.plot(Temperaturas, scv_sapt)
    plt.savefig("scv_e_scv_sapt.png")
