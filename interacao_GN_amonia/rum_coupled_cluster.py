from glob import glob

import numpy as np
import psi4
import re
import subprocess


psi4.set_memory('10 GB')
print('\n')
psi4.core.set_output_file('output.dat', False)

numpy_memory = 10

def cria_diretorio(gas, metodo, base):
    '''
    Essa função recebe três strings que correspondem ao gás nobre de interesse
    o método e a base que iremos executar os cálculos e cria um novo diretorio
    com as respectivas palavras.
    '''
    subprocess.run(['mkdir', f'{gas_nobre}_{metodo}_{base}'])

def cria_matriz(d):
    '''
    Essa função recebe as distâncias em angstrom que correspondem ao
    intervalo do cálculo e retorna uma matriz (array) com um conjunto de linhas
    e três colunas cujo os elementos são zeros.
    '''
    return np.zeros(((len(d)), 3))

def input_geo(geo, gas, d, i):
    '''
    Como parâmetros de entrada a função recebe a geometria do composto (geo)
    os símbolos dos gases nobres (gas) as distancias de cada interação (d) e
    um contador i e retorna a geometria do input.
    '''

    input = """
    0 1 """ + """\n""" + geo + """--
    0 1 """ + """\n""" + gas + """  """ + str(d[i]) + """  0.000000   0.00000000
    units angstrom
    symmetry c1
    """
    return input


def casando(ghost):
    def substitui(m):
        text = m.group()
        if text == 'N':
            return ghost + '(N)'
        if text == 'H':
            return ghost + '(H)'
    return substitui

def Eint(matriz):
    '''
    Essa função recebe uma matriz com valores correspondetes as Energias
    dos cálculos de estrutura eletrônica. A matriz (array) é composta por
    3 colunas, das quais temos:
    Col1 = Energia do Dimero
    Col2 = Eneriga do monomero A
    Col3 = Energia do monomero B
    A função retorna a energia de interação (Eint) considerando a relação:
    Eint = EAB - (EA + EB)
    '''
    EAB = [matriz[i, 0] for i in range(len(matriz))]
    EA = [matriz[i, 1] for i in range(len(matriz))]
    EB = [matriz[i, 2] for i in range(len(matriz))]

    soma_AB = [i+j for i, j in zip(EA, EB)]
    Eint = [i-j for i, j in zip(EAB, soma_AB)]

    return Eint

def cria_arquivo(nome, dist, eint1, eint2):
    with open(nome, 'w') as f:
        print(f'# Dist     Energia Total CCSD  Energia Total CCSD(T)', end='\n', file=f)
        for d, e1, e2 in zip(dist, np.around(eint1, 7), np.around(eint2, 7)):
            print(f'{d:5.2f} {e1:16.7f} {e2:16.7f}', end='\n', file=f)

geometrias_amonia = glob('*_sapt.xyz')
#print(geometrias_amonia)

metodos = ['ccsd(t)']
bases = ['jun-cc-pvdz']
#gases_nobres = ['He', 'Ne', 'Ar', 'Kr']
gases_nobres = ['He', 'Ne']

for gas_nobre in gases_nobres:
    for metodo in metodos:
        for base in bases:
            cria_diretorio(gas_nobre, metodo, base)
            for geo in geometrias_amonia:
                with open(geo, 'r') as f:
                    str_geo = f.read()

                distancias = np.arange(3.5, 11.1, 0.5)
                eccsd = cria_matriz(distancias)
                eccsdt = cria_matriz(distancias)

                for i, dist in enumerate(distancias):
                    # Construindo a geometria do Dimero
                    dimero = input_geo(str_geo, gas_nobre, distancias, i)
                   # constroi a molecula
                    psi4.geometry(dimero)
                    # calcula a energia
                    nivel = f'{metodo}/{base}'
                    psi4.energy(nivel)
                    # adiona as energias decompostas em suas respectivas listas
                    # As energias sao convertidas de Hartree para meV.
                    eccsd[i, 0] = psi4.variable('CCSD TOTAL ENERGY') * 27211.4
                    eccsdt[i, 0] = psi4.variable('CCSD(T) TOTAL ENERGY') * 27211.4

                    psi4.core.clean()

                    # Construindo a geometria do monomeroA
                    monomeroA = input_geo(str_geo, f'Gh({gas_nobre})', distancias, i)
                    #constroi a molecula
                    psi4.geometry(monomeroA)
                    # calcula a energia
                    nivel = f'{metodo}/{base}'
                    psi4.energy(nivel)
                    eccsd[i, 1] = psi4.variable('CCSD TOTAL ENERGY') * 27211.4
                    eccsdt[i, 1] = psi4.variable('CCSD(T) TOTAL ENERGY') * 27211.4
                    psi4.core.clean()

                    # Construindo a geometria do monomeroB
                    new_str_geo = re.sub(r'[a-zA-Z]', casando('Gh'), str_geo,
                                         flags=re.IGNORECASE)
                    monomeroB = input_geo(new_str_geo, gas_nobre, distancias, i)
                    #constroi a molecula
                    psi4.geometry(monomeroB)
                    # calcula a energia
                    nivel = f'{metodo}/{base}'
                    psi4.energy(nivel)
                    eccsd[i, 2] = psi4.variable('CCSD TOTAL ENERGY') * 27211.4
                    eccsdt[i, 2] = psi4.variable('CCSD(T) TOTAL ENERGY') * 27211.4
                    psi4.core.clean()

                    Eint_ccsd = Eint(eccsd)
                    Eint_ccsdt = Eint(eccsdt)

                sitio_inte = geo.replace('sapt.xyz', '').replace('amonia', '')
                nome_arq_out = metodo +  sitio_inte + base + '.dat'

                cria_arquivo(nome_arq_out, distancias, Eint_ccsd, Eint_ccsdt)

                proc2 = subprocess.run(['mv', nome_arq_out, f'{gas_nobre}_{metodo}_{base}'])
