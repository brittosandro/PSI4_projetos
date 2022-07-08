from glob import glob

import numpy as np
import psi4
import re
import subprocess


psi4.set_memory('10 GB')
print('\n')
psi4.core.set_output_file('output.dat', False)

numpy_memory = 10


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


geometrias_amonia = glob('*_sapt.xyz')
#print(geometrias_amonia)

#metodos = ['sapt0', 'sapt2', 'sapt2+']
#metodos = ['sapt0', 'sapt2', 'sapt2+', 'sapt2+(3)', 'sapt2+3', ]
#metodos = ['ccsd(t)']
metodos = ['ccsd(t)']
bases = ['jun-cc-pvdz']
#gases_nobres = ['He', 'Ne', 'Ar', 'Kr']
gases_nobres = ['He', 'Ne']

for geo in geometrias_amonia:
    with open(geo, 'r') as f:
        str_geo = f.read()

    for gas_nobre in gases_nobres:
        for metodo in metodos:
            for base in bases:
                proc1 = subprocess.run(['mkdir', f'{gas_nobre}_{metodo}_{base}'])

                distancias = np.arange(3.5, 11.1, 0.5)
                # listas de energias com zeros
                eccsd = np.zeros(((len(distancias)), 3))
                eccsdt = np.zeros(((len(distancias)), 3))

                for i, dist in enumerate(distancias):
                    # Construindo a geometria do Dimero
                    dimero = """
                    0 1 """ + """\n""" + str_geo + """--
                    0 1 """ + """\n""" + gas_nobre + """  """ + str(distancias[i]) + """  0.0000000      0.0000000000
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
                    eccsd[i, 0] = psi4.variable('CCSD TOTAL ENERGY') * 27211.4
                    eccsdt[i, 0] = psi4.variable('CCSD(T) TOTAL ENERGY') * 27211.4

                    psi4.core.clean()

                    # Construindo a geometria do monomeroA
                    monomeroA = """
                    0 1 """ + """\n""" + str_geo + """--
                    0 1 """ + """\n""" + f'Gh({gas_nobre})' + """  """ + str(distancias[i]) + """  0.0000000      0.0000000000
                    units angstrom
                    symmetry c1
                    """
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
                    monomeroB = """
                    0 1 """ + """\n""" + new_str_geo + """--
                    0 1 """ + """\n""" + gas_nobre + """  """ + str(distancias[i]) + """  0.0000000      0.0000000000
                    units angstrom
                    symmetry c1
                    """
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

                #time.sleep(2)
                #with open(nome_arq_out, 'w') as f:
                #    print(f'# Dist    Energia Eletrostatica       Energia Inducao        Energia Dispercao      Energia EXCH       Energia Sapt        Energia CT', end='\n', file=f)
                #    for d, el, ei, ed, eex, es in zip(distancias, np.around(eelst, 7), np.around(eind, 7), np.around(edisp, 7),
                #                                   np.around(eexch, 7), np.around(esapt, 7)):
                #        print(f'{d:5.2f}  {el:16.7f}   {ei:25.7f}   {ed:20.7f}   {eex:20.7f}   {es:15.7f}', end='\n', file=f)

                with open(nome_arq_out, 'w') as f:
                    print(f'# Dist     Energia Total CCSD  Energia Total CCSD(T)', end='\n', file=f)
                    for d, eccsdint, eccsdtint in zip(distancias, np.around(Eint_ccsd, 7), np.around(Eint_ccsdt, 7)):
                        print(f'{d:5.2f} {eccsdint:16.7f} {eccsdtint:16.7f}', end='\n', file=f)

                proc2 = subprocess.run(['mv', nome_arq_out, f'{gas_nobre}_{metodo}_{base}'])
