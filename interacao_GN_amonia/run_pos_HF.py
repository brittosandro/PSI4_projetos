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
    subprocess.run(['mkdir', f'{gas}_{metodo}_{base}'])

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

def cria_arquivo(nome, metodo, dist, eint1, eint2, eint3, eint4, eint5, eint6, eint7):
    if metodo == 'ccsd(t)':
        with open(nome, 'w') as f:
            print(f'# Distancia [angstrom]    |   Energia [meV]', end='\n', file=f)
            print(f'# -----------------------------------------', end='\n', file=f)
            print(f'# Dist   Eint CCSD(T)-NOCP  Eint CCSD(T)-CP', end='\n', file=f)
            print(f'# -----------------------------------------', end='\n', file=f)
            for d, e1, e2 in zip(dist, np.around(eint1, 6), np.around(eint2, 6)):
                print(f'{d:5.2f} {e1:16.9f} {e2:17.9f}', end='\n', file=f)
                
    if metodo == 'ccsd':
        with open(nome, 'w') as f:
            print(f'# Distancia [angstrom]    |   Energia [meV]', end='\n', file=f)
            print(f'# -----------------------------------------', end='\n', file=f)
            print(f'# Dist    Eint CCSD-NOCP    Eint CCSD-CP', end='\n', file=f)
            print(f'# -----------------------------------------', end='\n', file=f)
            for d, e1, e2 in zip(dist, np.around(eint1, 6), np.around(eint2, 6)):
                print(f'{d:6.2f} {e1:16.9f} {e2:17.9f}', end='\n', file=f)

    if metodo == 'mp2':
        with open(nome, 'w') as f:
            print(f'# Distancia [angstrom]    |   Energia [meV]', end='\n', file=f)
            print(f'# -----------------------------------------', end='\n', file=f)
            print(f'# Dist    Eint MP2-NOCP    Eint MP2-CP', end='\n', file=f)
            print(f'# -----------------------------------------', end='\n', file=f)
            for d, e1, e2 in zip(dist, np.around(eint1, 6), np.around(eint2, 6)):
                print(f'{d:6.2f} {e1:16.9f} {e2:17.9f}', end='\n', file=f)

    if metodo == 'mp4':
        with open(nome, 'w') as f:
            print(f'# Distancia [angstrom]    |   Energia [meV]', end='\n', file=f)
            print(f'# -----------------------------------------', end='\n', file=f)
            print(f'# Dist    Eint MP4-NOCP    Eint MP4-CP', end='\n', file=f)
            print(f'# -----------------------------------------', end='\n', file=f)
            for d, e1, e2 in zip(dist, np.around(eint1, 6), np.around(eint2, 6)):
                print(f'{d:6.2f} {e1:16.9f} {e2:15.9f}', end='\n', file=f)

    metodos_perturbativos = ('sapt0', 'sapt2', 'sapt2+', 'sapt2+(3)', 'sapt2+3')
    if metodo in metodos_perturbativos:
        with open(nome_arq_out, 'w') as f:
            print(f'#                        Distancia [angstrom]    |   Energia [meV]', end='\n', file=f)
            print(f'# ---------------------------------------------------------------------------------------------------------------', end='\n', file=f)
            print(f'# Dist    Energia Eletrostatica       Energia Inducao        Energia Dispercao      Energia EXCH     Energia Sapt', end='\n', file=f)
            print(f'# ---------------------------------------------------------------------------------------------------------------', end='\n', file=f)
            for d, el, ei, ed, eex, es in zip(dist, np.around(eint3, 7), np.around(eint4, 7),
                                              np.around(eint5, 7), np.around(eint6, 7), np.around(eint7, 7)):
                print(f'{d:6.2f}  {el:16.9f}   {ei:25.9f}   {ed:20.9f}   {eex:18.9f}   {es:13.9f}', end='\n', file=f)


def move_arquivo(nome, gas, metodo, base):
    subprocess.run(['mv', nome, f'{gas}_{metodo}_{base}'])


geometrias_amonia = glob('*_sapt.xyz')
#print(geometrias_amonia)

metodos = ['ccsd', 'ccsd(t)', 'mp2', 'mp4', 'sapt0',
           'sapt2', 'sapt2+', 'sapt2+(3)', 'sapt2+3']

bases = ['jun-cc-pvdz']
#gases_nobres = ['He', 'Ne', 'Ar', 'Kr']
gases_nobres = ['He',]

for gas_nobre in gases_nobres:
    for metodo in metodos:
        for base in bases:
            cria_diretorio(gas_nobre, metodo, base)
            nivel = f'{metodo}/{base}'
            for geo in geometrias_amonia:
                with open(geo, 'r') as f:
                    str_geo = f.read()

                distancias = np.arange(3.5, 4.1, 0.5)

                eccsd = cria_matriz(distancias)
                eccsdt = cria_matriz(distancias)

                en_sem_cp =  np.zeros((len(distancias)))
                en_com_cp = np.zeros((len(distancias)))

                eelst = np.zeros((len(distancias)))
                eind = np.zeros((len(distancias)))
                edisp = np.zeros((len(distancias)))
                eexch = np.zeros((len(distancias)))
                esapt = np.zeros((len(distancias)))

                for i, dist in enumerate(distancias):
                    # Construindo a geometria do Dimero
                    dimero = input_geo(str_geo, gas_nobre, distancias, i)
                    '''
                    #Essa é uma forma de expressar a energia do dimero considerando
                    #átomos fantasmas. É uma maneira ineficiente dado as formas
                    #internas que o PSI4 pode fazer.
                    # constroi a molecula
                    psi4.geometry(dimero)
                    # calcula a energia
                    psi4.energy(nivel)
                    # adiona as energias em suas respectivas matrizes
                    # As energias sao convertidas de Hartree para meV.
                    eccsd[i, 0] = psi4.variable('CCSD TOTAL ENERGY') * 27211.4
                    eccsdt[i, 0] = psi4.variable('CCSD(T) TOTAL ENERGY') * 27211.4
                    psi4.core.clean()

                    # Construindo a geometria do monomeroA
                    monomeroA = input_geo(str_geo, f'Gh({gas_nobre})', distancias, i)
                    #constroi a molecula
                    psi4.geometry(monomeroA)
                    # calcula a energia
                    psi4.energy(nivel)
                    eccsd[i, 1] = psi4.variable('CCSD TOTAL ENERGY') * 27211.4
                    eccsdt[i, 1] = psi4.variable('CCSD(T) TOTAL ENERGY') * 27211.4
                    psi4.core.clean()

                    # Construindo a geometria do monomeroB
                    Gh_str_geo = re.sub(r'[a-zA-Z]', casando('Gh'), str_geo,
                                         flags=re.IGNORECASE)
                    monomeroB = input_geo(Gh_str_geo, gas_nobre, distancias, i)
                    #constroi a molecula
                    psi4.geometry(monomeroB)
                    # calcula a energia
                    psi4.energy(nivel)
                    eccsd[i, 2] = psi4.variable('CCSD TOTAL ENERGY') * 27211.4
                    eccsdt[i, 2] = psi4.variable('CCSD(T) TOTAL ENERGY') * 27211.4
                    psi4.core.clean()
                    Eint_ccsd = Eint(eccsd)
                    Eint_ccsdt = Eint(eccsdt)
                    '''
                    if metodo == 'ccsd(t)':
                        psi4.geometry(dimero)
                        psi4.energy(nivel, bsse_type=['nocp', 'cp',])
                        en_sem_cp[i] = psi4.variable('NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY') * 27211.4
                        en_com_cp[i] = psi4.variable('CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY') * 27211.4
                        psi4.core.clean()
                    if metodo == 'ccsd':
                        psi4.geometry(dimero)
                        psi4.energy(nivel, bsse_type=['nocp', 'cp',])
                        en_sem_cp[i] = psi4.variable('NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY') * 27211.4
                        en_com_cp[i] = psi4.variable('CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY') * 27211.4
                        psi4.core.clean()
                    if metodo == 'mp2':
                        psi4.geometry(dimero)
                        psi4.energy(nivel, bsse_type=['nocp', 'cp',])
                        en_sem_cp[i] = psi4.variable('NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY') * 27211.4
                        en_com_cp[i] = psi4.variable('CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY') * 27211.4
                        psi4.core.clean()
                    if metodo == 'mp4':
                        psi4.geometry(dimero)
                        psi4.energy(nivel, bsse_type=['nocp', 'cp',])
                        en_sem_cp[i] = psi4.variable('NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY') * 27211.4
                        en_com_cp[i] = psi4.variable('CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY') * 27211.4
                        psi4.core.clean()
                    if metodo == 'sapt0':
                        psi4.geometry(dimero)
                        psi4.energy(nivel)
                        eelst[i] = psi4.variable('SAPT ELST ENERGY') * 27211.4
                        eind[i] = psi4.variable('SAPT IND ENERGY') * 27211.4
                        edisp[i] = psi4.variable('SAPT DISP ENERGY') * 27211.4
                        eexch[i] = psi4.variable('SAPT EXCH ENERGY') * 27211.4
                        esapt[i] = psi4.variable('SAPT TOTAL ENERGY') * 27211.4
                        psi4.core.clean()
                    if metodo == 'sapt2':
                        psi4.geometry(dimero)
                        psi4.energy(nivel)
                        eelst[i] = psi4.variable('SAPT ELST ENERGY') * 27211.4
                        eind[i] = psi4.variable('SAPT IND ENERGY') * 27211.4
                        edisp[i] = psi4.variable('SAPT DISP ENERGY') * 27211.4
                        eexch[i] = psi4.variable('SAPT EXCH ENERGY') * 27211.4
                        esapt[i] = psi4.variable('SAPT TOTAL ENERGY') * 27211.4
                        psi4.core.clean()
                    if metodo == 'sapt2+':
                        psi4.geometry(dimero)
                        psi4.energy(nivel)
                        eelst[i] = psi4.variable('SAPT ELST ENERGY') * 27211.4
                        eind[i] = psi4.variable('SAPT IND ENERGY') * 27211.4
                        edisp[i] = psi4.variable('SAPT DISP ENERGY') * 27211.4
                        eexch[i] = psi4.variable('SAPT EXCH ENERGY') * 27211.4
                        esapt[i] = psi4.variable('SAPT TOTAL ENERGY') * 27211.4
                        psi4.core.clean()
                    if metodo == 'sapt2+(3)':
                        psi4.geometry(dimero)
                        psi4.energy(nivel)
                        eelst[i] = psi4.variable('SAPT ELST ENERGY') * 27211.4
                        eind[i] = psi4.variable('SAPT IND ENERGY') * 27211.4
                        edisp[i] = psi4.variable('SAPT DISP ENERGY') * 27211.4
                        eexch[i] = psi4.variable('SAPT EXCH ENERGY') * 27211.4
                        esapt[i] = psi4.variable('SAPT TOTAL ENERGY') * 27211.4
                        psi4.core.clean()
                    if metodo == 'sapt2+3':
                        psi4.geometry(dimero)
                        psi4.energy(nivel)
                        eelst[i] = psi4.variable('SAPT ELST ENERGY') * 27211.4
                        eind[i] = psi4.variable('SAPT IND ENERGY') * 27211.4
                        edisp[i] = psi4.variable('SAPT DISP ENERGY') * 27211.4
                        eexch[i] = psi4.variable('SAPT EXCH ENERGY') * 27211.4
                        esapt[i] = psi4.variable('SAPT TOTAL ENERGY') * 27211.4
                        psi4.core.clean()



                sitio_inte = geo.replace('sapt.xyz', '').replace('amonia', '')
                nome_arq_out = metodo +  sitio_inte + base + '.dat'

                cria_arquivo(nome_arq_out, metodo, distancias, en_sem_cp, en_com_cp,
                             eelst, eind, edisp, eexch, esapt)
                move_arquivo(nome_arq_out, gas_nobre, metodo, base)
