from glob import glob
import numpy as np
import psi4
import re
import subprocess
import time


def cria_diretorio(gas, metodo='sem_metodo', base='sem_base'):
    '''
    Essa função recebe três strings que correspondem ao gás nobre de interesse
    o método e a base que iremos executar os cálculos e cria um novo diretorio
    com as respectivas palavras.
    '''
    if base == 'sem_base' and metodo == 'sem_metodo':
        subprocess.run(['mkdir', f'{gas}'])
    else:
        subprocess.run(['mkdir', f'{gas}_{metodo}_{base}'])


def cria_matriz(d):
    '''
    Essa função recebe as distâncias em angstrom que correspondem ao
    intervalo do cálculo e retorna uma matriz (array) com um conjunto de linhas
    e três colunas cujo os elementos são zeros.
    '''
    return np.zeros(((len(d)), 3))


def input_geo(geo, gas, d):
    '''
    Como parâmetros de entrada a função recebe a geometria do composto (geo)
    os símbolos dos gases nobres (gas) as distancias de cada interação (d) e
    um contador i e retorna a geometria do input.
    '''

    input = """
    0 1 """ + """\n""" + geo + """--
    0 1 """ + """\n""" + gas + """  """ + str(d) + """  0.000000   0.00000000
    units angstrom
    symmetry c1
    """
    return input


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

    if metodo == 'sherrill_gold_standard':
        with open(nome, 'w') as f:
            print(f'# Distancia [angstrom]    |   Energia [meV]', end='\n', file=f)
            print(f'# -----------------------------------------', end='\n', file=f)
            print(f'# Dist    Eint SGS-NOCP    Eint SGS-CP', end='\n', file=f)
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


def move_diretorio(gas, metodo, base):
    subprocess.run(['mv',  f'{gas}_{metodo}_{base}', f'{gas}/'])

# ---------------------------------------------------------------------
#  Iniciamos a execução deste script considerando as geometrias que
#  correspondem a molécula de amônia. Essas geometrias estão no
#  diretorio corrente e foram orientadas segundo o sítio de interação
#  de interesse (sítio que desejamos ver as interações).
#
#  Todos os métodos foram estabelecidos da lista métodos bem como as
#  bases e os gases que desejamos interagir em suas respectivas listas.
#
#  --------------------------------------------------------------------
#
#  Início Aqui! :)
#
#  --------------------------------------------------------------------

psi4.set_memory('15 GB')
psi4.set_num_threads(8)
print('\n')
psi4.core.set_output_file('output.dat', False)
psi4.set_options({'freeze_core': 'true'})
numpy_memory = 15

geometrias_amonia = glob('*_sapt.xyz')

metodos = ['ccsd', 'ccsd(t)', 'mp2', 'mp4', 'sapt0','sapt2', 'sapt2+',
           'sapt2+(3)', 'sapt2+3', 'sherrill_gold_standard']

bases = ['jun-cc-pvdz', 'jun-cc-pvdz', 'aug-cc-pvdz', 'aug-cc-pvtz',]

gases_nobres = ['He', 'Ne', 'Ar', 'Kr']

# int_ini define o ponto inicial da CEP.
int_ini = 3.2
# int_final define o ponto final da CEP
int_final = 8.6
# inc_ini define o incremento inicial da CEP e será modificado devido as
# mudanças na energia.
inc_ini = 0.2
# conta_min define um contador a partir do ponto mínimo da CEP.
conta_min = 0
# conta_pass_min conta a quantidade de passos com menor incremento próximo
# ao mínimo de energia na CEP
conta_pass_min = 14
# nova_dist1 define o intervalo de novas distâncias dentro de um intervalo
# devido ao fato da energia ir diminuindo.
nova_dist1 = 0

for gas_nobre in gases_nobres:
    cria_diretorio(gas_nobre)
    for metodo in metodos:
        for base in bases:
            if metodo != 'sherrill_gold_standard':
                cria_diretorio(gas_nobre, metodo, base)
            else:
                cria_diretorio(gas_nobre, metodo, base='')

            nivel = f'{metodo}/{base}'
            for geo in geometrias_amonia:
                with open(geo, 'r') as f:
                    str_geo = f.read()

                distancias = [int_ini]
                dist = int_ini

                en_sem_cp =  []
                en_com_cp = []

                eelst = []
                eind = []
                edisp = []
                eexch = []
                esapt = []

                while dist < int_final:
                    # inc_iniruindo a geometria do Dimero
                    dimero = input_geo(str_geo, gas_nobre, dist)

                    if metodo == 'sherrill_gold_standard':
                        n = len(distancias)
                        if n == 1:
                            psi4.geometry(dimero)
                            psi4.energy(metodo, bsse_type=['nocp', 'cp',])
                            en_sem_cp.append(psi4.variable('NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY') * 27211.4)
                            en_com_cp.append(psi4.variable('CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY') * 27211.4)
                            psi4.core.clean()
                            dist = round(int_ini + inc_ini, 2)
                            distancias.append(dist)
                        else:
                            psi4.geometry(dimero)
                            psi4.energy(metodo, bsse_type=['nocp', 'cp',])
                            en_sem_cp.append(psi4.variable('NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY') * 27211.4)
                            en_com_cp.append(psi4.variable('CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY') * 27211.4)
                            psi4.core.clean()
                            if en_com_cp[-1] - en_com_cp[-2] < 0 and conta_min <= 2:
                                nova_dist1 = round(dist+0.1, 2)
                            elif en_com_cp[-1] - en_com_cp[-2] < 0 and conta_min > 2 and conta_min < conta_pass_min:
                                nova_dist1 = round(dist+0.05, 2)
                            else:
                                nova_dist1 = round(dist + inc_ini, 2)
                            distancias.append(nova_dist1)
                            dist = nova_dist1
                            conta_min += 1

                    if metodo == 'ccsd(t)':
                        n = len(distancias)
                        if n == 1:
                            psi4.geometry(dimero)
                            psi4.energy(nivel, bsse_type=['nocp', 'cp',])
                            en_sem_cp.append(psi4.variable('NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY') * 27211.4)
                            en_com_cp.append(psi4.variable('CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY') * 27211.4)
                            psi4.core.clean()
                            dist = round(int_ini + inc_ini, 2)
                            distancias.append(dist)
                        else:
                            psi4.geometry(dimero)
                            psi4.energy(nivel, bsse_type=['nocp', 'cp',])
                            en_sem_cp.append(psi4.variable('NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY') * 27211.4)
                            en_com_cp.append(psi4.variable('CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY') * 27211.4)
                            psi4.core.clean()
                            if en_com_cp[-1] - en_com_cp[-2] < 0 and conta_min <= 2:
                                nova_dist1 = round(dist+0.1, 2)
                            elif en_com_cp[-1] - en_com_cp[-2] < 0 and conta_min > 2 and conta_min < conta_pass_min:
                                nova_dist1 = round(dist+0.05, 2)
                            else:
                                nova_dist1 = round(dist + inc_ini, 2)
                            distancias.append(nova_dist1)
                            dist = nova_dist1
                            conta_min += 1

                    if metodo == 'ccsd':
                        n = len(distancias)
                        if n == 1:
                            psi4.geometry(dimero)
                            psi4.energy(nivel, bsse_type=['nocp', 'cp',])
                            en_sem_cp.append(psi4.variable('NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY') * 27211.4)
                            en_com_cp.append(psi4.variable('CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY') * 27211.4)
                            psi4.core.clean()
                            dist = round(int_ini + inc_ini, 2)
                            distancias.append(dist)
                        else:
                            psi4.geometry(dimero)
                            psi4.energy(nivel, bsse_type=['nocp', 'cp',])
                            en_sem_cp.append(psi4.variable('NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY') * 27211.4)
                            en_com_cp.append(psi4.variable('CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY') * 27211.4)
                            psi4.core.clean()
                            if en_com_cp[-1] - en_com_cp[-2] < 0 and conta_min <= 2:
                                nova_dist1 = round(dist+0.1, 2)
                            elif en_com_cp[-1] - en_com_cp[-2] < 0 and conta_min > 2 and conta_min < conta_pass_min:
                                nova_dist1 = round(dist+0.05, 2)
                            else:
                                nova_dist1 = round(dist + inc_ini, 2)
                            distancias.append(nova_dist1)
                            dist = nova_dist1
                            conta_min += 1

                    if metodo == 'mp2':
                        n = len(distancias)
                        if n == 1:
                            psi4.geometry(dimero)
                            psi4.energy(nivel, bsse_type=['nocp', 'cp',])
                            en_sem_cp.append(psi4.variable('NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY') * 27211.4)
                            en_com_cp.append(psi4.variable('CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY') * 27211.4)
                            psi4.core.clean()
                            dist = round(int_ini + inc_ini, 2)
                            distancias.append(dist)
                        else:
                            psi4.geometry(dimero)
                            psi4.energy(nivel, bsse_type=['nocp', 'cp',])
                            en_sem_cp.append(psi4.variable('NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY') * 27211.4)
                            en_com_cp.append(psi4.variable('CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY') * 27211.4)
                            psi4.core.clean()
                            if en_com_cp[-1] - en_com_cp[-2] < 0 and conta_min <= 2:
                                nova_dist1 = round(dist+0.1, 2)
                            elif en_com_cp[-1] - en_com_cp[-2] < 0 and conta_min > 2 and conta_min < conta_pass_min:
                                nova_dist1 = round(dist+0.05, 2)
                            else:
                                nova_dist1 = round(dist + inc_ini, 2)
                            distancias.append(nova_dist1)
                            dist = nova_dist1
                            conta_min += 1

                    if metodo == 'mp4':
                        n = len(distancias)
                        if n == 1:
                            psi4.geometry(dimero)
                            psi4.energy(nivel, bsse_type=['nocp', 'cp',])
                            en_sem_cp.append(psi4.variable('NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY') * 27211.4)
                            en_com_cp.append(psi4.variable('CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY') * 27211.4)
                            psi4.core.clean()
                            dist = round(int_ini + inc_ini, 2)
                            distancias.append(dist)
                        else:
                            psi4.geometry(dimero)
                            psi4.energy(nivel, bsse_type=['nocp', 'cp',])
                            en_sem_cp.append(psi4.variable('NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY') * 27211.4)
                            en_com_cp.append(psi4.variable('CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY') * 27211.4)
                            psi4.core.clean()
                            if en_com_cp[-1] - en_com_cp[-2] < 0 and conta_min <= 2:
                                nova_dist1 = round(dist+0.1, 2)
                            elif en_com_cp[-1] - en_com_cp[-2] < 0 and conta_min > 2 and conta_min < conta_pass_min:
                                nova_dist1 = round(dist+0.05, 2)
                            else:
                                nova_dist1 = round(dist + inc_ini, 2)
                            distancias.append(nova_dist1)
                            dist = nova_dist1
                            conta_min += 1

                    if metodo == 'sapt0':
                        n = len(distancias)
                        if n == 1:
                            psi4.geometry(dimero)
                            psi4.energy(nivel)
                            eelst.append(psi4.variable('SAPT ELST ENERGY') * 27211.4)
                            eind.append(psi4.variable('SAPT IND ENERGY') * 27211.4)
                            edisp.append(psi4.variable('SAPT DISP ENERGY') * 27211.4)
                            eexch.append(psi4.variable('SAPT EXCH ENERGY') * 27211.4)
                            esapt.append(psi4.variable('SAPT TOTAL ENERGY') * 27211.4)
                            psi4.core.clean()
                            dist = round(int_ini + inc_ini, 2)
                            distancias.append(dist)
                        else:
                            psi4.geometry(dimero)
                            psi4.energy(nivel)
                            eelst.append(psi4.variable('SAPT ELST ENERGY') * 27211.4)
                            eind.append(psi4.variable('SAPT IND ENERGY') * 27211.4)
                            edisp.append(psi4.variable('SAPT DISP ENERGY') * 27211.4)
                            eexch.append(psi4.variable('SAPT EXCH ENERGY') * 27211.4)
                            esapt.append(psi4.variable('SAPT TOTAL ENERGY') * 27211.4)
                            psi4.core.clean()
                            if esapt[-1] - esapt[-2] < 0 and conta_min <= 2:
                                nova_dist1 = round(dist+0.1, 2)
                            elif esapt[-1] - esapt[-2] < 0 and conta_min > 2 and conta_min < conta_pass_min:
                                nova_dist1 = round(dist+0.05, 2)
                            else:
                                nova_dist1 = round(dist + inc_ini, 2)
                            distancias.append(nova_dist1)
                            dist = nova_dist1
                            conta_min += 1

                    if metodo == 'sapt2':
                            n = len(distancias)
                            if n == 1:
                                psi4.geometry(dimero)
                                psi4.energy(nivel)
                                eelst.append(psi4.variable('SAPT ELST ENERGY') * 27211.4)
                                eind.append(psi4.variable('SAPT IND ENERGY') * 27211.4)
                                edisp.append(psi4.variable('SAPT DISP ENERGY') * 27211.4)
                                eexch.append(psi4.variable('SAPT EXCH ENERGY') * 27211.4)
                                esapt.append(psi4.variable('SAPT TOTAL ENERGY') * 27211.4)
                                psi4.core.clean()
                                dist = round(int_ini + inc_ini, 2)
                                distancias.append(dist)
                            else:
                                psi4.geometry(dimero)
                                psi4.energy(nivel)
                                eelst.append(psi4.variable('SAPT ELST ENERGY') * 27211.4)
                                eind.append(psi4.variable('SAPT IND ENERGY') * 27211.4)
                                edisp.append(psi4.variable('SAPT DISP ENERGY') * 27211.4)
                                eexch.append(psi4.variable('SAPT EXCH ENERGY') * 27211.4)
                                esapt.append(psi4.variable('SAPT TOTAL ENERGY') * 27211.4)
                                psi4.core.clean()
                                if esapt[-1] - esapt[-2] < 0 and conta_min <= 2:
                                    nova_dist1 = round(dist+0.1, 2)
                                elif esapt[-1] - esapt[-2] < 0 and conta_min > 2 and conta_min < conta_pass_min:
                                    nova_dist1 = round(dist+0.05, 2)
                                else:
                                    nova_dist1 = round(dist + inc_ini, 2)
                                distancias.append(nova_dist1)
                                dist = nova_dist1
                                conta_min += 1

                    if metodo == 'sapt2+':
                            n = len(distancias)
                            if n == 1:
                                psi4.geometry(dimero)
                                psi4.energy(nivel)
                                eelst.append(psi4.variable('SAPT ELST ENERGY') * 27211.4)
                                eind.append(psi4.variable('SAPT IND ENERGY') * 27211.4)
                                edisp.append(psi4.variable('SAPT DISP ENERGY') * 27211.4)
                                eexch.append(psi4.variable('SAPT EXCH ENERGY') * 27211.4)
                                esapt.append(psi4.variable('SAPT TOTAL ENERGY') * 27211.4)
                                psi4.core.clean()
                                dist = round(int_ini + inc_ini, 2)
                                distancias.append(dist)
                            else:
                                psi4.geometry(dimero)
                                psi4.energy(nivel)
                                eelst.append(psi4.variable('SAPT ELST ENERGY') * 27211.4)
                                eind.append(psi4.variable('SAPT IND ENERGY') * 27211.4)
                                edisp.append(psi4.variable('SAPT DISP ENERGY') * 27211.4)
                                eexch.append(psi4.variable('SAPT EXCH ENERGY') * 27211.4)
                                esapt.append(psi4.variable('SAPT TOTAL ENERGY') * 27211.4)
                                psi4.core.clean()
                                if esapt[-1] - esapt[-2] < 0 and conta_min <= 2:
                                    nova_dist1 = round(dist+0.1, 2)
                                elif esapt[-1] - esapt[-2] < 0 and conta_min > 2 and conta_min < conta_pass_min:
                                    nova_dist1 = round(dist+0.05, 2)
                                else:
                                    nova_dist1 = round(dist + inc_ini, 2)
                                distancias.append(nova_dist1)
                                dist = nova_dist1
                                conta_min += 1

                    if metodo == 'sapt2+(3)':
                            n = len(distancias)
                            if n == 1:
                                psi4.geometry(dimero)
                                psi4.energy(nivel)
                                eelst.append(psi4.variable('SAPT ELST ENERGY') * 27211.4)
                                eind.append(psi4.variable('SAPT IND ENERGY') * 27211.4)
                                edisp.append(psi4.variable('SAPT DISP ENERGY') * 27211.4)
                                eexch.append(psi4.variable('SAPT EXCH ENERGY') * 27211.4)
                                esapt.append(psi4.variable('SAPT TOTAL ENERGY') * 27211.4)
                                psi4.core.clean()
                                dist = round(int_ini + inc_ini, 2)
                                distancias.append(dist)
                            else:
                                psi4.geometry(dimero)
                                psi4.energy(nivel)
                                eelst.append(psi4.variable('SAPT ELST ENERGY') * 27211.4)
                                eind.append(psi4.variable('SAPT IND ENERGY') * 27211.4)
                                edisp.append(psi4.variable('SAPT DISP ENERGY') * 27211.4)
                                eexch.append(psi4.variable('SAPT EXCH ENERGY') * 27211.4)
                                esapt.append(psi4.variable('SAPT TOTAL ENERGY') * 27211.4)
                                psi4.core.clean()
                                if esapt[-1] - esapt[-2] < 0 and conta_min <= 2:
                                    nova_dist1 = round(dist+0.1, 2)
                                elif esapt[-1] - esapt[-2] < 0 and conta_min > 2 and conta_min < conta_pass_min:
                                    nova_dist1 = round(dist+0.05, 2)
                                else:
                                    nova_dist1 = round(dist + inc_ini, 2)
                                distancias.append(nova_dist1)
                                dist = nova_dist1
                                conta_min += 1

                    if metodo == 'sapt2+3':
                        n = len(distancias)
                        if n == 1:
                            psi4.geometry(dimero)
                            psi4.energy(nivel)
                            eelst.append(psi4.variable('SAPT ELST ENERGY') * 27211.4)
                            eind.append(psi4.variable('SAPT IND ENERGY') * 27211.4)
                            edisp.append(psi4.variable('SAPT DISP ENERGY') * 27211.4)
                            eexch.append(psi4.variable('SAPT EXCH ENERGY') * 27211.4)
                            esapt.append(psi4.variable('SAPT TOTAL ENERGY') * 27211.4)
                            psi4.core.clean()
                            dist = round(int_ini + inc_ini, 2)
                            distancias.append(dist)
                        else:
                            psi4.geometry(dimero)
                            psi4.energy(nivel)
                            eelst.append(psi4.variable('SAPT ELST ENERGY') * 27211.4)
                            eind.append(psi4.variable('SAPT IND ENERGY') * 27211.4)
                            edisp.append(psi4.variable('SAPT DISP ENERGY') * 27211.4)
                            eexch.append(psi4.variable('SAPT EXCH ENERGY') * 27211.4)
                            esapt.append(psi4.variable('SAPT TOTAL ENERGY') * 27211.4)
                            psi4.core.clean()
                            if esapt[-1] - esapt[-2] < 0 and conta_min <= 2:
                                nova_dist1 = round(dist+0.1, 2)
                            elif esapt[-1] - esapt[-2] < 0 and conta_min > 2 and conta_min < conta_pass_min:
                                nova_dist1 = round(dist+0.05, 2)
                            else:
                                nova_dist1 = round(dist + inc_ini, 2)
                            distancias.append(nova_dist1)
                            dist = nova_dist1
                            conta_min += 1

                if metodo != 'sherrill_gold_standard':
                    sitio_inte = geo.replace('sapt.xyz', '').replace('amonia', '')
                    nome_arq_out = metodo +  sitio_inte + base + '.dat'
                    cria_arquivo(nome_arq_out, metodo, distancias, en_sem_cp, en_com_cp,
                                 eelst, eind, edisp, eexch, esapt)
                    move_arquivo(nome_arq_out, gas_nobre, metodo, base)
                else:
                    sitio_inte = geo.replace('sapt.xyz', '').replace('amonia', '')
                    nome_arq_out = metodo +  sitio_inte + '.dat'
                    cria_arquivo(nome_arq_out, metodo, distancias, en_sem_cp, en_com_cp,
                                 eelst, eind, edisp, eexch, esapt)
                    move_arquivo(nome_arq_out, gas_nobre, metodo, base='_')
            move_diretorio(gas_nobre, metodo, base)
