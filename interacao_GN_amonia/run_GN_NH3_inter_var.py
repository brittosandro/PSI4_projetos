from glob import glob
import numpy as np
import psi4
import re
import subprocess
import time


def calcula_energia(metodo, base, dimero, fator_conversao=1):
    '''
    Essa função calcula a energia de interação de certo sistema, portanto
    recebe como parâmetros o método do cálculo a base e o sistemas que neste
    caso chamamos de dímero. Caso queira mudar a unidade de medida padrão
    da energia que é hartree, deverá atribuir um valor para o fator de conversão.
    Se não for passado nenhum valor de fator de conversão
    esse será igual a 1, pois 1 é elemento neutro da multiplicação.

    A função retornará uma tupla com os valores de energia. Dependendo do
    método de cálculo a tupla irá retornar dois ou mais valores.

    Para métodos como ccsd, ccsd(t), mp2, mp4 ou sherrill_gold_standard a tupla
    é formada pela energia sem cálculo de counterpoise e com cálculo de
    conterpoise.

    Exemplo1: energia[ccsd/aug-cc-pvdz] = (energia_sem_counterpoisi,
                                           energia_com_conterpoise)

    O Exemplo1 mostra a energia sendo calculada com método ccsd na base aug-cc-
    pvdz e a tupla que é retornada pela função.

    No caso do método ser sapt a função retorna uma tupla com  as componentes
    energia eletrostática, energia de indução, energia de dispersão, energia
    de troca e energia sapt total.

    Exemplo2: energia[sapt0/aug-cc-pvdz] = (energia_eletrostática,
                                            energia_indução,
                                            energia_dispersão,
                                            energia_troca,
                                            energia_total)

    O Exemplo2 mostra a energia sendo calculada com método sapt0 na base aug-cc-
    pvdz e a tupla que a função retorna.
    '''

    metodo_nao_e_sapt = ('ccsd', 'ccsd(t)', 'mp2', 'mp4')

    if metodo == 'sherrill_gold_standard':
        psi4.geometry(dimero)
        psi4.energy(f'{metodo}', bsse_type=['nocp', 'cp',])

        e_s_cp = psi4.variable(
                'NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY') * fator_conversao

        e_c_cp = psi4.variable(
                'CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY') * fator_conversao
        psi4.core.clean()
        energia = e_s_cp, e_c_cp
    elif metodo in metodo_nao_e_sapt:
        psi4.geometry(dimero)
        psi4.energy(f'{metodo}/{base}', bsse_type=['nocp', 'cp',])

        e_s_cp = psi4.variable(
                'NOCP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY') * fator_conversao

        e_c_cp = psi4.variable(
                'CP-CORRECTED INTERACTION ENERGY THROUGH 2-BODY') * fator_conversao
        psi4.core.clean()
        energia = e_s_cp, e_c_cp
    else:
        psi4.geometry(dimero)
        psi4.energy(f'{metodo}/{base}')
        eels = psi4.variable('SAPT ELST ENERGY') * fator_conversao
        eind = psi4.variable('SAPT IND ENERGY') * fator_conversao
        edis = psi4.variable('SAPT DISP ENERGY') * fator_conversao
        exch = psi4.variable('SAPT EXCH ENERGY') * fator_conversao
        esap = psi4.variable('SAPT TOTAL ENERGY') * fator_conversao
        energia = eels, eind, edis, exch, esap

    return energia

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
    os símbolos dos gases nobres (gas) as distancias de cada interação (d).
    '''

    input = """
    0 1 """ + """\n""" + geo + """--
    0 1 """ + """\n""" + gas + """  """ + str(d) + """  0.000000   0.00000000
    units angstrom
    symmetry c1
    """
    return input


def cria_arquivo(nome, metodo, dist, eint1, eint2, eint3, eint4, eint5, eint6, eint7,
                 fator_conversao):
    unidade_energia = fator_conversao.split('2')[1]
    if metodo == 'ccsd(t)':
        with open(nome, 'w') as f:
            print(f'# Distancia [angstrom]    |   Energia [{unidade_energia}]', end='\n', file=f)
            print(f'# -----------------------------------------', end='\n', file=f)
            print(f'# Dist   Eint CCSD(T)-NOCP  Eint CCSD(T)-CP', end='\n', file=f)
            print(f'# -----------------------------------------', end='\n', file=f)
            for d, e1, e2 in zip(dist, np.around(eint1, 6), np.around(eint2, 6)):
                print(f'{d:5.2f} {e1:16.9f} {e2:17.9f}', end='\n', file=f)

    if metodo == 'ccsd':
        with open(nome, 'w') as f:
            print(f'# Distancia [angstrom]    |   Energia [{unidade_energia}]', end='\n', file=f)
            print(f'# -----------------------------------------', end='\n', file=f)
            print(f'# Dist    Eint CCSD-NOCP    Eint CCSD-CP', end='\n', file=f)
            print(f'# -----------------------------------------', end='\n', file=f)
            for d, e1, e2 in zip(dist, np.around(eint1, 6), np.around(eint2, 6)):
                print(f'{d:6.2f} {e1:16.9f} {e2:17.9f}', end='\n', file=f)

    if metodo == 'mp2':
        with open(nome, 'w') as f:
            print(f'# Distancia [angstrom]    |   Energia [{unidade_energia}]', end='\n', file=f)
            print(f'# -----------------------------------------', end='\n', file=f)
            print(f'# Dist    Eint MP2-NOCP    Eint MP2-CP', end='\n', file=f)
            print(f'# -----------------------------------------', end='\n', file=f)
            for d, e1, e2 in zip(dist, np.around(eint1, 6), np.around(eint2, 6)):
                print(f'{d:6.2f} {e1:16.9f} {e2:17.9f}', end='\n', file=f)

    if metodo == 'mp4':
        with open(nome, 'w') as f:
            print(f'# Distancia [angstrom]    |   Energia [{unidade_energia}]', end='\n', file=f)
            print(f'# -----------------------------------------', end='\n', file=f)
            print(f'# Dist    Eint MP4-NOCP    Eint MP4-CP', end='\n', file=f)
            print(f'# -----------------------------------------', end='\n', file=f)
            for d, e1, e2 in zip(dist, np.around(eint1, 6), np.around(eint2, 6)):
                print(f'{d:6.2f} {e1:16.9f} {e2:15.9f}', end='\n', file=f)

    if metodo == 'sherrill_gold_standard':
        with open(nome, 'w') as f:
            print(f'# Distancia [angstrom]    |   Energia [{unidade_energia}]', end='\n', file=f)
            print(f'# -----------------------------------------', end='\n', file=f)
            print(f'# Dist    Eint SGS-NOCP    Eint SGS-CP', end='\n', file=f)
            print(f'# -----------------------------------------', end='\n', file=f)
            for d, e1, e2 in zip(dist, np.around(eint1, 6), np.around(eint2, 6)):
                print(f'{d:6.2f} {e1:16.9f} {e2:15.9f}', end='\n', file=f)

    metodos_perturbativos = ('sapt0', 'sapt2', 'sapt2+', 'sapt2+(3)', 'sapt2+3')
    if metodo in metodos_perturbativos:
        with open(nome_arq_out, 'w') as f:
            print(f'#                        Distancia [angstrom]    |   Energia [{unidade_energia}]', end='\n', file=f)
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

metodos = ['ccsd', 'ccsd(t)', 'mp2', 'mp4',]

#metodos = ['sherrill_gold_standard']

#metodos = ['sapt0', 'sapt2', 'sapt2+', 'sapt2+(3)']

bases = ['jun-cc-pvdz',]

#gases_nobres = ['He', 'Ne', 'Ar', 'Kr']
gases_nobres = ['He',]

# int_ini define o ponto inicial da CEP.
int_ini = 3.2

# int_final define o ponto final da CEP
int_final = 4.5

# inc_ini define o incremento inicial da CEP e será modificado devido as
# mudanças na energia.
inc_ini = 0.2

# conta_min_energia é um contador para saber quantas vezes a energia é
# minimizada a medida que a energia é calcula.
conta_min_energia = 0

# nova_dist1 define o intervalo de novas distâncias dentro de um cert intervalo
# devido ao fato da energia diminuir.
nova_dist1 = 0

# Caso queiramos apresentar um fator de conversão de energia que seja
# conveniente com nosso estudo poderemos listar esse fator pelas variáveis
# abaixo.
# No caso do estudo da interação entre gases nobres e amônia é muito conveniente
# utilizarmos mili eletron volts (meV), portanto devemos transformar de
# hartree para meV. Assim a variável é hartree2meV, como o fator de conversão
# é 27211.399
hartree2meV = 27211.399
# fator_conv é a variavel que irá receber o fator de conversão do seu interesse
fator_conv = hartree2meV
# tipo_conversao_energia revebe a conversão que buscamos executar, neste caso
# queremos converter de hartree para meV. É importante para apresentar os dados
# na função que cria os arquivos
tipo_conversao_energia = 'hartree2meV'

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

                while dist <= int_final:
                    # inc_iniruindo a geometria do Dimero
                    dimero = input_geo(str_geo, gas_nobre, dist)

                    if metodo == 'sherrill_gold_standard':
                        n = len(distancias)
                        if n == 1:
                            en1 = calcula_energia(metodo, base, dimero, fator_conv)[0]
                            en2 = calcula_energia(metodo, base, dimero, fator_conv)[1]
                            en_sem_cp.append(en1)
                            en_com_cp.append(en2)
                            dist = round(int_ini + inc_ini, 2)
                            distancias.append(dist)
                        else:
                            en1 = calcula_energia(metodo, base, dimero, fator_conv)[0]
                            en2 = calcula_energia(metodo, base, dimero, fator_conv)[1]
                            en_sem_cp.append(en1)
                            en_com_cp.append(en2)
                            if en_com_cp[-1] - en_com_cp[-2] < 0:
                                nova_dist1 = round(dist+0.1, 2)
                                conta_min_energia += 1
                            elif en_com_cp[-1] - en_com_cp[-2] < 0 and conta_min_energia >= 2:
                                nova_dist1 = round(dist+0.05, 2)
                                conta_min_energia += 1
                            else:
                                nova_dist1 = round(dist + inc_ini, 2)
                            distancias.append(nova_dist1)
                            dist = nova_dist1

                    if metodo == 'ccsd(t)':
                        n = len(distancias)
                        if n == 1:
                            en1 = calcula_energia(metodo, base, dimero, fator_conv)[0]
                            en2 = calcula_energia(metodo, base, dimero, fator_conv)[1]
                            en_sem_cp.append(en1)
                            en_com_cp.append(en2)
                            dist = round(int_ini + inc_ini, 2)
                            distancias.append(dist)
                        else:
                            en1 = calcula_energia(metodo, base, dimero, fator_conv)[0]
                            en2 = calcula_energia(metodo, base, dimero, fator_conv)[1]
                            en_sem_cp.append(en1)
                            en_com_cp.append(en2)
                            if en_com_cp[-1] - en_com_cp[-2] < 0:
                                nova_dist1 = round(dist+0.1, 2)
                                conta_min_energia += 1
                            elif en_com_cp[-1] - en_com_cp[-2] < 0 and conta_min_energia >= 2:
                                nova_dist1 = round(dist+0.05, 2)
                                conta_min_energia += 1
                            else:
                                nova_dist1 = round(dist + inc_ini, 2)
                            distancias.append(nova_dist1)
                            dist = nova_dist1

                    if metodo == 'ccsd':
                        n = len(distancias)
                        if n == 1:
                            en1 = calcula_energia(metodo, base, dimero, fator_conv)[0]
                            en2 = calcula_energia(metodo, base, dimero, fator_conv)[1]
                            en_sem_cp.append(en1)
                            en_com_cp.append(en2)
                            dist = round(int_ini + inc_ini, 2)
                            distancias.append(dist)
                        else:
                            en1 = calcula_energia(metodo, base, dimero, fator_conv)[0]
                            en2 = calcula_energia(metodo, base, dimero, fator_conv)[1]
                            en_sem_cp.append(en1)
                            en_com_cp.append(en2)
                            if en_com_cp[-1] - en_com_cp[-2] < 0:
                                nova_dist1 = round(dist+0.1, 2)
                                conta_min_energia += 1
                            elif en_com_cp[-1] - en_com_cp[-2] < 0 and conta_min_energia >= 2:
                                nova_dist1 = round(dist+0.05, 2)
                                conta_min_energia += 1
                            else:
                                nova_dist1 = round(dist + inc_ini, 2)
                            distancias.append(nova_dist1)
                            dist = nova_dist1

                    if metodo == 'mp2':
                        n = len(distancias)
                        if n == 1:
                            en1 = calcula_energia(metodo, base, dimero, fator_conv)[0]
                            en2 = calcula_energia(metodo, base, dimero, fator_conv)[1]
                            en_sem_cp.append(en1)
                            en_com_cp.append(en2)
                            dist = round(int_ini + inc_ini, 2)
                            distancias.append(dist)
                        else:
                            en1 = calcula_energia(metodo, base, dimero, fator_conv)[0]
                            en2 = calcula_energia(metodo, base, dimero, fator_conv)[1]
                            en_sem_cp.append(en1)
                            en_com_cp.append(en2)
                            if en_com_cp[-1] - en_com_cp[-2] < 0:
                                nova_dist1 = round(dist+0.1, 2)
                                conta_min_energia += 1
                            elif en_com_cp[-1] - en_com_cp[-2] < 0 and conta_min_energia >= 2:
                                nova_dist1 = round(dist+0.05, 2)
                                conta_min_energia += 1
                            else:
                                nova_dist1 = round(dist + inc_ini, 2)
                            distancias.append(nova_dist1)
                            dist = nova_dist1

                    if metodo == 'mp4':
                        n = len(distancias)
                        if n == 1:
                            en1 = calcula_energia(metodo, base, dimero, fator_conv)[0]
                            en2 = calcula_energia(metodo, base, dimero, fator_conv)[1]
                            en_sem_cp.append(en1)
                            en_com_cp.append(en2)
                            dist = round(int_ini + inc_ini, 2)
                            distancias.append(dist)
                        else:
                            en1 = calcula_energia(metodo, base, dimero, fator_conv)[0]
                            en2 = calcula_energia(metodo, base, dimero, fator_conv)[1]
                            en_sem_cp.append(en1)
                            en_com_cp.append(en2)
                            if en_com_cp[-1] - en_com_cp[-2] < 0:
                                nova_dist1 = round(dist+0.1, 2)
                                conta_min_energia += 1
                            elif en_com_cp[-1] - en_com_cp[-2] < 0 and conta_min_energia >= 2:
                                nova_dist1 = round(dist+0.05, 2)
                                conta_min_energia += 1
                            else:
                                nova_dist1 = round(dist + inc_ini, 2)
                            distancias.append(nova_dist1)
                            dist = nova_dist1

                    if metodo == 'sapt0':
                        n = len(distancias)
                        if n == 1:
                            eel = calcula_energia(metodo, base, dimero, fator_conv)[0]
                            eid = calcula_energia(metodo, base, dimero, fator_conv)[1]
                            edi = calcula_energia(metodo, base, dimero, fator_conv)[2]
                            exh = calcula_energia(metodo, base, dimero, fator_conv)[3]
                            esap = calcula_energia(metodo, base, dimero, fator_conv)[4]

                            eelst.append(eel)
                            eind.append(eid)
                            edisp.append(edi)
                            eexch.append(exh)
                            esapt.append(esap)

                            dist = round(int_ini + inc_ini, 2)
                            distancias.append(dist)
                        else:
                            eel = calcula_energia(metodo, base, dimero, fator_conv)[0]
                            eid = calcula_energia(metodo, base, dimero, fator_conv)[1]
                            edi = calcula_energia(metodo, base, dimero, fator_conv)[2]
                            exh = calcula_energia(metodo, base, dimero, fator_conv)[3]
                            esap = calcula_energia(metodo, base, dimero, fator_conv)[4]

                            eelst.append(eel)
                            eind.append(eid)
                            edisp.append(edi)
                            eexch.append(exh)
                            esapt.append(esap)

                            if esapt[-1] - esapt[-2] < 0:
                                nova_dist1 = round(dist+0.1, 2)
                                conta_min_energia += 1
                            elif esapt[-1] - esapt[-2] < 0 and conta_min_energia >= 2:
                                nova_dist1 = round(dist+0.05, 2)
                                conta_min_energia += 1
                            else:
                                nova_dist1 = round(dist + inc_ini, 2)
                            distancias.append(nova_dist1)
                            dist = nova_dist1

                    if metodo == 'sapt2':
                            n = len(distancias)
                            if n == 1:
                                eel = calcula_energia(metodo, base, dimero, fator_conv)[0]
                                eid = calcula_energia(metodo, base, dimero, fator_conv)[1]
                                edi = calcula_energia(metodo, base, dimero, fator_conv)[2]
                                exh = calcula_energia(metodo, base, dimero, fator_conv)[3]
                                esap = calcula_energia(metodo, base, dimero, fator_conv)[4]

                                eelst.append(eel)
                                eind.append(eid)
                                edisp.append(edi)
                                eexch.append(exh)
                                esapt.append(esap)

                                dist = round(int_ini + inc_ini, 2)
                                distancias.append(dist)
                            else:
                                eel = calcula_energia(metodo, base, dimero, fator_conv)[0]
                                eid = calcula_energia(metodo, base, dimero, fator_conv)[1]
                                edi = calcula_energia(metodo, base, dimero, fator_conv)[2]
                                exh = calcula_energia(metodo, base, dimero, fator_conv)[3]
                                esap = calcula_energia(metodo, base, dimero, fator_conv)[4]

                                eelst.append(eel)
                                eind.append(eid)
                                edisp.append(edi)
                                eexch.append(exh)
                                esapt.append(esap)

                                if esapt[-1] - esapt[-2] < 0:
                                    nova_dist1 = round(dist+0.1, 2)
                                    conta_min_energia += 1
                                elif esapt[-1] - esapt[-2] < 0 and conta_min_energia >= 2:
                                    nova_dist1 = round(dist+0.05, 2)
                                    conta_min_energia += 1
                                else:
                                    nova_dist1 = round(dist + inc_ini, 2)
                                distancias.append(nova_dist1)
                                dist = nova_dist1

                    if metodo == 'sapt2+':
                            n = len(distancias)
                            if n == 1:
                                eel = calcula_energia(metodo, base, dimero, fator_conv)[0]
                                eid = calcula_energia(metodo, base, dimero, fator_conv)[1]
                                edi = calcula_energia(metodo, base, dimero, fator_conv)[2]
                                exh = calcula_energia(metodo, base, dimero, fator_conv)[3]
                                esap = calcula_energia(metodo, base, dimero, fator_conv)[4]

                                eelst.append(eel)
                                eind.append(eid)
                                edisp.append(edi)
                                eexch.append(exh)
                                esapt.append(esap)

                                dist = round(int_ini + inc_ini, 2)
                                distancias.append(dist)
                            else:
                                eel = calcula_energia(metodo, base, dimero, fator_conv)[0]
                                eid = calcula_energia(metodo, base, dimero, fator_conv)[1]
                                edi = calcula_energia(metodo, base, dimero, fator_conv)[2]
                                exh = calcula_energia(metodo, base, dimero, fator_conv)[3]
                                esap = calcula_energia(metodo, base, dimero, fator_conv)[4]

                                eelst.append(eel)
                                eind.append(eid)
                                edisp.append(edi)
                                eexch.append(exh)
                                esapt.append(esap)

                                if esapt[-1] - esapt[-2] < 0:
                                    nova_dist1 = round(dist+0.1, 2)
                                    conta_min_energia += 1
                                elif esapt[-1] - esapt[-2] < 0 and conta_min_energia >= 2:
                                    nova_dist1 = round(dist+0.05, 2)
                                    conta_min_energia += 1
                                else:
                                    nova_dist1 = round(dist + inc_ini, 2)
                                distancias.append(nova_dist1)
                                dist = nova_dist1

                    if metodo == 'sapt2+(3)':
                            n = len(distancias)
                            if n == 1:
                                eel = calcula_energia(metodo, base, dimero, fator_conv)[0]
                                eid = calcula_energia(metodo, base, dimero, fator_conv)[1]
                                edi = calcula_energia(metodo, base, dimero, fator_conv)[2]
                                exh = calcula_energia(metodo, base, dimero, fator_conv)[3]
                                esap = calcula_energia(metodo, base, dimero, fator_conv)[4]

                                eelst.append(eel)
                                eind.append(eid)
                                edisp.append(edi)
                                eexch.append(exh)
                                esapt.append(esap)

                                dist = round(int_ini + inc_ini, 2)
                                distancias.append(dist)
                            else:
                                eel = calcula_energia(metodo, base, dimero, fator_conv)[0]
                                eid = calcula_energia(metodo, base, dimero, fator_conv)[1]
                                edi = calcula_energia(metodo, base, dimero, fator_conv)[2]
                                exh = calcula_energia(metodo, base, dimero, fator_conv)[3]
                                esap = calcula_energia(metodo, base, dimero, fator_conv)[4]

                                eelst.append(eel)
                                eind.append(eid)
                                edisp.append(edi)
                                eexch.append(exh)
                                esapt.append(esap)

                                if esapt[-1] - esapt[-2] < 0:
                                    nova_dist1 = round(dist+0.1, 2)
                                    conta_min_energia += 1
                                elif esapt[-1] - esapt[-2] < 0 and conta_min_energia >= 2:
                                    nova_dist1 = round(dist+0.05, 2)
                                    conta_min_energia += 1
                                else:
                                    nova_dist1 = round(dist + inc_ini, 2)
                                distancias.append(nova_dist1)
                                dist = nova_dist1

                    if metodo == 'sapt2+3':
                        n = len(distancias)
                        if n == 1:
                            eel = calcula_energia(metodo, base, dimero, fator_conv)[0]
                            eid = calcula_energia(metodo, base, dimero, fator_conv)[1]
                            edi = calcula_energia(metodo, base, dimero, fator_conv)[2]
                            exh = calcula_energia(metodo, base, dimero, fator_conv)[3]
                            esap = calcula_energia(metodo, base, dimero, fator_conv)[4]

                            eelst.append(eel)
                            eind.append(eid)
                            edisp.append(edi)
                            eexch.append(exh)
                            esapt.append(esap)

                            dist = round(int_ini + inc_ini, 2)
                            distancias.append(dist)
                        else:
                            eel = calcula_energia(metodo, base, dimero, fator_conv)[0]
                            eid = calcula_energia(metodo, base, dimero, fator_conv)[1]
                            edi = calcula_energia(metodo, base, dimero, fator_conv)[2]
                            exh = calcula_energia(metodo, base, dimero, fator_conv)[3]
                            esap = calcula_energia(metodo, base, dimero, fator_conv)[4]

                            eelst.append(eel)
                            eind.append(eid)
                            edisp.append(edi)
                            eexch.append(exh)
                            esapt.append(esap)

                            if esapt[-1] - esapt[-2] < 0:
                                nova_dist1 = round(dist+0.1, 2)
                                conta_min_energia += 1
                            elif esapt[-1] - esapt[-2] < 0 and conta_min_energia >= 2:
                                nova_dist1 = round(dist+0.05, 2)
                                conta_min_energia += 1
                            else:
                                nova_dist1 = round(dist + inc_ini, 2)
                            distancias.append(nova_dist1)
                            dist = nova_dist1

                conta_min_energia = 0

                if metodo != 'sherrill_gold_standard':
                    sitio_inte = geo.replace('sapt.xyz', '').replace('amonia', '')
                    nome_arq_out = metodo +  sitio_inte + base + '.dat'
                    cria_arquivo(nome_arq_out, metodo, distancias, en_sem_cp, en_com_cp,
                                 eelst, eind, edisp, eexch, esapt, tipo_conversao_energia)
                    move_arquivo(nome_arq_out, gas_nobre, metodo, base)
                else:
                    sitio_inte = geo.replace('sapt.xyz', '').replace('amonia', '')
                    nome_arq_out = metodo +  sitio_inte + '.dat'
                    cria_arquivo(nome_arq_out, metodo, distancias, en_sem_cp, en_com_cp,
                                 eelst, eind, edisp, eexch, esapt, tipo_conversao_energia)
                    move_arquivo(nome_arq_out, gas_nobre, metodo, base='')

            move_diretorio(gas_nobre, metodo, base)
