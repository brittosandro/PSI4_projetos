###############################################################################
#
# Esse script extrai dados de arquivos de energia total sapt em um diretório.
#
###############################################################################
#
#    O input deve ser colocado em um diretório em que estejam os diretórios
# das respectivas moléculas que foram realizados os cálculos sapt. Por
# exemplo: Existe um diretório chamado calc_Gases_Nobres_Amonia, e dentro
# deste diretório existem vários subdiretórios chamados de Ar, He, Ne. Então
# esse script deverá ser adicionado dentro de calc_Gases_Nobres_Amonia. O que
# o script irá fazer e entrar em cada um dos subdiretórios (Ar, He, Ne)
# ler os arquivos com os mais distintos cálculos sapt e mostrar dicionários
# com todas as moléculas que estiverem no diretório.
#
###############################################################################

import numpy as np
import os


def caminho(diretorio, arquivo):
    '''
    Essa função recebe o diretório e o arquivo com os dados de energia SAPT e
    retorna o caminho de onde os arquivos estão.
    '''
    return os.path.join(diretorio, arquivo)

def pega_metodo_base(nome_sitio):
    '''
    Essa função recebe o nome de um certo sítio de interação
    e retorna o método e a base que as interações foram calculadas,
    no formato: metodo/base.
    '''
    mb =  caminho(dir, arq).replace('./', '').split('/')[1].replace('_'
                                   + nome_sitio + '_', '/').replace('.dat', '')
    return mb

def distancia_energia_sapt(caminho):
    '''
    Essa função recebe um caminho de um arquivo com os dados de distancia
    e energia SAPT para diferentes bases e métodos e retorna uma tupla cujos
    elementos são arrays de distância e energia SAPT.
    '''
    dados =  np.loadtxt(caminho, comments='#')
    distancia = dados[:, 0]
    esapt = dados[:, 5]
    return distancia, esapt

def indice_para_menor_energia_sapt(esapt):
    '''
    Essa função recebe o array de energia SAPT e retorna o valor de índice
    associado a menor energia do array.
    '''
    return list(esapt).index(min(esapt))

def ordena_dicionario(d):
    '''
    Essa função recebe um dicionario cuja chave designa o método e a base do
    cálculo e o valor é uma tupla com a distancia de equilibrio e a energia
    sapt total. A função retorna um dicionario em ordem crescente de energia
    sapt. Por exemplo: {
    'sapt0/aug-cc-pvtz': (3.8, -14.387572),
    'sapt2+/aug-cc-pvtz': (3.8, -14.3323951),
    'sapt2+3/aug-cc-pvtz': (3.8, -14.0726924)
    }
    '''
    dic_ordenado = {}
    for item in sorted(d, key=d.get):
        dic_ordenado[item] = d[item]
    return dic_ordenado

molecula1 = 'Ar'
molecula2 = 'He'
molecula3 = 'Ne'
molecula4 = 'Kr'

sitio1 = 's1'
sitio2 = 's2'
sitio3 = 's3'

sitio1_mol1 = {}
sitio2_mol1 = {}
sitio3_mol1 = {}

sitio1_mol2 = {}
sitio2_mol2 = {}
sitio3_mol2 = {}

sitio1_mol3 = {}
sitio2_mol3 = {}
sitio3_mol3 = {}

sitio1_mol4 = {}
sitio2_mol4 = {}
sitio3_mol4 = {}

diretorio_corrente = '.'
for dir, subdirs, arqs in os.walk(diretorio_corrente):
    if 'sapt' in dir:
        for arq in arqs:
            if 'sapt' in arq:
                if sitio1 in caminho(dir, arq):
                    if molecula1 in caminho(dir, arq):
                        metodo_base = pega_metodo_base(sitio1)
                        dist, esapt = distancia_energia_sapt(caminho(dir, arq))
                        ind_menor_en = indice_para_menor_energia_sapt(esapt)
                        sitio1_mol1[metodo_base] = dist[ind_menor_en], min(esapt)
                        sitio1_mol1_ord = ordena_dicionario(sitio1_mol1)

                    if molecula2 in caminho(dir, arq):
                        metodo_base = pega_metodo_base(sitio1)
                        dist, esapt = distancia_energia_sapt(caminho(dir, arq))
                        ind_menor_en = indice_para_menor_energia_sapt(esapt)
                        sitio1_mol2[metodo_base] = dist[ind_menor_en], min(esapt)
                        sitio1_mol2_ord = ordena_dicionario(sitio1_mol2)

                    if molecula3 in caminho(dir, arq):
                        metodo_base = pega_metodo_base(sitio1)
                        dist, esapt = distancia_energia_sapt(caminho(dir, arq))
                        ind_menor_en = indice_para_menor_energia_sapt(esapt)
                        sitio1_mol3[metodo_base] = dist[ind_menor_en], min(esapt)
                        sitio1_mol3_ord = ordena_dicionario(sitio1_mol3)

                    if molecula4 in caminho(dir, arq):
                        metodo_base = pega_metodo_base(sitio1)
                        dist, esapt = distancia_energia_sapt(caminho(dir, arq))
                        ind_menor_en = indice_para_menor_energia_sapt(esapt)
                        sitio1_mol4[metodo_base] = dist[ind_menor_en], min(esapt)
                        sitio1_mol4_ord = ordena_dicionario(sitio1_mol4)

                if sitio2 in caminho(dir, arq):
                    if molecula1 in caminho(dir, arq):
                        metodo_base = pega_metodo_base(sitio2)
                        dist, esapt = distancia_energia_sapt(caminho(dir, arq))
                        ind_menor_en = indice_para_menor_energia_sapt(esapt)
                        sitio2_mol1[metodo_base] = dist[ind_menor_en], min(esapt)
                        sitio2_mol1_ord = ordena_dicionario(sitio2_mol1)

                    if molecula2 in caminho(dir, arq):
                        metodo_base = pega_metodo_base(sitio2)
                        dist, esapt = distancia_energia_sapt(caminho(dir, arq))
                        ind_menor_en = indice_para_menor_energia_sapt(esapt)
                        sitio2_mol2[metodo_base] = dist[ind_menor_en], min(esapt)
                        sitio2_mol2_ord = ordena_dicionario(sitio2_mol2)

                    if molecula3 in caminho(dir, arq):
                        metodo_base = pega_metodo_base(sitio2)
                        dist, esapt = distancia_energia_sapt(caminho(dir, arq))
                        ind_menor_en = indice_para_menor_energia_sapt(esapt)
                        sitio2_mol3[metodo_base] = dist[ind_menor_en], min(esapt)
                        sitio2_mol3_ord = ordena_dicionario(sitio2_mol3)

                    if molecula4 in caminho(dir, arq):
                        metodo_base = pega_metodo_base(sitio2)
                        dist, esapt = distancia_energia_sapt(caminho(dir, arq))
                        ind_menor_en = indice_para_menor_energia_sapt(esapt)
                        sitio2_mol4[metodo_base] = dist[ind_menor_en], min(esapt)
                        sitio2_mol4_ord = ordena_dicionario(sitio2_mol4)

                if sitio3 in caminho(dir, arq):
                    if molecula1 in caminho(dir, arq):
                        metodo_base = pega_metodo_base(sitio3)
                        dist, esapt = distancia_energia_sapt(caminho(dir, arq))
                        ind_menor_en = indice_para_menor_energia_sapt(esapt)
                        sitio3_mol1[metodo_base] = dist[ind_menor_en], min(esapt)
                        sitio3_mol1_ord = ordena_dicionario(sitio3_mol1)

                    if molecula2 in caminho(dir, arq):
                        metodo_base = pega_metodo_base(sitio3)
                        dist, esapt = distancia_energia_sapt(caminho(dir, arq))
                        ind_menor_en = indice_para_menor_energia_sapt(esapt)
                        sitio3_mol2[metodo_base] = dist[ind_menor_en], min(esapt)
                        sitio3_mol2_ord = ordena_dicionario(sitio3_mol2)

                    if molecula3 in caminho(dir, arq):
                        metodo_base = pega_metodo_base(sitio3)
                        dist, esapt = distancia_energia_sapt(caminho(dir, arq))
                        ind_menor_en = indice_para_menor_energia_sapt(esapt)
                        sitio3_mol3[metodo_base] = dist[ind_menor_en], min(esapt)
                        sitio3_mol3_ord = ordena_dicionario(sitio3_mol3)

                    if molecula4 in caminho(dir, arq):
                        metodo_base = pega_metodo_base(sitio3)
                        dist, esapt = distancia_energia_sapt(caminho(dir, arq))
                        ind_menor_en = indice_para_menor_energia_sapt(esapt)
                        sitio3_mol4[metodo_base] = dist[ind_menor_en], min(esapt)
                        sitio3_mol4_ord = ordena_dicionario(sitio3_mol4)


print(f'Sitio 1 {molecula1}')
print(sitio1_mol1_ord)
print()

print(f'Sitio 2 {molecula1}')
print(sitio2_mol1_ord)
print()

print(f'Sitio 3 {molecula1}')
print(sitio3_mol1_ord)
print('-'*95)

print(f'Sitio 1 {molecula2}')
print(sitio1_mol2_ord)
print()

print(f'Sitio 2 {molecula2}')
print(sitio2_mol2_ord)
print()

print(f'Sitio 3 {molecula2}')
print(sitio3_mol2_ord)
print('-'*95)

print(f'Sitio 1 {molecula3}')
print(sitio1_mol3_ord)
print()

print(f'Sitio 2 {molecula3}')
print(sitio2_mol3_ord)
print()

print(f'Sitio 3 {molecula3}')
print(sitio3_mol3_ord)
print('-'*95)

print(f'Sitio 1 {molecula4}')
print(sitio1_mol4_ord)
print()

print(f'Sitio 2 {molecula4}')
print(sitio2_mol4_ord)
print()

print(f'Sitio 3 {molecula4}')
print(sitio3_mol4_ord)
print('-'*95)
