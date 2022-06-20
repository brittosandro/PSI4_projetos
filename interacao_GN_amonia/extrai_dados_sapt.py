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
                                   + sitio1 + '_', '/').replace('.dat', '')
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

molecula1 = 'Ar'

sitio1 = 's1'
sitio2 = 's2'
sitio3 = 's3'

sitio1_mol1 = {}
sitio2_mol1 = {}
sitio3_mol1 = {}

diretorio_corrente = '.'
for dir, subdirs, arqs in os.walk(diretorio_corrente):
    for arq in arqs:
        if sitio1 in caminho(dir, arq):
            if molecula1 in caminho(dir, arq):
                metodo_base = pega_metodo_base(sitio1)
                dist, esapt = distancia_energia_sapt(caminho(dir, arq))
                ind_menor_en = indice_para_menor_energia_sapt(esapt)
                sitio1_mol1[metodo_base] = dist[ind_menor_en], min(esapt)

        if sitio2 in caminho(dir, arq):
            if molecula1 in caminho(dir, arq):
                metodo_base = caminho(dir, arq).replace('./', '').split('/')[1].replace('_' + sitio2 +'_', '/').replace('.dat', '')
                metodo_base = pega_metodo_base(sitio1)
                dist, esapt = distancia_energia_sapt(caminho(dir, arq))
                ind_menor_en = indice_para_menor_energia_sapt(esapt)
                sitio2_mol1[metodo_base] = dist[ind_menor_en], min(esapt)

        if sitio3 in caminho(dir, arq):
            if molecula1 in caminho(dir, arq):
                metodo_base = caminho(dir, arq).replace('./', '').split('/')[1].replace('_' + sitio3 + '_', '/').replace('.dat', '')
                metodo_base = pega_metodo_base(sitio1)
                dist, esapt = distancia_energia_sapt(caminho(dir, arq))
                ind_menor_en = indice_para_menor_energia_sapt(esapt)
                sitio3_mol1[metodo_base] = dist[ind_menor_en], min(esapt)

print(f'Sitio 1 {molecula1}')
print('Dicionario NÃO ORDENADO')
print(sitio1_mol1)
print()


print('Dicionario ORDENADO')
sitio1_mol1_ordenado = {}
for item in sorted(sitio1_mol1, key=sitio1_mol1.get):
    sitio1_mol1_ordenado[item] = sitio1_mol1[item]
print(sitio1_mol1_ordenado)
print()


'''
print(f'Sitio 2 {molecula1}')
print(sitio2_mol1)
print()

print(f'Sitio 3 {molecula1}')
print(sitio3_mol1)
print()
'''

