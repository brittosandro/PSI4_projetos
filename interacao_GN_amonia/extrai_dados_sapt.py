import numpy as np
import os


def caminho(diretorio, arquivo):
    return os.path.join(diretorio, arquivo)

sitio1_mol1 = {}
sitio2_mol1 = {}
sitio3_mol1 = {}

sitio1 = 's1'
sitio2 = 's2'
sitio3 = 's3'

molecula1 = 'Ar'

diretorio_corrente = '.'
for dir, subdirs, arqs in os.walk(diretorio_corrente):
    for arq in arqs:
        if sitio1 in caminho(dir, arq):
            if molecula1 in caminho(dir, arq):
                metodo_base = caminho(dir, arq).replace('./', '').split('/')[1].replace('_' + sitio1 + '_', '/').replace('.dat', '')
                #print(metodo_base)
                dados = np.loadtxt(caminho(dir, arq), comments='#')
                dist = dados[:, 0]
                esapt = dados[:, 5]
                indice_menor_energia_sapt = list(esapt).index(min(esapt))
                #print('Distancia    Minimo SAPT   indice SAPT')
                #print(f'{dist[indice_menor_energia_sapt]}         {min(esapt)}         {indice_menor_energia_sapt}')
                #print()
                sitio1_mol1[metodo_base] = dist[indice_menor_energia_sapt], min(esapt)

        if sitio2 in caminho(dir, arq):
            if molecula1 in caminho(dir, arq):
                metodo_base = caminho(dir, arq).replace('./', '').split('/')[1].replace('_' + sitio2 +'_', '/').replace('.dat', '')
                #print(metodo_base)
                dados = np.loadtxt(caminho(dir, arq), comments='#')
                dist = dados[:, 0]
                esapt = dados[:, 5]
                indice_menor_energia_sapt = list(esapt).index(min(esapt))
                #print('Distancia    Minimo SAPT   indice SAPT')
                #print(f'{dist[indice_menor_energia_sapt]}         {min(esapt)}         {indice_menor_energia_sapt}')
                #print()
                sitio2_mol1[metodo_base] = dist[indice_menor_energia_sapt], min(esapt)

        if sitio3 in caminho(dir, arq):
            if molecula1 in caminho(dir, arq):
                metodo_base = caminho(dir, arq).replace('./', '').split('/')[1].replace('_' + sitio3 + '_', '/').replace('.dat', '')
                #print(metodo_base)
                dados = np.loadtxt(caminho(dir, arq), comments='#')
                dist = dados[:, 0]
                esapt = dados[:, 5]
                indice_menor_energia_sapt = list(esapt).index(min(esapt))
                #print('Distancia    Minimo SAPT   indice SAPT')
                #print(f'{dist[indice_menor_energia_sapt]}         {min(esapt)}         {indice_menor_energia_sapt}')
                #print()
                sitio3_mol1[metodo_base] = dist[indice_menor_energia_sapt], min(esapt)

print(f'Sitio 1 {molecula1}')
print(sitio1_mol1)
print()

print(f'Sitio 2 {molecula1}')
print(sitio2_mol1)
print()

print(f'Sitio 3 {molecula1}')
print(sitio3_mol1)
print()
