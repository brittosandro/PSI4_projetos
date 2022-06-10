import numpy as np
import os


def caminho(diretorio, arquivo):
    return os.path.join(diretorio, arquivo)

sitio1_Ar = {}
sitio2_Ar = {}
sitio3_Ar = {}

diretorio_corrente = '.'
for dir, subdirs, arqs in os.walk(diretorio_corrente):
    for arq in arqs:
        if 's1' in caminho(dir, arq):
            if 'Ar' in caminho(dir, arq):
                metodo_base = caminho(dir, arq).replace('./', '').split('/')[1].replace('_s1_', '/').replace('.dat', '')
                print(metodo_base)
                dados = np.loadtxt(caminho(dir, arq), comments='#')
                dist = dados[:, 0]
                esapt = dados[:, 5]
                indice_menor_energia_sapt = list(esapt).index(min(esapt))
                print('Distancia    Minimo SAPT   indice SAPT')
                print(f'{dist[indice_menor_energia_sapt]}         {min(esapt)}         {indice_menor_energia_sapt}')
                print()
                sitio1_Ar[metodo_base] = dist[indice_menor_energia_sapt], min(esapt)

        if 's2' in caminho(dir, arq):
            if 'Ar' in caminho(dir, arq):
                metodo_base = caminho(dir, arq).replace('./', '').split('/')[1].replace('_s2_', '/').replace('.dat', '')
                print(metodo_base)
                dados = np.loadtxt(caminho(dir, arq), comments='#')
                dist = dados[:, 0]
                esapt = dados[:, 5]
                indice_menor_energia_sapt = list(esapt).index(min(esapt))
                print('Distancia    Minimo SAPT   indice SAPT')
                print(f'{dist[indice_menor_energia_sapt]}         {min(esapt)}         {indice_menor_energia_sapt}')
                print()
                sitio2_Ar[metodo_base] = dist[indice_menor_energia_sapt], min(esapt)

        if 's3' in caminho(dir, arq):
            if 'Ar' in caminho(dir, arq):
                metodo_base = caminho(dir, arq).replace('./', '').split('/')[1].replace('_s3_', '/').replace('.dat', '')
                print(metodo_base)
                dados = np.loadtxt(caminho(dir, arq), comments='#')
                dist = dados[:, 0]
                esapt = dados[:, 5]
                indice_menor_energia_sapt = list(esapt).index(min(esapt))
                print('Distancia    Minimo SAPT   indice SAPT')
                print(f'{dist[indice_menor_energia_sapt]}         {min(esapt)}         {indice_menor_energia_sapt}')
                print()
                sitio3_Ar[metodo_base] = dist[indice_menor_energia_sapt], min(esapt)

print('Sitio 1 Ar')
print(sitio1_Ar)
print()

print('Sitio 2 Ar')
print(sitio2_Ar)
print()

print('Sitio 3 Ar')
print(sitio3_Ar)
print()
