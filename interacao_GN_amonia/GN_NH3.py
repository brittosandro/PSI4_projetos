from glob import glob
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt
import numpy as np
import psi4
import seaborn as sns
import subprocess
import time


geometrias_amonia = glob('*_sapt.xyz')
#print(geometrias_amonia)

#metodos = ['sapt0', 'sapt2', 'sapt2+']
metodos = ['sapt0', 'sapt2', 'sapt2+', 'sapt2+(3)', 'sapt2+3', ]
bases = ['jun-cc-pvdz']
gases_nobres = ['He', 'Ne', 'Ar', 'Kr']
#gases_nobres = ['He']

for geo in geometrias_amonia:
    with open(geo, 'r') as f:
        str_geo = f.read()

    for gas_nobre in gases_nobres:
        for metodo in metodos:
            for base in bases:
                proc1 = subprocess.run(['mkdir', f'{gas_nobre}_{metodo}_{base}'])

                distancias = np.arange(1.4, 10.1, 0.1)
                # listas de energias com zeros
                eelst = np.zeros((len(distancias)))
                eind = np.zeros((len(distancias)))
                edisp = np.zeros((len(distancias)))
                eexch = np.zeros((len(distancias)))
                esapt = np.zeros((len(distancias)))
                #ect = np.zeros((len(distancias)))

                for i, dist in enumerate(distancias):
                    # Construindo a geometria do Dimero
                    dimero = """
                    0 1 """ + """\n""" + str_geo + """--
                    0 1 """ + """\n""" + gas_nobre + """  """ + str(distancias[i]) + """     0.0000000      0.0000000000
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
                    #ect[i] = psi4.variable('SAPT CT ENERGY') * 27211.4
                    eelst[i] = psi4.variable('SAPT ELST ENERGY') * 27211.4
                    eind[i] = psi4.variable('SAPT IND ENERGY') * 27211.4
                    edisp[i] = psi4.variable('SAPT DISP ENERGY') * 27211.4
                    eexch[i] = psi4.variable('SAPT EXCH ENERGY') * 27211.4
                    esapt[i] = psi4.variable('SAPT TOTAL ENERGY') * 27211.4

                    psi4.core.clean()

                sitio_inte = geo.replace('sapt.xyz', '').replace('amonia', '')
                nome_arq_out = metodo +  sitio_inte + base + '.dat'


                time.sleep(2)
                with open(nome_arq_out, 'w') as f:
                    print(f'# Dist    Energia Eletrostatica       Energia Inducao        Energia Dispercao      Energia EXCH       Energia Sapt        Energia CT', end='\n', file=f)
                    for d, el, ei, ed, eex, es in zip(distancias, np.around(eelst, 7), np.around(eind, 7), np.around(edisp, 7),
                                                   np.around(eexch, 7), np.around(esapt, 7)):
                        print(f'{d:5.2f}  {el:16.7f}   {ei:25.7f}   {ed:20.7f}   {eex:20.7f}   {es:15.7f}', end='\n', file=f)


                proc2 = subprocess.run(['mv', nome_arq_out, f'{gas_nobre}_{metodo}_{base}'])
