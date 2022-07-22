from glob import glob
import numpy as np
import psi4
import subprocess
import time


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
    no_com
    no_reorient
    """
    return input


def cria_diretorio(gas, funcional='sem_funcional', base='sem_base'):
    '''
    Essa função recebe três strings que correspondem ao gás nobre de interesse
    o método e a base que iremos executar os cálculos e cria um novo diretorio
    com as respectivas palavras.
    '''
    if base == 'sem_base' and funcional == 'sem_funcional':
        subprocess.run(['mkdir', f'{gas}'])
    else:
        subprocess.run(['mkdir', f'{gas}_{funcional}_{base}'])


def cria_arquivo(nome, funcional, dist, eint3, eint4, eint5, eint6, eint7):
    with open(nome_arq_out, 'w') as f:
            print(f'# ---------------------------------------------------------------------------------------------------------------', end='\n', file=f)
            print(f'#                        Distancia [angstrom]    |   Energia [meV]', end='\n', file=f)
            print(f'# ---------------------------------------------------------------------------------------------------------------', end='\n', file=f)
            print(f'# Dist    Energia Eletrostatica       Energia Inducao        Energia Dispercao      Energia EXCH     Energia Sapt', end='\n', file=f)
            print(f'# ---------------------------------------------------------------------------------------------------------------', end='\n', file=f)
            for d, el, ei, ed, eex, es in zip(dist, np.around(eint3, 7), np.around(eint4, 7),
                                              np.around(eint5, 7), np.around(eint6, 7), np.around(eint7, 7)):
                print(f'{d:6.2f}  {el:16.9f}   {ei:25.9f}   {ed:20.9f}   {eex:18.9f}   {es:13.9f}', end='\n', file=f)


def move_arquivo(nome, gas, funcional, base):
    subprocess.run(['mv', nome, f'{gas}_{funcional}_{base}'])


def move_diretorio(gas, funcional, base):
    subprocess.run(['mv',  f'{gas}_{funcional}_{base}', f'{gas}/'])


def grac(ecat, eneu, ehomo):
    EI = ecat - eneu
    return EI - ehomo


def calcula_CRAC(funcional, base, molecula):

    psi4.set_options({'freeze_core': 'true',
                      'basis': base,
                      'reference': 'uhf',
                      })
    #geometrias para os calculos
    geo_neutro = """
        0 1 """ + """\n""" + molecula + """
        units angstrom
        symmetry c1
        """
    geo_cation = """
        1 2  """ + """\n""" + molecula + """
        units angstrom
        symmetry c1
        """

    psi4.core.set_output_file('output_neutro.dat', False)
    geo_neutro = psi4.geometry(geo_neutro)
    en_neutro, wfn_neutro = psi4.energy(funcional, molecule=geo_neutro, return_wfn=True)
    HOMO = wfn_neutro.epsilon_a_subset("AO", "ALL").np[wfn_neutro.nalpha()-1]
    #LUMO = wfn_neutro.epsilon_a_subset("AO", "ALL").np[wfn_neutro.nalpha()]

    print('\n')

    #print(f'HOMO = {HOMO}')
    #print(f'LUMO = {LUMO}')
    en_neu = psi4.variable('DFT TOTAL ENERGY')
    print(f'Energia Neutro: {en_neu}')
    psi4.core.clean()

    psi4.core.set_output_file('output_cation.dat', False)
    geo_cation = psi4.geometry(geo_cation)
    en_cation = psi4.energy(funcional, molecule=geo_cation)
    en_cat = psi4.variable('DFT TOTAL ENERGY')
    print(f'Energia Cation: {en_cat}')
    psi4.core.clean()

    g = grac(en_cat, en_neu, HOMO)

    print(f'Funcional = {funcional}')
    print(f'\nGRAC = {g} para Base {base}')

    return g


#psi4.set_memory('3 GB')
#psi4.set_num_threads(2)
#psi4.core.set_output_file('output.dat', False)
#psi4.set_options({'freeze_core': 'true'})

geometrias_amonia = glob('*_sapt.xyz')
#funcionais = ['b3lyp', 'X3LYP', 'cam-b3lyp', 'pbe0', 'WB97X-D', 'WB97X-D3', 'B3LYP-D3MBJ', 'WB97X-D3BJ', 'CAM-B3LYP-D3BJ', 'PBE0-D3BJ',]
funcionais = ['WB97X-D', ]
bases = ['aug-cc-pvtz']
gases_nobres = ['He']

psi4.core.set_output_file('output.dat', False)
for gas_nobre in gases_nobres:
    cria_diretorio(gas_nobre)
    for funcional in funcionais:
        for base in bases:
            cria_diretorio(gas_nobre, funcional, base)
            for geo in geometrias_amonia:
                with open(geo, 'r') as f:
                    str_geo = f.read()

                grac_a = calcula_CRAC(funcional, base, str_geo)
                grac_b = calcula_CRAC(funcional, base, gas_nobre)

                distancias = np.arange(2.8, 6.6, 0.2)
                eelst = np.zeros((len(distancias)))
                eind = np.zeros((len(distancias)))
                edisp = np.zeros((len(distancias)))
                eexch = np.zeros((len(distancias)))
                esapt = np.zeros((len(distancias)))

                for i, dist in enumerate(distancias):
                    psi4.set_options({'freeze_core': 'true',
                                      'reference': 'rhf',
                                      'basis': base,
                                      'scf_type': 'df',
                                      #'maxiter': 10000,
                                      #'D_CONVERGENCE': 1e-6,
                                      #'E_CONVERGENCE': 1e-6,
                                      'SAPT_DFT_FUNCTIONAL': funcional,
                                      'sapt_dft_grac_shift_a': grac_a,
                                      'sapt_dft_grac_shift_b': grac_b})

                    dimero = input_geo(str_geo, gas_nobre, distancias, i)
                    psi4.geometry(dimero)
                    psi4.energy('sapt(dft)')
                    eelst[i] = psi4.variable('SAPT ELST ENERGY') * 27211.4
                    eind[i] = psi4.variable('SAPT IND ENERGY') * 27211.4
                    edisp[i] = psi4.variable('SAPT DISP ENERGY') * 27211.4
                    eexch[i] = psi4.variable('SAPT EXCH ENERGY') * 27211.4
                    esapt[i] = psi4.variable('SAPT TOTAL ENERGY') * 27211.4
                    psi4.core.clean()
                    #time.sleep(5)
                    sitio_inte = geo.replace('sapt.xyz', '').replace('amonia', '')
                    nome_arq_out = funcional +  sitio_inte + base + '.dat'

                    cria_arquivo(nome_arq_out, funcional, distancias, eelst, eind,
                                 edisp, eexch, esapt)
                    move_arquivo(nome_arq_out, gas_nobre, funcional, base)
            move_diretorio(gas_nobre, funcional, base)
