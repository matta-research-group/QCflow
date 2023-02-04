import cclib

def get_data(file):
    '''
    When provied the name of the file, it creates an cclib object

    file: str, the name of the .log file
    '''
    data = cclib.ccopen(file).parse()
    return data

def cal_reorg(name, ion):
    '''
    when provided the name of the oligamer and the type of ion, function returns the reorganisation energy(eV).

    name: str, name of the oligomer in the oligomer dictionary (e.g.'0_2')

    ion: str, type of ion:
        anion = 'an'
        cation = 'cat'
    '''

    #lower case charge,upper case geometry, i.e. n = neutral charge, N = neutral geometry
    if ion == 'cat':
        nN = get_data(f'{name}/{name}_pop_opt_n.log')
        cN = get_data(f'{name}/{name}_ver_c.log')
        cC = get_data(f'{name}/{name}_opt_c.log')
        nC = get_data(f'{name}/{name}_cat_geo_n.log')

        EnN = nN.scfenergies[nN.optstatus==4][0]
        EcN = cN.scfenergies[0]
        EcC = cC.scfenergies[cC.optstatus==4][0]
        EnC = nC.scfenergies[0]

        reorg_en = (EcN-EnN)+(EnC-EcC)

        #reorg_en = (cN.scfenergies[0]-nN.scfenergies[nN.optstatus==4][0])+(nC.scfenergies[0]-cC.scfenergies[cC.optstatus==4][0])

    if ion == 'an':
        nN = get_data(f'{name}/{name}_pop_opt_n.log')
        aN = get_data(f'{name}/{name}_ver_a.log')
        aA = get_data(f'{name}/{name}_opt_a.log')
        nA = get_data(f'{name}/{name}_an_geo_n.log')

        EnN = nN.scfenergies[nN.optstatus==4][0]
        EaN = aN.scfenergies[0]
        EaA = aA.scfenergies[aA.optstatus==4][0]
        EnA = nA.scfenergies[0]

        reorg_en = (EaN-EnN)+(EnA-EaA)

    return reorg_en

#Personal Note: too many vatiables, need to be reduced to increase the speed of calculation

def cal_HOMO(name):
    '''
    when provided the name of the neutral oligamer, function returns the HOMO energy(eV).

    name: str, name of the oligomer in the oligomer dictionary (e.g.'0_2')
    '''
    data = get_data(f'{name}/{name}_pop_opt_n.log')
    HOMO = data.moenergies[0][data.homos[0]]
    return HOMO

def cal_LUMO(name):
    '''
    when provided the name of the neutral oligamer, function returns the LUMO energy(eV).

    name: str, name of the oligomer in the oligomer dictionary (e.g.'0_2')
    '''
    data = get_data(f'{name}/{name}_pop_opt_n.log')
    LUMO = data.moenergies[0][data.homos[0]+1]
    return LUMO

def cal_gap(name):
    '''
    when provided the name of the neutral oligamer, function returns the energy difference between HOMO
    and LUMO (i.e. HOMO-LUMO gap) in eV.

    name: str, name of the oligomer in the oligomer dictionary (e.g.'0_2')
    '''
    gap = cal_LUMO(name) - cal_HOMO(name)
    return gap

def cal_IP(name,kind):
    '''
    when provided the name of the oligamer, and the type of the IP, function returns
    the ionisation energy.

    name: str, name of the oligomer in the oligomer dictionary (e.g.'0_2')

    kind: str, type of IP:
        vertical = 'ver'
        adiabatic = 'adi'
    '''
    if kind == 'ver':
        nN = get_data(f'{name}/{name}_pop_opt_n.log')
        cN = get_data(f'{name}/{name}_ver_c.log')

        EnN = nN.scfenergies[nN.optstatus==4][0]
        EcN = cN.scfenergies[0]

        vert_IP = EcN - EnN
        return vert_IP

    if kind == 'adi':
        nN = get_data(f'{name}/{name}_pop_opt_n.log')
        cC = get_data(f'{name}/{name}_opt_c.log')

        EnN = nN.scfenergies[nN.optstatus==4][0]
        EcC = cC.scfenergies[cC.optstatus==4][0]
        adiabatic_IP = EcC - EnN
        return adiabatic_IP

def cal_EA(name,kind):
    '''
    when provided the name of the oligamer, and the type of the EA, function returns the energy for Electron
    Affinity.

    name: str, name of the oligomer in the oligomer dictionary (e.g.'0_2')

    kind: str, type of EA:
        vertical = 'ver'
        adiabatic = 'adi'
    '''
    if kind == 'ver':
        nA = get_data('0_2/0_2_an_geo_n.log')
        aA = get_data('0_2/0_2_opt_a.log')

        EnA = nA.scfenergies[0]
        EaA = aA.scfenergies[aA.optstatus==4][0]

        vertical_EA = EaA - EnA
        return vertical_EA

    if kind == 'adi':
        aA = get_data('0_2/0_2_opt_a.log')
        nN = get_data('0_2/0_2_pop_opt_n.log')

        EaA = aA.scfenergies[aA.optstatus==4][0]
        EnN = EnN = nN.scfenergies[nN.optstatus==4][0]

        adiabatic_EA = EaA - EnN
        return adiabatic_EA
