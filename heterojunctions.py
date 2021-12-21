import numpy as np

# Variable declaration
h = 4.135e-15  # Planck's constant, Ev * c
e = 4.803e-10
m_0 = 9.1e-28  # Rest mass of an electron, g
erg = 1.6e-12  # 1 eV

# Input variables ######################################################
a_Si = 5.431e-8  # Lattice constant, cm
epsilon_Si = 11.7  # Dielectric constant
E_g_Si = 1.17  # eV

a_Ge = 5.658e-8  # Lattice constant, cm
epsilon_Ge = 16.2  # Dielectric constant
E_g_Ge = 0.65  # eV

a_GaAs = 5.65e-8
epsilon_GaAs = 10.9
E_g_GaAs = 4.97

a_AlAs = 5.66e-8
epsilon_AlAs = 8.16
E_g_AlAs = 5.82

a_GaP = 5.4505e-8
epsilon_GaP = 11.1
E_g_GaP = 2.26

a_InP = 5.8687e-8
epsilon_InP = 12.5
E_g_InP = 1.344

a_GaSb = 6.09e-8
epsilon_GaSb = 15.7
E_g_GaSb = 0.726

a_InSb = 6.479e-8
epsilon_InSb = 16.8
E_g_InSb = 0.726

def plasm_freq(a):
    return np.sqrt(4 * np.pi * (32 / pow(a, 3)) * pow(e, 2) / m_0)  # СГС


def energy(a, epsilon):
    return h * plasm_freq(a) / np.sqrt(epsilon - 1) * erg / 1e-11


def main():
    Si = energy(a_Si, epsilon_Si)
    Ge = energy(a_Ge, epsilon_Ge)
    print(f'\nE_g_Si* = {Si} eV')
    print(f'E_g_Ge* = {Ge} eV')
    print(f'E_v* = E_v = {np.abs(Si - Ge)/2} eV')
    # GaAs = energy(a_GaAs, epsilon_GaAs)
    # AlAs = energy(a_AlAs, epsilon_AlAs)
    # print(f'\nE_g_GaAs* = {GaAs} eV')
    # print(f'E_g_AlAs* = {AlAs} eV')
    # print(f'E_v* = E_v = {np.abs(GaAs - AlAs)/2} eV')

    # GaP = energy(a_GaP, epsilon_GaP)
    # InP = energy(a_InP, epsilon_InP)
    # print(f'\nE_g_GaP* = {GaP} eV')
    # print(f'E_g_InP* = {InP} eV')
    # print(f'E_v* = E_v = {np.abs(GaP - InP) / 2} eV')

    # GaP = energy(a_GaP, epsilon_GaP)
    # print(f'\nE_g_GaP* = {GaP} eV')
    # GaAs = energy(a_GaAs, epsilon_GaAs)
    # print(f'E_g_GaAs* = {GaAs} eV')
    # print(f'E_v* = E_v = {np.abs(GaP - GaAs) / 2} eV')

    GaSb = energy(a_GaSb, epsilon_GaSb)
    print(f'\nE_g_GaSb* = {GaSb} eV')
    InSb = energy(a_InSb, epsilon_InSb)
    print(f'E_g_InSbs* = {InSb} eV')
    print(f'E_v* = E_v = {np.abs(GaSb - InSb) / 2} eV')


if __name__ == '__main__':
    main()
