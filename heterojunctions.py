import numpy as np

# Variable declaration
h = 6.626e-027  # Planck's constant, Erg*s
e = 4.803e-10
m_0 = 9.1e-028  # Rest mass of an electron, g
erg = 1e+11  # in 1 erg ??? strange

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


def plasm_freq(a):
    return np.sqrt(4 * np.pi * (32 / pow(a, 3)) * pow(e, 2) / m_0)


def energy(a, epsilon):
    return h * plasm_freq(a) / np.sqrt(epsilon - 1) * erg


def main():
    Si = energy(a_Si, epsilon_Si)
    Ge = energy(a_Ge, epsilon_Ge)
    print(f'\nE_g_Si* = {Si} eV')
    print(f'E_g_Ge* = {Ge} eV')
    print(f'E_v* = E_v = {np.abs(Si - Ge)/2} eV')
    GaAs = energy(a_GaAs, epsilon_GaAs)
    AlAs = energy(a_AlAs, epsilon_AlAs)
    print(f'\nE_g_GaAs* = {GaAs} eV')
    print(f'E_g_AlAs* = {AlAs} eV')
    print(f'E_v* = E_v = {np.abs(GaAs - AlAs)/2} eV')


if __name__ == '__main__':
    main()
