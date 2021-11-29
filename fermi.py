import numpy as np

# Variable declaration
k = 1.38e-023  # Boltzmann constant, J/K
h = 6.626e-034  # Planck's constant, Js
eV = 1.6e-019  # Joule equivalent of 1 eV
m_0 = 9.1e-031  # Rest mass of an electron, kg

# Input variables ######################################################
m_e = 0.063 * m_0  # Effective mass of electron, kg
m_h = 0.51 * m_0  # Effective mass of hole, kg

E_g = 1.42  # Energy gap, eV  // E_c - E_v = E_g => E_v = 0 => E_c = E_g
N_d = 1e+17  # Donor concentration, cm^-3
N_a = 1e+18  # Acceptor concentration, cm^-3

T = 300  # Room temperature, kelvin

# E_ion_d = 0.006  # Energy ionization for donor, eV
# E_d = E_g - E_ion_d  # Donor energy

# E_ion_a = 0.03  # Energy ionization for acceptor, eV
# E_a = E_ion_a  # Acceptor energy


def intrinsic():
    E_f = E_g / 2 + 3 / 4 * k * T * np.log(m_h / m_e)
    return E_f


def concentration():
    N_c = pow((2 * m_e * k * T * np.pi), 3 / 2) * 2 / pow((h), 3) / 1e+6 # из м -> см
    N_v = pow((2 * m_h * k * T * np.pi), 3 / 2) * 2 / pow((h), 3) / 1e+6
    n_i = np.sqrt(N_c * N_v) * np.exp(-E_g / (2 * k * T / eV))
    return n_i, N_c, N_v


def doped():
    E_i = intrinsic()

    n_i, N_c, N_v = concentration()

    print('\nSemiconductor D-type')
    E_f = k * T * np.log(N_d / n_i) / eV + E_i
    print(f'E_f = {E_f} eV')

    print('\nSemiconductor A-type')
    E_f = - k * T * np.log(N_a / n_i) / eV + E_i
    print(f'E_f = {E_f} eV\n')

    print(f'\tn_i = {n_i} cm^-3\n\tN_c = {N_c} cm^-3\n\tN_v = {N_v} cm^-3')


def main():
    print(f'\nIntrinsic fermi-level = {intrinsic()} eV')
    doped()


if __name__ == '__main__':
    main()
