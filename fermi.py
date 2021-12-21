import numpy as np

# Variable declaration
k = 1.38e-023  # Boltzmann constant, J/K
h = 6.626e-034  # Planck's constant, Js
eV = 1.6e-019  # Joule equivalent of 1 eV
m_0 = 9.1e-031  # Rest mass of an electron, kg
ksi = 8.87e-12  # V/m^2

# Input variables ######################################################
N_d = 5e17  # Donor concentration, cm^-3
N_a = 8e16  # Acceptor concentration, cm^-3
T = 250  # Room temperature, kelvin

m_e_GaP = 1.12 * m_0  # Effective mass of electron, kg
m_h_GaP = 0.79 * m_0  # Effective mass of hole, kg
E_g_GaP = 2.26  # Energy gap, eV  // E_c - E_v = E_g => E_v = 0 => E_c = E_g
ksi_GaP = 11.1
E_i_GaP = 3.8  # Electron affinity Xe, eV

m_e_InP = 0.8 * m_0
m_h_InP = 0.6 * m_0
E_g_InP = 1.344
ksi_InP = 12.5
E_i_InP = 4.38

m_e_Si = 0.36 * m_0
m_h_Si = 0.81 * m_0
E_g_Si = 1.12
ksi_Si = 11.7
E_i_Si = 4.38

E_ion_d_P = 0.045  # Energy ionization for donor, eV
E_ion_a_P = 0.045


def E_d(E_g, E_ion_d):
    return E_g - E_ion_d


def E_a(E_ion_a):
    return E_ion_a


# def phi(E_g, N_c, E_f_d, N_v, E_f_a, ksi_material):
#     phi_n = - (eV * (N_d - N_c * np.exp((E_f_d - E_g) / (k * T)))) / (ksi * ksi_material)
#     phi_p = - (eV * (N_a - N_v * np.exp((- E_f_a + E_g) / (k * T)))) / (ksi * ksi_material)
#     phi_full = phi_p + phi_n
#     return phi_n, phi_p, phi_full


def opz(delta_phi, ksi_material):
    opz_n = np.sqrt(2 * delta_phi * ksi * ksi_material / (eV * N_d * 1e+6) * (N_a / (N_a + N_d)))
    opz_p = np.sqrt(2 * delta_phi * ksi * ksi_material / (eV * N_a * 1e+6) * (N_d / (N_a + N_d)))
    opz_full = opz_n + opz_p
    return opz_n, opz_p, opz_full


def W(delta_phi, ksi_material):
    return pow(2 * delta_phi * ksi * ksi_material / (eV * N_d * 1e+6), 1 / 2)  # cm^-3 -> m^-3


def intrinsic(E_g, m_h, m_e):
    E_f = E_g / 2 + 3 / 4 * k * T * np.log(m_h / m_e)
    return E_f


def concentration(E_g, m_h, m_e):
    N_c = pow((2 * m_e * k * T * np.pi), 3 / 2) * 2 / pow((h), 3) / 1e+6  # из м -> см
    N_v = pow((2 * m_h * k * T * np.pi), 3 / 2) * 2 / pow((h), 3) / 1e+6
    n_i_2 = np.sqrt(N_c * N_v) * np.exp(-E_g / (2 * k * T / eV))
    return n_i_2, N_c, N_v


def p_n_concentration(N_c, N_v, E_f_d, E_f_a, E_f):
    n_n = N_c * np.exp((-E_f_d + E_f) / (k * T / eV))
    p_p = N_v * np.exp((E_f_a - E_f) / (k * T / eV))
    return n_n, p_p


def doped_concentration(E_f_d, E_f_a):
    N_a_minus = N_a / (1 + 4 * np.exp((E_a - E_f_a) / (k * T)))
    N_d_plus = N_d / (1 + 1 / 2 * np.exp((E_d - E_f_d) / (k * T)))
    return N_a_minus, N_d_plus


def doped(E_g, m_h, m_e):
    E_i = intrinsic(E_g, m_h, m_e)
    n_i, N_c, N_v = concentration(E_g, m_h, m_e)
    E_f_d = k * T * np.log(N_d / n_i) / eV + E_i
    E_f_a = - k * T * np.log(N_a / n_i) / eV + E_i
    delta_phi = E_f_d - E_f_a

    return E_f_d, E_f_a, delta_phi


def energy_lvl_at_H(m, ksi_material):
    return m * 13.6 / (m_0 * pow(ksi_material, 2))


def main():
    E_i = intrinsic(E_g_Si, m_h_Si, m_e_Si)
    print(f'\nIntrinsic fermi-level = {E_i} eV')

    E_f_d, E_f_a, delta_phi = doped(E_g_Si, m_h_Si, m_e_Si)
    print(f'E_f_d = {E_f_d} eV')
    print(f'E_f_a = {E_f_a} eV')
    print(f'Delta PHI = {delta_phi}\n')
    E_H_e, E_H_h = energy_lvl_at_H(m_e_Si, ksi_Si), energy_lvl_at_H(m_h_Si, ksi_Si)
    print(f'E_H_e = {E_g_Si - E_H_e}, E_H_h = {E_H_h}\n')

    n_i_2, N_c, N_v = concentration(E_g_Si, m_h_Si, m_e_Si)
    print(f'\tn_i^2 = {n_i_2} cm^-3\n\tN_c = {N_c} cm^-3\n\tN_v = {N_v} cm^-3\n')

    # n_n, p_p = p_n_concentration(N_c, N_v, E_f_d, E_f_a, E_i)
    # print(f'\tn_n = {n_n} cm^-3\n\tp_p = {p_p} cm^-3\n')

    # delta_phi = 0.15
    # w = W(delta_phi, ksi_Si)
    # print(f'W = {w} m')

    # phi_n, phi_p, phi_full = phi(N_c, E_f_d, N_v, E_f_a, ksi_Si)
    # print(f'PHI: p - {phi_p}, n - {phi_n}, full - {phi_full}')
    #
    opz_n, opz_p, opz_full = opz(delta_phi, ksi_Si)
    print(f'Wp = {opz_p * 100} cm, Wn = {opz_n * 100} cm, W = {opz_full * 100} cm')


if __name__ == '__main__':
    main()
