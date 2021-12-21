import numpy as np

# h = 4.135e-15  # Planck's constant, Ev * c
# e = 4.803e-10
# m_0 = 9.1e-28  # Rest mass of an electron, g
# erg = 1.6e-12  # 1 eV
# k = 1.38e-16

k = 1.38e-023  # Boltzmann constant, J/K
h = 6.626e-034  # Planck's constant, Js
e = 1.6e-019  # Joule equivalent of 1 eV
m_0 = 9.1e-031  # Rest mass of an electron, kg

T = 250

N_d = 5e17  # Donor concentration, cm^-3
N_a = 8e16  # Acceptor concentration, cm^-3

mu_n = 1400
mu_p = 450

m_n = 0.36*m_0
m_p = 0.81*m_0

def D(mu):
    return mu*k*T/e

def tau(m, mu):
    return m*mu/e

def L(D, tau):
    return np.sqrt(D*tau)

print(f'D_n = {D(mu_n)}, tau_n = {tau(m_n, mu_n)}, L_n = {L(D(mu_n), tau(m_n, mu_n))}')
print(f'D_p = {D(mu_p)}, tau_p = {tau(m_p, mu_p)}, L_p = {L(D(mu_p), tau(m_p, mu_p))}')

def I(D, L, N):
    return D*N/L*e

print(I(D(mu_n), L(D(mu_n), tau(m_n, mu_n)), N_a) + I(D(mu_p), L(D(mu_p), tau(m_p, mu_p)), N_d))
