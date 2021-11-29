import numpy as np
from fompy import constants, materials, models
from scipy.optimize import fsolve


def check_parameters(parameters, Nc, Nv):
    message = 'ok'

    if (parameters['E_as'] > parameters['E_gap']):
        message = 'Error! The surface states do NOT fall into the energy gap of the semiconductor'
    if (parameters['E_out'] * 3.3 * 1e-5 > parameters['N_as'] * constants.e):
        message = 'Error! The external field is larger than the field created by surface acceptors'
    if (parameters['N_d0'] > Nc or parameters['N_d0'] > Nv):
        message = 'Error! The donor concentration is greater than the concentration of the intrinsic densities of states in the valence and conduction bands'

    return message


def equation_for_phi_left(x, parameters):
    x_erg = x * constants.eV
    return np.sqrt(parameters['epsilon'] * x_erg * parameters['N_d0'] / (2 * np.pi * pow(constants.e, 2)))


def equation_for_phi_right(x, parameters):
    x_erg = x * constants.eV
    return (parameters['N_as'] * (1 / (1 + np.exp((parameters['E_as'] + x_erg - parameters['E_f'])
                                                  / (constants.k * parameters['T'])))) + parameters['E_out'] / (
                        4 * np.pi * constants.e))


def equation_for_phi(x, parameters):
    return equation_for_phi_left(x, parameters) - equation_for_phi_right(x, parameters)


def W(phi, parameters):
    phi_erg = phi * constants.eV
    return np.sqrt(parameters['epsilon'] * phi_erg / (parameters['N_d0'] * 2 * np.pi * pow(constants.e, 2)))


def solve_equation_find_phi(parameters):
    x_0 = 0.001
    phi = fsolve(equation_for_phi, x_0, args=parameters)
    return phi[0]


def data_for_graph(phi, W, parameters):  # phi [eV], W [cm]
    N = 30  # knot number
    h = W * 2 / N  # step

    # parabola: ax ^ 2 + bx + c
    c = phi
    a = c / (W ** 2)
    b = -2 * W * a

    x_s = []  # Coordinate
    E_f_s = []  # Fermi energy
    E_v_s = []  # Valence band ceiling
    E_c_s = []  # Conduction band bottom
    E_d_s = []  # Donor energy
    E_as_s = []  # Energy of surface acceptors

    E_f = parameters['E_f'] / constants.eV
    E_gap = parameters['E_gap'] / constants.eV
    E_d = parameters['E_d'] / constants.eV
    E_as = parameters['E_as'] / constants.eV

    for i in range(N + 1):
        x_s.append(i * h)
        E_f_s.append(E_f)
        E_as_s.append(E_as)

        if x_s[i] > W:  # flat zone
            E_v_s.append(0)
            E_c_s.append(E_gap)
            E_d_s.append(E_d)

        else:  # zone with parabola: flat + parabola bend
            bend = a * x_s[i] ** 2 + b * x_s[i] + c
            E_v_s.append(bend)
            E_c_s.append(E_gap + bend)
            E_d_s.append(E_d + bend)

    return x_s, E_f_s, E_v_s, E_c_s, E_d_s, E_as_s


def calculate(parameters):
    semiconductor = models.Semiconductor(parameters['m_e'] * constants.me, parameters['m_h'] * constants.me,
                                         parameters['E_gap'] * constants.eV, eps=parameters['epsilon'], chi=None)
    T = parameters['T']
    message = check_parameters(parameters, semiconductor.Nc(T), semiconductor.Nv(T))
    results = dict(message='', x_s=0, E_f_s=0, E_v_s=0, E_c_s=0, E_d_s=0, E_as_s=0, phi=0, W=0)
    if message == 'ok':
        results['message'] = message

        # СГС
        parameters['E_gap'] = parameters['E_gap'] * constants.eV
        parameters['E_d'] = parameters['E_d'] * constants.eV
        parameters['E_as'] = parameters['E_as'] * constants.eV
        parameters['E_out'] = parameters['E_out'] * 3.3 * 1e-5

        try:
            parameters['E_f'] = semiconductor.fermi_level(T)
            results['phi'] = solve_equation_find_phi(parameters)  # eV
            results['W'] = W(results['phi'], parameters)  # cm
            results['x_s'], results['E_f_s'], results['E_v_s'], results['E_c_s'], results['E_d_s'], \
            results['E_as_s'] = data_for_graph(results['phi'], results['W'], parameters)
        except Exception:
            message = 'Error! Incorrect data'
            results['message'] = message
            return results

    else:
        results['message'] = message
    return results
