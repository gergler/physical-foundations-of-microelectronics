import numpy as np
import matplotlib.pyplot as plt

# constants
h_Plank_SI = 1.054 * 10 ** (-34)  # Дж*с
m_electron_SI = 9.109 * 10 ** (-31)  # кг
lattice_parameter_Silicium_SI = 0.5 * (10 ** (-9))  # м
eV = 1.602 * 10 ** (-19)  # Дж

# m = m_electron_SI
# a = lattice_parameter_Silicium_SI
# b = 1/10*a
# U_0_eV = 30
# U_0 = U_0_eV*one_eV_SI

m = 0.49 * m_electron_SI
a = 10 * (10 ** (-9))
b = 5 * (10 ** (-9))
U_0_eV = 0.58
U_0 = U_0_eV * eV

N = 100
E_0 = U_0 - U_0 / N  # first non-zero approximation
accuracy = U_0 / 10000000  # accuracy with to find energy_beg, energy_end


def dichotomy(func, args, left_border, right_border, epsilon):
    a = min(left_border, right_border)
    b = max(left_border, right_border)
    left_func_val = func(a, args)
    iteration_counter = 0
    max_iteration_count = 100

    if (left_func_val * func(b, args) > 0):
        print("Not correct params for dichotomy method")
        return
    while ((b - a) > epsilon and iteration_counter < max_iteration_count):
        iteration_counter += 1
        middle_func_val = func((b + a) / 2, args)
        if (left_func_val * middle_func_val <= 0):
            b = (b + a) / 2
        else:
            left_func_val = middle_func_val
            a = (b + a) / 2
    if (iteration_counter >= max_iteration_count):
        print("Too much iterations")
        return
    return (b + a) / 2


# Periodic potential:
# U = 0,    n*(a+b) < x < n*(a+b) + a
# U = U_0,  n*(a+b) + a <= (n+1)*(a + b)
# b - тонкий высокий барьер, а - межбарьерное расстояние


# left part of equation
def left_equation_part(a, b, alpha, beta):
    return ((pow(alpha, 2) - beta * beta) / 2 / beta / alpha * np.sinh(alpha * b) * np.sin(beta * a) + np.cosh(
        alpha * b) * np.cos(beta * a))


def energy_equation(E, args):
    U_0, a, b, m, wave_number_k = args
    constant = pow(2 * m / pow(h_Plank_SI, 2), 1 / 2)
    alpha = constant * pow((U_0 - E), 1 / 2)
    beta = constant * pow(E, 1 / 2)
    return left_equation_part(a, b, alpha, beta) - np.cos(wave_number_k * (a + b))


def energy(a, b, U_0, m, E_0, accuracy):
    constant = pow(2 * m / pow(h_Plank_SI, 2), 1 / 2)
    E_start, E_end = E_0, E_0
    step = U_0 / 100000
    cur_E = E_start
    process_started = False
    while (True):
        alpha = constant * pow((U_0 - cur_E), 1 / 2)
        beta = constant * pow(cur_E, 1 / 2)
        if (not process_started):
            if (left_equation_part(a, b, alpha, beta) < 1.0):
                cur_E -= step
            else:
                process_started = True
        elif (left_equation_part(a, b, alpha, beta) > 1.0):
            cur_E += step
        elif (left_equation_part(a, b, alpha, beta) < 1.0):
            if (step < accuracy):
                E_start = cur_E
                break
            step /= 2.
            cur_E -= step
        else:
            E_start = cur_E
            break
    minus_one_achieved = False
    while (True):
        alpha = constant * pow((U_0 - cur_E), 1 / 2)
        beta = constant * pow(cur_E, 1 / 2)
        if (left_equation_part(a, b, alpha, beta) > -1.0):
            if (minus_one_achieved and step < accuracy):
                E_end = cur_E
                break
            cur_E += step
        elif (left_equation_part(a, b, alpha, beta) < -1.0):
            minus_one_achieved = True
            step /= 2.
            cur_E -= step
        else:
            E_end = cur_E
            break
    return E_start, E_end


def main():
    energy_beg, energy_end = energy(a, b, U_0, m, E_0, accuracy)
    print("For first Brullien zone")
    print("E_min [J] = ", energy_beg)
    print("E_max [J] = ", energy_end)

    necessary_k = 0.5 * np.pi / (a + b)
    found_energy = dichotomy(energy_equation, (U_0, a, b, m, necessary_k), energy_beg, energy_end, accuracy)
    print("Found E [J] = ", found_energy)
    print("Found E [eV] = ", found_energy / eV)

    k_N = 100
    k_step = 2 * np.pi / (a + b) / k_N
    k_arr = np.arange(-1 * np.pi / (a + b), np.pi / (a + b) + k_step, k_step)
    energy_arr = [dichotomy(energy_equation, (U_0, a, b, m, k), energy_beg, energy_end, accuracy) for k in k_arr]

    for i in range(len(energy_arr)):
        if (energy_arr[i] == None):
            if (i == 0):
                energy_arr[i] = energy_arr[i + 1]
            elif (i == k_N):
                energy_arr[i] = energy_arr[i - 1]
            else:
                energy_arr[i] = (energy_arr[i - 1] + energy_arr[i + 1]) / 2

    energy_arr = [E / eV for E in energy_arr]
    kl_arr = [cur_k * (a + b) for cur_k in k_arr]

    fig, axes = plt.subplots(nrows=1, ncols=1)

    axes.plot(kl_arr, energy_arr, lw=2, color='b', alpha=1)
    axes.set_xlabel("k*(a+b)")
    axes.set_ylabel("Energy")
    axes.set_title("Energy of the first Brullien zone")

    fig.set_figwidth(6)
    fig.set_figheight(6)
    plt.show()


if __name__ == '__main__':
    main()
