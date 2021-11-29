import tkinter as tk
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import calculations

matplotlib.use("TkAgg")

height = 850
width = 1500

coef_phys_parameters = {'N_d0': 1e12, 'N_as': 1e13, 'E_out': 1e4}


def draw_figure(ax, canvas, results):
    ax.cla()
    ax.set_xlabel("x [cm]")
    ax.set_ylabel("E [eV]")
    ax.grid()

    if results['message'] == 'ok':
        ax.plot(results['x_s'], results['E_f_s'], label='Fermi Energy')
        ax.plot(results['x_s'], results['E_v_s'], label='Valence Band')
        ax.plot(results['x_s'], results['E_c_s'], label='Conduction Band')
        ax.plot(results['x_s'], results['E_d_s'], label='Donor Energy')
        ax.plot(results['x_s'], results['E_as_s'], label='Acceptor Energy')
        ax.axhline(results['phi'], c='k', linestyle='dashed')
        ax.axvline(results['W'], c='k', linestyle='dashed')
        ax.legend(fontsize=10, loc='right')

    output_info(results)
    canvas.draw()


def selection():
    global mat
    if var1.get() == 1:
        mat = 'Si'
    elif var2.get() == 1:
        mat = 'Ge'
    elif var3.get() == 1:
        mat = 'GaAs'
    else:
        mat = 'custom'


def get_calculated_values():
    args = {
        "E_gap": float(E_g.get()),  # Band gap [eV]
        "epsilon": float(epsilon.get()),  # Dielectric permittivity
        "m_h": float(m_h.get()),  # Effective hole mass [m_0]
        "m_e": float(m_e.get()),  # Effective electron mass [m_0]
        "E_d": float(E_d.get()),  # Donors level [eV]
        "N_d0": float(N_d0.get()) * coef_phys_parameters['N_d0'],  # Concentration of donors 10^27 [cm^(-3)]
        "E_as": float(E_as.get()),  # Surface acceptors level [eV]
        "N_as": float(N_as.get()) * coef_phys_parameters['N_as'],  # Concentration of surface acceptors 10^27 [cm^(-3)]
        "T": float(T.get()),  # Temperature [K]
        "E_out": float(E_out.get()) * coef_phys_parameters['E_out'],  # External electric field
        "material": mat  # Material
    }

    return calculations.calculate(args)


def GaAs_calculated():
    GaAs_args = {
        "E_gap": 1.424,  # Band gap [eV]
        "epsilon": 12.9,  # Dielectric permittivity
        "m_h": 0.53,  # Effective hole mass [m_0]
        "m_e": 0.063,  # Effective electron mass [m_0]
        "E_d": 1e-2,  # Donors level [eV]
        "N_d0": 10 * coef_phys_parameters['N_d0'],  # Concentration of donors 10^27 [cm^(-3)]
        "E_as": 1.424 / 2,  # Surface acceptors level [eV]
        "N_as": 10 * coef_phys_parameters['N_as'],  # Concentration of surface acceptors 10^27 [cm^(-3)]
        "T": 300,  # Temperature [K]
        "E_out": 1,  # External electric field
        "material": 'Si'  # Material
    }

    clear()

    E_g.insert(0, "1.423")
    epsilon.insert(0, "12.9")
    m_h.insert(0, "0.53")
    m_e.insert(0, "0.063")
    E_d.insert(0, "0.01")
    N_d0.insert(0, "10")
    E_as.insert(0, "0.712")
    N_as.insert(0, "10")
    T.insert(0, "300")
    E_out.insert(0, "1")

    return calculations.calculate(GaAs_args)


def Ge_calculated():
    Ge_args = {
        "E_gap": 0.661,  # Band gap [eV]
        "epsilon": 16.2,  # Dielectric permittivity
        "m_h": 0.34,  # Effective hole mass [m_0]
        "m_e": 0.22,  # Effective electron mass [m_0]
        "E_d": 1e-2,  # Donors level [eV]
        "N_d0": 10 * coef_phys_parameters['N_d0'],  # Concentration of donors 10^27 [cm^(-3)]
        "E_as": 0.661 / 2,  # Surface acceptors level [eV]
        "N_as": 10 * coef_phys_parameters['N_as'],  # Concentration of surface acceptors 10^27 [cm^(-3)]
        "T": 300,  # Temperature [K]
        "E_out": 1,  # External electric field
        "material": 'Si'  # Material
    }

    clear()

    E_g.insert(0, "0.661")
    epsilon.insert(0, "16.2")
    m_h.insert(0, "0.34")
    m_e.insert(0, "0.22")
    E_d.insert(0, "0.01")
    N_d0.insert(0, "10")
    E_as.insert(0, "0.3305")
    N_as.insert(0, "10")
    T.insert(0, "300")
    E_out.insert(0, "1")

    return calculations.calculate(Ge_args)


def Si_calculated():
    Si_args = {
        "E_gap": 1.12,  # Band gap [eV]
        "epsilon": 11.7,  # Dielectric permittivity
        "m_h": 0.81,  # Effective hole mass [m_0]
        "m_e": 0.36,  # Effective electron mass [m_0]
        "E_d": 1e-2,  # Donors level [eV]
        "N_d0": 10 * coef_phys_parameters['N_d0'],  # Concentration of donors 10^27 [cm^(-3)]
        "E_as": 1.12 / 2,  # Surface acceptors level [eV]
        "N_as": 10 * coef_phys_parameters['N_as'],  # Concentration of surface acceptors 10^27 [cm^(-3)]
        "T": 300,  # Temperature [K]
        "E_out": 1,  # External electric field
        "material": 'Si'  # Material
    }

    clear()

    E_g.insert(0, "1.12")
    epsilon.insert(0, "11.7")
    m_h.insert(0, "0.81")
    m_e.insert(0, "0.36")
    E_d.insert(0, "0.01")
    N_d0.insert(0, "10")
    E_as.insert(0, "0.56")
    N_as.insert(0, "10")
    T.insert(0, "300")
    E_out.insert(0, "1")

    return calculations.calculate(Si_args)


def output_info(results):
    if results['message'] != 'ok':
        tk.Label(window, text=results['message'], fg='orange').grid(row=16, column=5, columnspan=15, rowspan=15,
                                                                    stick='we')
    else:
        tk.Label(window,
                 text=f"fermi level: {results['E_f_s'][0]:.4f} [eV]\nbending of zone PHI: {results['phi']:.4f} [eV]\nspace charge region W: {results['W']:.4f} [cm]",
                 fg='black', bg='white').grid(row=16, column=5, columnspan=15, rowspan=15, stick='we')


def start():
    if mat == 'custom':
        draw_figure(ax, canvas, get_calculated_values())
    if mat == 'Si':
        draw_figure(ax, canvas, Si_calculated())
    if mat == 'Ge':
        draw_figure(ax, canvas, Ge_calculated())
    if mat == 'GaAs':
        draw_figure(ax, canvas, GaAs_calculated())


def clear():
    E_g.delete(0, tk.END)
    E_d.delete(0, tk.END)
    N_d0.delete(0, tk.END)
    E_as.delete(0, tk.END)
    N_as.delete(0, tk.END)
    m_e.delete(0, tk.END)
    m_h.delete(0, tk.END)
    epsilon.delete(0, tk.END)
    T.delete(0, tk.END)
    E_out.delete(0, tk.END)


def help():
    help_text = """HOW TO USE THIS PROGRAM:\n*This program illustrates Fermi level pinning by surface acceptors*\n\n\n- To start the program, choose the type of semiconductor: Si, Ge, GaAs or CUSTOM\n\n*****for CUSTOM set parameters at left blanks, for known types of semiconductor, the data will be set automatically*****\n\n- The START button will configure all the parameters and will draw graph\n\n- If you want to clear all parameters, use the CLEAR button"""
    message_win = tk.Tk()
    message_win.config(bg='lightcyan')
    message_win.title('HELP')
    message_win.geometry(f"1150x300+300+200")
    tk.Label(message_win, text=help_text, fg='black', bg='lightcyan', font="Arial 14", justify=tk.CENTER).pack()
    message_win.mainloop()


##################################################################################
window = tk.Tk()
photo = tk.PhotoImage(file='icon.png')
window.iconphoto(False, photo)
window.config(bg='white')
window.title('Pinning of Fermi Level')
window.geometry(f"{width}x{height}+100+100")
window.resizable(True, True)

text = [['forbidden zone E_g: ', '[eV]'], ['dielectric constant ε: ', ' '], ['effective hole masses m_h: ', 'm_0'],
        ['effective electron masses m_e: ', 'm_0'], ['donor level position E_d: ', '[eV]'],
        ['donor concentration N_d0: ', '10^12[cm^(-3)]'],
        ['energy level position E_as: ', '[eV]'], ['surface acceptor density N_as: ', '10^13[cm^(-3)]'],
        ['temperature T: ', '[K]'], ['external electric field E_out:', '10^4 [V/m]']]

for i in range(len(text)):
    tk.Label(window, text=text[i][0], bg='white', ).grid(row=i, column=1, stick='w')
    tk.Label(window, text=text[i][1], bg='white').grid(row=i, column=3, stick='w')
    window.grid_rowconfigure(i, minsize=30)

E_g = tk.Entry(window)
E_g.insert(0, "5")
E_g.grid(row=0, column=2)

epsilon = tk.Entry(window)
epsilon.insert(0, "12.5")
epsilon.grid(row=1, column=2)

m_h = tk.Entry(window)
m_h.insert(0, "0.5")
m_h.grid(row=2, column=2)

m_e = tk.Entry(window)
m_e.insert(0, "0.5")
m_e.grid(row=3, column=2)

E_d = tk.Entry(window)
E_d.insert(0, "0.5")
E_d.grid(row=4, column=2)

N_d0 = tk.Entry(window)
N_d0.insert(0, "50")
N_d0.grid(row=5, column=2)

E_as = tk.Entry(window)
E_as.insert(0, "2.5")
E_as.grid(row=6, column=2)

N_as = tk.Entry(window)
N_as.insert(0, "50")
N_as.grid(row=7, column=2)

T = tk.Entry(window)
T.insert(0, "300")
T.grid(row=8, column=2)

E_out = tk.Entry(window)
E_out.insert(0, "20")
E_out.grid(row=9, column=2)

tk.Label(window, text='\n\nmaterial: ', bg='white').grid(row=10, column=1, stick='w')
var1 = tk.IntVar()
var2 = tk.IntVar()
var3 = tk.IntVar()
var4 = tk.IntVar()
tk.Checkbutton(window, text='Si', variable=var1, onvalue=1, offvalue=0, command=selection).grid(row=10, column=2,
                                                                                                stick='w')
tk.Checkbutton(window, text='Ge', variable=var2, onvalue=1, offvalue=0, command=selection).grid(row=11, column=2,
                                                                                                stick='w')
tk.Checkbutton(window, text='GaAs', variable=var3, onvalue=1, offvalue=0, command=selection).grid(row=10, column=3,
                                                                                                  stick='w')
tk.Checkbutton(window, text='custom', variable=var4, onvalue=1, offvalue=0, command=selection).grid(row=11, column=3,
                                                                                                    stick='w')
window.grid_columnconfigure(0, minsize=30)
window.grid_rowconfigure(10, minsize=30)
window.grid_rowconfigure(11, minsize=30)

tk.Button(window, text='START', bg='peachpuff', command=start).grid(row=len(text) + 2, column=1, columnspan=3,
                                                                    stick='we')
tk.Button(window, text='CLEAR', bg='turquoise', command=clear).grid(row=len(text) + 3, column=1, columnspan=3,
                                                                    stick='we')
tk.Button(window, text='HELP', bg='palegreen', command=help).grid(row=len(text) + 4, column=1, columnspan=3,
                                                                  stick='we')

fig = Figure(figsize=(8, 8))
ax = fig.add_subplot()
ax.grid()
ax.set_xlabel("x [cm]")
ax.set_ylabel("E [eV]")

canvas = FigureCanvasTkAgg(fig)
canvas.get_tk_widget().grid(row=0, column=5, columnspan=15, rowspan=15)
window.grid_columnconfigure(4, minsize=200)

window.mainloop()
