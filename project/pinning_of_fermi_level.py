import tkinter as tk
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import calculations

matplotlib.use("TkAgg")

coef_phys_parameters = {'N_d0': 1e16, 'N_as': 1e17, 'E_out': 1e4}


def draw_figure(ax, canvas, results):
    ax.cla()
    ax.set_xlabel("x [cm]")
    ax.set_ylabel("E [eV]")
    ax.grid()
    if results['message'] == 'ok':
        ax.plot(results['x_s'], results['E_c_s'], c='red', label='Conduction Band')
        ax.plot(results['x_s'], results['E_f_s'], c='darkorange', label='Fermi Energy')
        ax.plot(results['x_s'][0:round(len(results['x_s'])/2)], results['E_as_s'][0:round(len(results['E_as_s'])/2)], c='green', label='Acceptor Energy')
        ax.plot(results['x_s'], results['E_d_s'], c='blue', label='Donor Energy')
        ax.plot(results['x_s'], results['E_v_s'], c='purple', label='Valence Band')
        ax.axhline(results['phi'], c='k', linestyle='dashed')
        ax.axvline(results['W'], c='k', linestyle='dashed')
        ax.legend(fontsize=10, loc='right')
    output_info(results)
    canvas.draw()
def get_calculated_values():
    args = {
        "E_gap": float(E_g.get()),  # Band gap [eV]
        "epsilon": float(epsilon.get()),  # Dielectric permittivity
        "m_h": float(m_h.get()),  # Effective hole mass [m_0]
        "m_e": float(m_e.get()),  # Effective electron mass [m_0]
        "E_d": float(E_d.get()),  # Donors level [eV]
        "N_d0": float(N_d0.get()) * coef_phys_parameters['N_d0'],  # Concentration of donors [cm^(-3)]
        "E_as": float(E_as.get()),  # Surface acceptors level [eV]
        "N_as": float(N_as.get()) * coef_phys_parameters['N_as'],  # Concentration of surface acceptors [cm^(-3)]
        "T": float(T.get()),  # Temperature [K]
        "E_out": float(E_out.get()) * coef_phys_parameters['E_out'],  # External electric field
    }
    return calculations.calculate(args)
def GaAs_calculated():
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
def Ge_calculated():
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
def Si_calculated():
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
def output_info(results):
    if results['message'] != 'ok':
        tk.Label(window, text=results['message'], fg='darkorange', bg='white', font="Arial 12", justify=tk.LEFT).grid(row=0, column=5, columnspan=15, rowspan=2,
                                                                    stick='we')
    else:
        tk.Label(window,
                 text=f"fermi level:\t\t{results['E_f_s'][0]:.6f} [eV]\nbending of zone PHI:\t\t{results['phi']:.6f} [eV]\nspace charge region W:\t{results['W']:.8f} [cm]",
                 fg='black', bg='white', font="Arial 12", justify=tk.LEFT).grid(row=0, column=5, columnspan=15, rowspan=2, stick='we')
def start():
    draw_figure(ax, canvas, get_calculated_values())
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
    help_text = """This program illustrates Fermi level pinning by surface acceptors \n\n\n1) Enter the required values in the forms on the left or click on the button: Si, Ge, GaAs, which will automatically set the values \n\n2) Press the START button to draw graphs \n\n3) Press the CLEAR button to clear all forms \n\n4) Press HELP for quick reference"""
    message_win = tk.Tk()
    message_win.config(bg='lightcyan')
    message_win.title('HELP')
    message_win.geometry(f"1150x300+300+200")
    tk.Label(message_win, text=help_text, fg='black', bg='lightcyan', font="Arial 14", justify=tk.LEFT).pack()
    message_win.mainloop()
##################################################################################
window = tk.Tk()
window.config(bg='white')
window.title('Pinning of Fermi Level')
window.geometry('{}x{}'.format(window.winfo_screenwidth(), window.winfo_screenheight()))
window.resizable(True, True)
text = [['forbidden zone E_g: ', '[eV]'], ['dielectric constant ε: ', ' '], ['effective hole masses m_h: ', 'm_0'],
        ['effective electron masses m_e: ', 'm_0'], ['donor level position E_d: ', '[eV]'],
        ['donor concentration N_d: ', '10^16[cm^(-3)]'],
        ['energy level position E_as: ', '[eV]'], ['surface acceptor density N_as: ', '10^17[cm^(-2)]'],
        ['temperature T: ', '[K]'], ['external electric field E_out:', '10^4 [V/m]']]
for i in range(len(text)):
    tk.Label(window, text=text[i][0], bg='white', font="Arial 10").grid(row=i, column=1, stick='w')
    tk.Label(window, text=text[i][1], bg='white', font="Arial 10").grid(row=i, column=3, stick='w')
    #window.grid_rowconfigure(i, minsize=20)
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
N_d0.insert(0, "10")
N_d0.grid(row=5, column=2)
E_as = tk.Entry(window)
E_as.insert(0, "2.5")
E_as.grid(row=6, column=2)
N_as = tk.Entry(window)
N_as.insert(0, "10")
N_as.grid(row=7, column=2)
T = tk.Entry(window)
T.insert(0, "300")
T.grid(row=8, column=2)
E_out = tk.Entry(window)
E_out.insert(0, "1")
E_out.grid(row=9, column=2)
tk.Label(window, text='semiconductor: ', bg='white').grid(row=10, column=1, stick='w')
tk.Button(window, text='Si', bg='white', command=Si_calculated, font="Arial 10").grid(row=10, column=2, columnspan=2, stick='we')
tk.Button(window, text='Ge', bg='white', command=Ge_calculated, font="Arial 10").grid(row=11, column=2, columnspan=2, stick='we')
tk.Button(window, text='GaAs', bg='white', command=GaAs_calculated, font="Arial 10").grid(row=11, column=1, stick='we')
window.grid_columnconfigure(0, minsize=30)
tk.Button(window, text='START', bg='peachpuff', command=start, font="Arial 10").grid(row=12, column=1, columnspan=3,
                                                                    stick='we')
tk.Button(window, text='CLEAR', bg='turquoise', command=clear, font="Arial 10").grid(row=13, column=1, columnspan=3,
                                                                    stick='we')
tk.Button(window, text='HELP', bg='palegreen', command=help, font="Arial 10").grid(row=14, column=1, columnspan=3,
                                                                  stick='we')
fig = Figure(figsize=(11, 8))
ax = fig.add_subplot()
ax.grid()
ax.set_xlabel("x [cm]")
ax.set_ylabel("E [eV]")
canvas = FigureCanvasTkAgg(fig)
canvas.get_tk_widget().grid(row=0, column=5, columnspan=15, rowspan=16)
window.grid_columnconfigure(4, minsize=20)
window.mainloop()
