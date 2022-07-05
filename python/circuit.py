import numpy as np
import matplotlib.pyplot as plt


def solver(T, sim_time):
    a_0 = 0.0
    b_0 = 0.0

    steps = int(sim_time / T)

    a = np.empty(steps)
    b = np.empty(steps)

    a_rk = np.empty(steps)
    b_rk = np.empty(steps)

    a[0] = a_rk[0] = a_0
    b[0] = b_rk[0] = b_0

    for step in range(steps - 1):

        time = T * step

        # Euler
        a[step + 1] = a[step] + T * a_dot(time, a[step], b[step])
        b[step + 1] = b[step] + T * b_dot(time, a[step], b[step])

        # Runge-Kutta
        a_k1 = a_dot(time, a_rk[step], b_rk[step])
        b_k1 = b_dot(time, a_rk[step], b_rk[step])
        a_k2 = a_dot(time + (T/2.0), a_rk[step]+T*(a_k1/2.0), b_rk[step]+T*(b_k1/2.0))
        b_k2 = b_dot(time + (T/2.0), a_rk[step]+T*(a_k1/2.0), b_rk[step]+T*(b_k1/2.0))
        a_k3 = a_dot(time + (T/2.0), a_rk[step]+T*(a_k2/2.0), b_rk[step]+T*(b_k2/2.0))
        b_k3 = b_dot(time + (T/2.0), a_rk[step]+T*(a_k2/2.0), b_rk[step]+T*(b_k2/2.0))
        a_k4 = a_dot(time + T, a_rk[step]+T*a_k3, b_rk[step]+T*b_k3)
        b_k4 = b_dot(time + T, a_rk[step]+T*a_k3, b_rk[step]+T*b_k3)

        a_rk[step + 1] = a_rk[step] + (1.0/6.0) * T * (a_k1 + 2.0 * a_k2 + 2.0 * a_k3 + a_k4)
        b_rk[step + 1] = b_rk[step] + (1.0/6.0) * T * (b_k1 + 2.0 * b_k2 + 2.0 * b_k3 + b_k4)

        print(f"{100*(step)/float(steps-1):.2f}%")

    return a_rk, b_rk, a, b


def a_dot(t, a, b):
    ## Equation
    # v = i*R + 1/C *1/s * i + L * s * i

    v = 15.0 #volts
    r = 10.0 #ohms
    l = 0.15 #henry
    c = 50e-6 # farads
    return b


def b_dot(t, a, b):
    ## Equation
    # v = i*R + 1/C *1/s * i + L * s * i

    v = 15.0 #volts
    r = 10.0 #ohms
    l = 0.15 #henry
    c = 50e-6 # farads
    return -(r/l)*b-(1/(l*c))*a+(1/l)*v#*np.sin(2*np.pi*50*t)


def main():
    T = 1e-4
    sim_time = 1.0 #s

    a_rk, b_rk, a, b = solver(T, sim_time)

    # plt.plot(b)
    plt.plot(b_rk)
    plt.grid()
    plt.show()


if __name__ == "__main__":
    main()