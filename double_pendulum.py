import matplotlib.pyplot as plt
import numpy as np
import pygame as pg


class IVProblem:
    def __init__(self, rhs, y0, interval, name='IVP'):
        """
        :param rhs: 'right hand side' function of the ordinary differential equation f(t,y)
        :param y0: array with initial values
        :param interval: start and end value of the interval of independent variables; often initial and end time
        :param name: descriptive name of the problem
        """
        self.rhs = rhs
        self.y0 = y0
        self.t0, self.tend = interval
        self.name = name


class IVPsolver:
    """
    IVP solver class for explicit one-step discretization methods
    with constant step size
    """

    def __init__(self, problem, discretization, stepsize):
        self.problem = problem
        self.discretization = discretization
        self.stepsize = stepsize
        self.ts = self.problem.t0  # for pygame animation
        self.ys = self.problem.y0  # for pygame animation
        self.t_list = [self.ts]  # for graph
        self.y_list = [self.ys]  # for graph

    def solve(self):
        while self.ts <= self.problem.tend:
            self.step()
        return self.t_list, self.y_list

    def step(self, append_to_list=True):
        self.ts, self.ys = self.discretization(self.problem.rhs, self.ts, self.ys,
                                               self.stepsize)
        if append_to_list:
            self.t_list.append(self.ts)
            self.y_list.append(self.ys)


def parameters_double_pendulum():
    l = 1
    g = 9.81
    return np.array([l, g])


def rhs(t, y):
    l, g = parameters_double_pendulum()
    z = np.zeros(4)
    z[0] = y[2]
    z[1] = y[3]
    M = np.array([[2, np.cos(y[1] - y[0])], [np.cos(y[1] - y[0]), 1]])
    v = np.array([-2 * g / l * np.sin(y[0]) + np.sin(y[1] - y[0]) * y[3] ** 2,
                  -g / l * np.sin(y[1]) - np.sin(y[1] - y[0]) * y[2] ** 2])
    z[2:4] = np.linalg.solve(M, v)
    return z


def expliciteuler(rhs, ts, ys, h):
    return ts + h, ys + h * rhs(ts, ys)


def rungekutta4(rhs, ts, ys, h):
    k1 = h * rhs(ts, ys)
    k2 = h * rhs(ts + h / 2., ys + k1 / 2.)
    k3 = h * rhs(ts + h / 2., ys + k2 / 2.)
    k4 = h * rhs(ts + h, ys + k3)
    return ts + h, ys + (k1 + 2 * k2 + 2 * k3 + k4) / 6.


def draw_pendulum(screen, p1, p2, m):
    pg.draw.line(screen, (255, 255, 255), p1, p2, 3)
    pg.draw.circle(screen, (255, 0, 0), p2, m)


def main():
    initial_angle1 = np.pi / 2
    initial_angle2 = 0
    double_pendulum = IVProblem(rhs, np.array([initial_angle1, initial_angle2, 0., 0.]), [0., 10.],
                                'mathem. pendulum')
    pendulum_Euler = IVPsolver(double_pendulum, expliciteuler, 0.001)
    pendulum_RK4 = IVPsolver(double_pendulum, rungekutta4, 0.001)
    tRK4, yRK4 = pendulum_RK4.solve()

    l, g = parameters_double_pendulum()

    # x-position-velocity diagrams
    pos, yRK4_transp = polar_to_cartesian(l, yRK4)
    dx1dt = l * yRK4_transp[2] * np.cos(yRK4_transp[0])
    dx2dt = dx1dt + l * yRK4_transp[3] * np.cos(yRK4_transp[1])
    plt.subplot(1, 2, 1), plt.plot(pos[0][0], dx1dt), \
    plt.title('Pendulum result with RK4'), \
    plt.xlabel('x_1'), plt.ylabel('dx_1/dt')
    plt.subplot(1, 2, 2), plt.plot(pos[1][0], dx2dt), \
    plt.title('Pendulum result with Runge Kutta'), \
    plt.xlabel('x_2'), plt.ylabel('dx_2/dt')

    plt.show()

    # animation double pendulum
    breite, höhe = 900, 600
    pg.init()
    screen = pg.display.set_mode([breite, höhe])
    screen2 = pg.Surface([breite, höhe])
    clock = pg.time.Clock()
    weitermachen = True
    first_run = True
    pendulum_RK4 = IVPsolver(double_pendulum, rungekutta4, 0.001)

    # inspired by https://github.com/Gravitar64/A-beautiful-code-in-Python/blob/master/Teil_12%20Doppelpendel.py
    while weitermachen:
        clock.tick(60)
        for event in pg.event.get():  # close animation
            if event.type == pg.QUIT or (event.type == pg.KEYDOWN and event.key == pg.K_ESCAPE):
                weitermachen = False
        screen.fill((0, 0, 0))

        screen.blit(screen2, (0, 0))

        for i in range(20):
            pendulum_RK4.step(False)

        if not first_run:
            pos_old = pos
        pos, _ = polar_to_cartesian(l, pendulum_RK4.ys[:2])
        pos = np.array(pos)

        scale = np.array([100, -100])
        translate = np.array([450, 200])
        draw_pendulum(screen, scale * pos[0] + translate, scale * pos[1] + translate, 10)  # pendulum 2
        draw_pendulum(screen, np.array([0, 0]) + translate, scale * pos[0] + translate, 10)  # pendulum 1
        if not first_run:
            pg.draw.line(screen2, (0, 255, 0), scale * pos[1] + translate, scale * pos_old[1] + translate, 2)
        first_run = False
        pg.display.flip()  # 2 Bildschirme, damit das Bild nicht flackert
    pg.quit()


def polar_to_cartesian(l, yRK4):
    yRK4_transp = np.transpose(yRK4)
    x1 = l * np.sin(yRK4_transp[0])
    y1 = -l * np.cos(yRK4_transp[0])
    x2 = x1 + l * np.sin(yRK4_transp[1])
    y2 = y1 - l * np.cos(yRK4_transp[1])
    return [(x1, y1), (x2, y2)], yRK4_transp


if __name__ == '__main__':
    main()
