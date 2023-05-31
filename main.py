import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.integrate import odeint

def movement(y, t, M, m, c, g, R):
    #y[0, 1, 2, 3] = phi, psi, phi', psi'
    #dy[] = dt
    dy = np.zeros_like(y)
    dy[0] = y[2]
    dy[1] = y[3]
    a11 = 1
    a12 = np.cos(y[0])
    b1 = -1*2*(c / m) * (1 - np.cos((y[0]+y[1])/2)) * np.sin((y[0]+y[1])/2) - (g / R) * np.sin(y[0])
    a21 = np.cos(y[0])
    a22 = 1 + 2*(M / m)
    b2 = -1*2*(c / m) * (1 - np.cos((y[0]+y[1])/2)) * np.sin((y[0]+y[1])/2) + y[2]*y[2] * np.sin(y[0])
    dy[2] = (b1 * a22 - b2 * a12)/(a11 * a22 - a21 * a12)
    dy[3] = (a11 * b2 - a21 * b1)/(a11 * a22 - a21 * a12)
    return dy

M = 10
m = 100
c = 1
g = 9.81
R = 1

def rotate(origin, point, angle):
    x = origin[0]
    y = origin[1]
    px = point[0]
    py = point[1]
    nx = x + np.cos(angle) * (px - x) - np.sin(angle) * (py - y)
    ny = y + np.sin(angle) * (px - x) + np.cos(angle) * (py - y)
    return np.array([nx, ny])

def distance(x1, x2, y1, y2):
    return ((x1 - x2)**2 + (y1 - y2)**2)**0.5

t = np.linspace(0, 10, 1001)

r = 0.8 #internal radius
#R = 1 #external radius



###
phi0 = np.pi/3  # 1 3pi/2
psi0 = 0.6  # 1
dphi0 = 0
dpsi0 = 0
Y0 = [phi0, psi0, dphi0, dpsi0]

t_fin = 10
Nt = 1001
t2 = np.linspace(0, t_fin, Nt)
Y = odeint(movement, Y0, t2, (M, m, c, g, R))

phi = Y[:, 0]
psi = Y[:, 1]
dphi = Y[:, 2]
dpsi = Y[:, 3]
ddphi = [movement(y, t, M, m, c, g, R)[2] for y, t in zip(Y, t)]
ddpsi = [movement(y, t, M, m, c, g, R)[3] for y, t in zip(Y, t)]
###

#speed = 2 #rotation speed
x0 = 2 - phi * R #origin C position
y0 = 2
#xA, yA = rotate([x0, y0], [x0, y0 + R], speed*t)
#xB, yB = rotate([x0, y0], [x0+(r + (R - r)/2)/(2**0.5), y0 - (r + (R - r)/2)/(2**0.5)], speed*(np.sin(t)))
xA = []
xB = []
yA = []
yB = []
for i in range(Nt):
    xA.append(x0[i] + np.sin(psi[i]))
    yA.append(y0 + np.cos(psi[i]))
    xB.append(x0[i] + np.sin(np.pi - phi[i]))
    yB.append(y0 + np.cos(np.pi - phi[i]))

#making N
N = []
for i in range(Nt):
    N.append(m * (g * np.cos(phi[i]) + R * (dphi[i]*dphi[i] - ddpsi[i] * np.sin(phi[i]))) + 2 * R * c * (1 - np.cos((phi[i]+psi[i])/2))*np.cos((phi[i]+psi[i])/2))

fig_graphs = plt.figure(figsize=[13, 7])
ax_graphs = fig_graphs.add_subplot(2, 2, 1)
ax_graphs.plot(t, phi, color='blue')
ax_graphs.set_title("phi(t)")
ax_graphs.set(xlim = [0, t_fin])
ax_graphs.grid(True)

ax_graphs = fig_graphs.add_subplot(2, 2, 2)
ax_graphs.plot(t, psi, color='red')
ax_graphs.set_title("psi(t)")
ax_graphs.set(xlim = [0, t_fin])
ax_graphs.grid(True)

ax_graphs = fig_graphs.add_subplot(2, 2, 3)
ax_graphs.plot(t, dphi, color='orange')
ax_graphs.set_title("phi'(t)")
ax_graphs.set(xlim = [0, t_fin])
ax_graphs.grid(True)

ax_graphs = fig_graphs.add_subplot(2, 2, 4)
ax_graphs.plot(t, N, color='green')
ax_graphs.set_title("N (t)")
ax_graphs.set(xlim = [0, t_fin])
ax_graphs.grid(True)


#making external circle
Alpha = np.linspace(0, 7, 50)
xCircBig = R * np.sin(Alpha)
yCircBig = R * np.cos(Alpha)
#making internal circle
xCircSmall = r * np.sin(Alpha)
yCircSmall = r * np.cos(Alpha)
#making ground
xGround = [0,0,10]
yGround= [1,1,1]

#making spring pattern
n = 13
h = 0.1
#xSpring = np.linspace(0, 1, 2*n + 1)
#ySpring = np.zeros(2*n + 1)
#ss = 0
#for i in range(2*n + 1):
#    ySpring[i] = h * np.sin(ss)
#    ss += 3.14/2

xSpringARR = [] #2d array of all X positions of spring
ySpringARR = []
for i in range(Nt):
    spX = np.linspace(0, distance(xA[i], xB[i], yA[i], yB[i]), 2*n+1)  # correct length spring
    spY = np.zeros(2*n+1)
    ss = 0
    for j in range(2*n+1):
        spY[j] = h*np.sin(ss)
        ss += 3.14 / 2
    # now we have full not rotated, but correct length spring
    for k in range(2*n+1):
        if yB[i] > yA[i] and xB[i] < xA[i]:
            spX[k], spY[k] = rotate([0, 0], [spX[k], spY[k]], np.pi - np.math.atan(abs(yB[i]-yA[i])/abs(xA[i]-xB[i])))
        elif yB[i] > yA[i] and xB[i] > xA[i]: #=
            spX[k], spY[k] = rotate([0, 0], [spX[k], spY[k]], np.math.atan(abs(yB[i] - yA[i]) / abs(xB[i] - xA[i])))
        elif yB[i] < yA[i] and xB[i] > xA[i]:
            spX[k], spY[k] = rotate([0, 0], [spX[k], spY[k]], -np.math.atan(abs(yA[i] - yB[i]) / abs(xB[i] - xA[i])))
        elif yB[i] < yA[i] and xB[i] < xA[i]: #=
            spX[k], spY[k] = rotate([0, 0], [spX[k], spY[k]], np.pi + np.math.atan(abs(yA[i] - yB[i]) / abs(xA[i] - xB[i])))
        elif yB[i] == yA[i]:
            if xA[i] <= xB[i]:
                spX[k], spY[k] = rotate([0, 0], [spX[k], spY[k]], 0)
            else:
                spX[k], spY[k] = rotate([0, 0], [spX[k], spY[k]], np.pi)
        elif xB[i] == xA[i]:
            if yB[i] >= yA[i]:
                spX[k], spY[k] = rotate([0, 0], [spX[k], spY[k]], np.pi/2)
            else:
                spX[k], spY[k] = rotate([0, 0], [spX[k], spY[k]], -np.pi/2)
    # now we have rotated spring. last step - put it in correct location
    for k in range(2*n+1):
        spX[k] += xA[i]
        spY[k] += yA[i]
    # adding
    xSpringARR.append(spX)
    ySpringARR.append(spY)

#generate string
for i in range(Nt):
    orientsX = np.linspace(xA[i], xB[i], 2*n+1)
    orientsY = np.linspace(yA[i], yB[i], 2*n+1)


fig = plt.figure(figsize=[9, 9])
ax = fig.add_subplot(1, 1, 1)
ax.axis('equal')
ax.set(xlim=[0, 4], ylim=[0,4])

xaxis = np.array([5, 0])
spring = ax.plot(xSpringARR[0], ySpringARR[0], color=[0,0,0])[0]

A = ax.plot(xA[0], yA[0], 'o', color=[1, 0, 0])[0]
B = ax.plot(xB[0], yB[0], 'o', color=[0, 0, 1], markersize=25)[0]
circeBig = ax.plot(xCircBig+x0[0], yCircBig+y0, color=[0,0,0])[0]
circeSmall = ax.plot(xCircSmall+x0[0], yCircSmall+y0, color=[0,0,0])[0]
ax.plot(xGround, yGround, color='black', linewidth=3)

def kadr(i):
    A.set_data(xA[i], yA[i])
    B.set_data(xB[i], yB[i])
    circeBig.set_data(xCircBig+x0[i], yCircBig+y0)
    circeSmall.set_data(xCircSmall+x0[i], yCircSmall+y0)
    spring.set_data(xSpringARR[i], ySpringARR[i])
    return [A, B, circeBig, circeSmall, spring]

kino = FuncAnimation(fig, kadr, interval=t[1]-t[0], frames=len(t))

plt.show()
