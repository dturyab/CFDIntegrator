from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation,FFMpegWriter
import os

PATH_TO_FFMPEG = 'C:\\ffmpeg\\bin\\ffmpeg.exe'

dens = []
vel = []
pres = []
xs = []
ts = []

cur_path = os.path.dirname(__file__)
anim_relat_path = "res\\SODANIM_T0.1_N800.txt"
file_path = os.path.join(cur_path, anim_relat_path)

with open(anim_relat_path) as f:
    j = 1
    model = f.readline()
    alpha = float(model[model.find("ALPHA"):].split('=')[1])
    timestep = f.readline()
    while timestep:
        t, p = timestep.split(":")
        ts.append(float(t))

        dens_t, vel_t, pres_t = [], [], []
        p = p.split(';')
        for i in range(len(p)):
            x, rho, imp, E = list(map(float, p[i].split(',')))
            if j == 1:
                xs.append(x)
            dens_t.append(rho)
            vel_t.append(imp/rho)
            pres_t.append((E/rho-imp/rho/rho*imp/2.)*rho*(alpha-1))

        dens.append(dens_t)
        vel.append(vel_t)
        pres.append(pres_t)
        timestep = f.readline()
        j += 1


fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(5,7))
fig.subplots_adjust(hspace=0.25)
x1 = 1
ax1.set_xlim(0, x1)
ax1.set_ylabel("Density")
ax2.set_xlim(0, x1)
ax2.set_ylabel("Velocity")
ax3.set_xlim(0, x1)
ax3.set_ylabel("Pressure")
ax3.set_xlabel("X")

line1,  = ax1.plot([], [], lw=1)
line2,  = ax2.plot([], [], lw=1)
line3,  = ax3.plot([], [], lw=1)

def init():
    line1.set_data([], [])
    line2.set_data([], [])
    line3.set_data([], [])
    return line1, line2, line3

def animate(i):
    x = xs

    y = dens[i]
    ax1.set_ylim(min(y), 1.05*max(y))
    line1.set_data(x, y)

    y = vel[i]
    ax2.set_ylim(min(y), 1.05*max(y))
    line2.set_data(x, y)

    y = pres[i]
    ax3.set_ylim(min(y), 1.05*max(y))
    line3.set_data(x, y)

    return line1, line2, line3

anim = FuncAnimation(fig, animate, init_func=init,
                     frames=500, interval=30, blit=True)

plt.rcParams['animation.ffmpeg_path'] = PATH_TO_FFMPEG
FFwriter=FFMpegWriter(fps=25, extra_args=['-vcodec', 'libx264'])
anim.save('shocktubeAnim.mp4', writer=FFwriter)