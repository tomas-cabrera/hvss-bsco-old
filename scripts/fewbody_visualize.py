import numpy as np
import polars as pl
import matplotlib.pyplot as plt

###############################################################################

###############################################################################

###############################################################################

traj = pl.read_csv("trajectory0.dat")
traj = traj[1:-1]
traj.sort("t")
print(traj)

plt.plot(np.arange(traj.height), traj[:, "t"].to_numpy())
plt.show()
plt.close()

fig, ax = plt.subplots(1, 3, figsize=(18, 6))
colors = ("r","g","b")
cmaps = ("Reds_r", "Greens_r", "Blues_r")
cp = (
    ("x","y"),
    ("y","z"),
    ("x","z"),
)
zoom = 20
for i in range(3):
    for j in range(3):
        ax[i].plot(
            traj[:, "{}{}".format(cp[i][0], j)].to_numpy(),
            traj[:, "{}{}".format(cp[i][1], j)].to_numpy(),
            color=colors[j],
            label="Object {}".format(j),
        )
        ax[i].scatter(
            traj[0, "{}{}".format(cp[i][0], j)],
            traj[0, "{}{}".format(cp[i][1], j)],
            color=colors[j])
        ax[i].scatter(
            traj[-1, "{}{}".format(cp[i][0], j)],
            traj[-1, "{}{}".format(cp[i][1], j)],
            color=colors[j],
            marker="*",
        )
    ax[i].set_xlim((-zoom, zoom))
    ax[i].set_ylim((-zoom, zoom))
    ax[i].set_xlabel(cp[i][0])
    ax[i].set_ylabel(cp[i][1])
    ax[i].legend()
plt.show()
plt.close()

fig = plt.figure()
ax = fig.add_subplot(projection="3d")
for j in range(3):
    ax.plot(
        traj[:, "{}{}".format("x", j)].to_numpy()[:,0],
        traj[:, "{}{}".format("y", j)].to_numpy()[:,0],
        traj[:, "{}{}".format("z", j)].to_numpy()[:,0],
        label="Object {}".format(j),
        color=colors[j],
        alpha=0.5,
    )
    ax.scatter(
        traj[0, "{}{}".format("x", j)],
        traj[0, "{}{}".format("y", j)],
        traj[0, "{}{}".format("z", j)],
        color=colors[j])
    ax.scatter(
        traj[-1, "{}{}".format("x", j)],
        traj[-1, "{}{}".format("y", j)],
        traj[-1, "{}{}".format("z", j)],
        color=colors[j],
        marker="*",
    )
ax.set_xlim((-zoom, zoom))
ax.set_ylim((-zoom, zoom))
ax.set_zlim((-zoom, zoom))
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.legend()
plt.show()
plt.close()

