import matplotlib.pyplot as plt
from prose import Observation, viz

files = snakemake.input

fig = plt.figure(None, (6*len(files), 4))

for i, f in enumerate(files):
    ax = plt.subplot(1, len(files), i+1)
    obs = Observation(f)
    obs.plot()
    if i == 0:
        ylim = ax.get_ylim()
    plt.title(f.split("/")[-1])
    noises = obs.noise_stats(verbose=False)
    viz.corner_text("\n".join([f"{n}: {v:.2e}" for n, v in noises.items()]))
    ax.set_ylim(ylim)

fig.suptitle(snakemake.wildcards.dataset)

plt.savefig(snakemake.output[0])