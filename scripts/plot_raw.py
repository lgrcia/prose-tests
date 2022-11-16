import matplotlib.pyplot as plt
from prose import Observation

files = snakemake.input

for i, f in enumerate(files):
    obs = Observation(f)
    v = f.split("/")[-1][-10:-1]
    plt.plot(obs.time, obs.raw_flux, label=f"...{v}")

plt.legend()
plt.title(snakemake.wildcards.dataset)
plt.savefig(snakemake.output[0])