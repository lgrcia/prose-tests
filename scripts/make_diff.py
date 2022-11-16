from prose import Observation
import yaml

obs = Observation(snakemake.input[0])
target = yaml.load(open(snakemake.input[1], "r"), Loader=yaml.FullLoader)["target"]

obs.target = target
obs.broeg2005()
obs.save(snakemake.output[0])