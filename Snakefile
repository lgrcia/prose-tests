from pathlib import Path
from glob import glob
import yaml

configfile: "config.yaml"

datasets = [Path(p).stem for p in glob("/Users/lgrcia/code/prose-tests/datasets/*")]

versions = [
    "f2495763976d692b5cf4a86592d9f9489e56ef7f", # master
    "31aff30dbbf2d6cf9e340a70acc76c23b283d49b", # 2.3.0
    ]

for version in versions:
    base = yaml.load(open("envs/prose_base.yaml", "r"), Loader=yaml.FullLoader)
    prose = base["dependencies"][-1]["pip"][0]
    base["dependencies"][-1]["pip"][0] = prose.replace("version", str(version))
    yaml.dump(base, open(f"envs/prose_{version}.yaml", "w"))

rule all:
    input:
        expand("data/products/{dataset}/{version}.phot", version=versions, dataset=datasets),
        expand("figures/{dataset}/diff_photometry.png", dataset=datasets),
        expand("figures/{dataset}/raw_photometry.png", dataset=datasets)


rule reduce_dataset:
    input: 
        "datasets/{dataset}",
        "envs/prose_{version}.yaml"
    output: 
        temp("data/products/{dataset}/{version}_nodiff.phot"),
        "data/products/{dataset}/{version}_processing_time.yaml"
    conda: "envs/prose_{version}.yaml" 
    script: "scripts/reduction.py"

rule make_diff:
    input: 
        "data/products/{dataset}/{version}_nodiff.phot",
        "datasets/{dataset}/target.yaml"
    output: "data/products/{dataset}/{version}.phot"
    conda: "envs/prose_{version}.yaml" 
    script: "scripts/make_diff.py"

rule plot_diff_photometry:
    input: expand("data/products/{{dataset}}/{version}.phot", version=versions)
    output: "figures/{dataset}/diff_photometry.png"
    conda: "envs/prose_master.yaml" 
    script: "scripts/plot_diff.py"
    
rule plot_raw_photometry:
    input: expand("data/products/{{dataset}}/{version}.phot", version=versions)
    output: "figures/{dataset}/raw_photometry.png"
    conda: "envs/prose_master.yaml" 
    script: "scripts/plot_raw.py"
    

