rule sampling:
    input:
        "src/scripts/sampling.py"
    output:
        "src/data/Fig3.npy"
    conda:
        "environment.yml"
    script:
        "src/scripts/sampling.py"