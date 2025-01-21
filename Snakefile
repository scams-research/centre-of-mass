rule sampling:
    input:
        "src/scripts/sampling.py"
    output:
        "src/data/Fig3.npy"
    conda:
        "environment.yml"
    cache: True
    script:
        "src/scripts/sampling.py"