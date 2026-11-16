# TADA — Topology-Accommodating Direction Assignment Algorithm

TADA resolves nematic head–tail ambiguity by growing an oriented director field through a distance-sorted neighbor graph. This repository provides the current standalone algorithm and the complete pipeline required to run TADA through [OVITO Nematics](https://github.com/Sapta98/Nematics).

## Paper

The original TADA paper is available in this repository: [TADA: The Topology-Accommodating Direction Assignment Algorithm for Liquid Crystals](paper/TADA.pdf).

Saptarshi Saha, Amit Acharya, and Gerald J. Wang, *Journal of Chemical Theory and Computation* **2023**, *19*, 324–332. [https://doi.org/10.1021/acs.jctc.2c00826](https://doi.org/10.1021/acs.jctc.2c00826)

## Standalone quick start

```bash
python -m pip install -r requirements.txt
python TADA.py path/to/director_field.txt --cutoff 4.5 --seed 0 --output output/tada.txt
```

Input is a whitespace-delimited `x y z nx ny nz` table. Output contains the aligned directors and `assignment_order`.

## Full TADA pipeline

```bash
ovitos -m pip install --user git+https://github.com/Sapta98/Nematics.git
ovitos pipelines/run_tada.py examples/Plushalf output/plushalf_tada.dump --cutoff 4.5
```

The pipeline:

1. imports a LAMMPS dump;
2. converts raw scalar-first quaternion properties to nematic directors;
3. runs TADA;
4. constructs the antagonistic-pair discontinuity layer;
5. runs Discrete Layer Analysis (DLA); and
6. exports the aligned field and defect properties.

Pass `--no-dla` to stop after layer construction. The OVITO desktop GUI exposes the same workflow through **Add modification → Nematic Defect Analysis → TADA**.

## Representative examples

- `examples/Plushalf` — planar +1/2 defect.
- `examples/Comet` — three-dimensional comet configuration.

For additional 2D/3D examples and GUI documentation, see [OVITO Nematics](https://github.com/Sapta98/Nematics).

## License

MIT. See [LICENSE](LICENSE).
