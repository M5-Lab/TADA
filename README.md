# TADA — Traversal-based Alignment of Directors Algorithm

TADA resolves nematic head–tail ambiguity by growing an oriented director field through a distance-sorted neighbor graph. This repository provides the current standalone algorithm and the complete pipeline required to run TADA through [OVITO Nematics](https://github.com/Sapta98/Nematics).

## Standalone quick start

```bash
python -m pip install -r requirements.txt
python TADA.py examples/directors.txt --cutoff 4.5 --seed 0 --output output/tada.txt
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
- `examples/directors.txt` — compact standalone smoke-test input.

For additional 2D/3D examples and GUI documentation, see [OVITO Nematics](https://github.com/Sapta98/Nematics).

## License

MIT. See [LICENSE](LICENSE).
