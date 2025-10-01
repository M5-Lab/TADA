# TADA - Topological Analysis of Defect Alignments

A Python implementation of the TADA (Topological Analysis of Defect Alignments) algorithm for analyzing defect alignments in nematic liquid crystal systems using parallel processing and batch optimization.

## Overview

TADA is a sophisticated algorithm designed to analyze the topological structure of defect alignments in nematic liquid crystal systems. It uses a "growth" algorithm that propagates director field assignments from seed particles outward, creating modified director fields that reveal the underlying topological structure.

## Features

- **Parallel Processing**: Utilizes multiple CPU cores for efficient computation
- **Batch Processing**: Processes large datasets in manageable batches
- **Optimized Neighbor Finding**: Pre-computes neighbor lists once for all batches
- **Modular Design**: Clean class-based architecture for easy extension
- **CLI Interface**: Simple command-line interface for batch processing
- **Comprehensive Output**: Saves modified director fields and analysis summaries

## Installation

### Prerequisites

```bash
pip install numpy scipy tqdm joblib
```

### Required Python Packages

- `numpy`: Numerical computations
- `scipy`: Spatial data structures (KDTree)
- `tqdm`: Progress bars
- `joblib`: Parallel processing

## Usage

### Command Line Interface

```bash
# Basic usage
python3 TADA.py --long-axis 3.0 --file data.txt

# Advanced usage with custom settings
python3 TADA.py --long-axis 4.5 --file geometry.txt --output results --cores 8 --batch-size 200
```

#### CLI Arguments

- `--long-axis` (required): Long axis length of uniaxial mesogen (radius = 1.5 × long_axis)
- `--file` (required): Input file containing positions and directors
- `--output` (optional): Output directory for results (default: 'tada_output')
- `--cores` (optional): Number of CPU cores to use (default: -1 for all available)
- `--batch-size` (optional): Number of realizations per batch (default: 100)

### Programmatic Usage

```python
from TADA import TADAAnalyzer

# Create analyzer
analyzer = TADAAnalyzer(
    long_axis=3.0,
    file_path='data.txt',
    output_dir='results',
    n_cores=8,
    batch_size=100
)

# Run analysis
analyzer.run_analysis()
```

## Input Format

The input file should contain space-separated values with the following format:

```
x y z nx ny nz
x y z nx ny nz
...
```

Where:
- `x y z`: Particle positions (3D coordinates)
- `nx ny nz`: Director vectors (3D unit vectors)

### Example Input File

```
0.0 0.0 0.0 1.0 0.0 0.0
1.0 0.0 0.0 0.0 1.0 0.0
0.0 1.0 0.0 0.0 0.0 1.0
...
```

## Output Files

### Batch Results

TADA saves results in compressed `.npz` files for each batch:

- `director_field_start_XXXXXX_to_XXXXXX.npz`: Contains:
  - `director_field`: Modified director fields for the batch
  - `start_indices`: Starting particle indices used as seeds
  - `start_seeds`: Same as start_indices (for systematic processing)

### Summary File

- `batch_summary.txt`: Contains analysis parameters and batch information

### Example Output Structure

```
tada_output/
├── director_field_start_000000_to_000099.npz
├── director_field_start_000100_to_000199.npz
├── director_field_start_000200_to_000299.npz
└── batch_summary.txt
```

## Algorithm Details

### Core TADA Algorithm

1. **Neighbor Finding**: Uses KDTree to find all particles within the search radius
2. **Growth Algorithm**: Propagates director assignments from seed particles
3. **Parallel Processing**: Each batch processes multiple seed particles simultaneously
4. **Director Field Modification**: Creates modified director fields based on topological analysis

### Key Parameters

- **Search Radius**: `1.5 × long_axis` (automatically calculated)
- **Batch Size**: Number of realizations processed per parallel batch
- **Realizations**: Equal to the number of particles (each particle as a seed)

### Optimization Features

- **Pre-computed Neighbors**: Neighbor lists calculated once and reused across all batches
- **Memory Management**: Results saved immediately to free memory
- **Parallel Batching**: Multiple CPU cores process different batches simultaneously

## Class Reference

### TADAAnalyzer

Main class for TADA analysis.

#### Constructor

```python
TADAAnalyzer(long_axis, file_path, output_dir='tada_output', n_cores=-1, batch_size=100)
```

#### Methods

- `load_data()`: Load positions and directors from file
- `find_neighbors()`: Pre-compute neighbor lists
- `run_analysis()`: Execute complete analysis workflow
- `print_parameters()`: Display analysis parameters

#### Internal Methods

- `neighbors()`: Create sorted neighbor list with padding
- `growth()`: Growth algorithm for director field propagation
- `propagate()`: Parallel TADA algorithm implementation
- `batched()`: Process a batch of realizations
- `tada()`: Main analysis method with parallel processing

## Examples

### Example 1: Basic Analysis

```bash
python3 TADA.py --long-axis 3.0 --file liquid_crystal_data.txt
```

### Example 2: High-Performance Analysis

```bash
python3 TADA.py \
    --long-axis 4.5 \
    --file large_system.txt \
    --output high_perf_results \
    --cores 16 \
    --batch-size 500
```

### Example 3: Programmatic Usage

```python
from TADA import TADAAnalyzer

# Analyze multiple systems
systems = ['system1.txt', 'system2.txt', 'system3.txt']
long_axes = [3.0, 4.5, 2.8]

for system, long_axis in zip(systems, long_axes):
    analyzer = TADAAnalyzer(
        long_axis=long_axis,
        file_path=system,
        output_dir=f'results_{system.split(".")[0]}',
        n_cores=8,
        batch_size=100
    )
    analyzer.run_analysis()
```

## Performance Considerations

### Memory Usage

- **Neighbor Lists**: Stored once per analysis (not per batch)
- **Batch Processing**: Results saved immediately to minimize memory usage
- **Parallel Processing**: Memory usage scales with number of cores

### Computational Complexity

- **Neighbor Finding**: O(N log N) - done once
- **Growth Algorithm**: O(N × average_neighbors) per realization
- **Parallel Speedup**: Near-linear scaling with available cores

### Recommended Settings

- **Small Systems** (< 1000 particles): `batch_size=50`, `cores=4`
- **Medium Systems** (1000-10000 particles): `batch_size=100`, `cores=8`
- **Large Systems** (> 10000 particles): `batch_size=200`, `cores=16`

## Troubleshooting

### Common Issues

1. **Memory Errors**: Reduce batch_size or number of cores
2. **Slow Performance**: Increase batch_size (if memory allows)
3. **File Not Found**: Check file path and format
4. **Invalid Parameters**: Ensure long_axis > 0

### Performance Tips

- Use SSD storage for better I/O performance
- Monitor memory usage with large datasets
- Adjust batch_size based on available RAM
- Use appropriate number of cores for your system

## License

MIT License

Copyright (c) 2025 TADA Implementation

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do to, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## References

- TADA algorithm for topological analysis of defect alignments
- Parallel processing optimization for large-scale nematic liquid crystal systems
- Efficient neighbor finding using spatial data structures (KDTree)