import os
from scipy.spatial import KDTree
import numpy as np
from tqdm import tqdm
from joblib import Parallel, delayed
import argparse
import sys


class TADAAnalyzer:
	"""
	TADA (Topological Analysis of Defect Alignments) analyzer.
	
	This class encapsulates all parameters and methods for TADA analysis,
	with radius automatically calculated from long_axis.
	"""
	
	def __init__(self, long_axis, file_path, output_dir='tada_output', n_cores=-1, batch_size=100):
		"""
		Initialize TADA analyzer.
		
		Args:
			long_axis (float): Long axis length of uniaxial mesogen
			file_path (str): Path to input file containing positions and directors
			output_dir (str): Directory to save results
			n_cores (int): Number of CPU cores to use (-1 for all available)
			batch_size (int): Number of realizations per batch
		"""
		self.long_axis = long_axis
		self.file_path = file_path
		self.output_dir = output_dir
		self.n_cores = n_cores
		self.batch_size = batch_size
		
		# Calculate derived parameters
		self.radius = 1.5 * self.long_axis  # Search radius = 1.5 * long_axis
		
		# Data storage
		self.positions = None
		self.directors = None
		self.neighbors_list = None
		
		# Validation
		self._validate_parameters()
		
	def _validate_parameters(self):
		"""Validate input parameters."""
		if self.long_axis <= 0:
			raise ValueError("long_axis must be a positive number")
		if not os.path.exists(self.file_path):
			raise FileNotFoundError(f"File '{self.file_path}' does not exist")
		if self.n_cores == 0:
			raise ValueError("n_cores must be non-zero")
		if self.batch_size <= 0:
			raise ValueError("batch_size must be a positive integer")
	
	def load_data(self):
		"""
		Load positions and directors from file.
		
		Expected file format: x y z nx ny nz (space-separated)
		"""
		print(f"Reading data from {self.file_path}...")
		with open(self.file_path, 'r') as file:
			data = file.readlines()
		
		data = np.array([list(map(float, dat.strip().split())) for dat in data])
		self.positions = data[:, :3]
		self.directors = data[:, 3:]
		print(f"Loaded {len(self.positions)} particles")
		
	def find_neighbors(self):
		"""Find and store neighbors within the search radius."""
		print(f"Finding neighbors within radius {self.radius}...")
		self.neighbors_list = self.neighbors(self.positions, self.radius)
		print(f"Found neighbor list for {len(self.positions)} particles")
		
	def print_parameters(self):
		"""Print all parameters used in the analysis."""
		print("\n" + "="*50)
		print("TADA ANALYSIS PARAMETERS")
		print("="*50)
		print(f"Long axis length: {self.long_axis}")
		print(f"Search radius: {self.radius} (1.5 × long_axis)")
		print(f"Input file: {self.file_path}")
		print(f"Output directory: {self.output_dir}")
		print(f"CPU cores: {self.n_cores if self.n_cores > 0 else 'all available'}")
		print(f"Batch size: {self.batch_size}")
		print(f"Number of particles: {len(self.positions) if self.positions is not None else 'Not loaded'}")
		print("="*50)
		
	def neighbors(self, positions, radius):
		"""Create sorted neighbor list with padding - compatible with older Python versions"""
		# Make KDTree for positions
		KD = KDTree(positions)
		
		# Find neighbors for each position (compatible with older scipy versions)
		neighbors_list = KD.query_ball_point(positions, radius)
		
		# Find maximum number of neighbors at each coordinate
		max_neighbors = max(len(neighbors) for neighbors in neighbors_list)
		
		# Pad neighbor list
		padded_neighbors = np.full((len(neighbors_list), max_neighbors), -1)  # Using -1 as a placeholder
		for i, indices in enumerate(neighbors_list):
			padded_neighbors[i, :len(indices)] = indices
		
		return padded_neighbors

	def growth(self, LAM, count, reference, directors, neighbors, center, realizations):
		"""Growth algorithm for parallel TADA"""
		# Fetch LAM for neighbors of centers
		LAM0 = np.where(neighbors[center] == -1, -1, LAM[np.arange(realizations)[:, None], neighbors[center]])

		# Loop over all mesogens in neighborhood, that haven't been assigned a vector, for each center
		while np.sum(np.sum(LAM0 == 0, axis=1)) != 0:
			
			# Check how many centers' neighbors are unassigned
			remaining = np.where(np.sum(LAM0 == 0, axis=1) != 0)[0]
			
			# Which mesogen is closest to center and unassigned
			mask_0 = (LAM0 == 0)
			indices_0 = np.where(np.any(mask_0, axis=1), np.argmax(mask_0, axis=1), -1)

			# Denote closest unassigned mesogen to central mesogen as Mesogen Under Consideration for Assignment (MUCA)
			muca = (neighbors[center, indices_0])
			
			# Fetch LAM for neighbors of muca
			LAM1 = np.where(neighbors[muca] == -1, -1, LAM[np.arange(realizations)[:, None], neighbors[muca]])

			# Which mesogen is closest to muca and assigned
			mask_1 = (LAM1 > 0)
			indices_1 = np.where(np.any(mask_1, axis=1), np.argmax(mask_1, axis=1), -1)
			
			# Find closest assigned mesogen to MUCA
			closest = neighbors[muca, indices_1]
			
			# Calculate dot products
			directions = (np.einsum('ij,ij->i', reference[np.arange(realizations), closest], directors[np.arange(realizations), muca]))
			
			# Assign vector to reference field
			indices = np.where(directions >= 0)
			reference[np.arange(realizations)[indices], muca[indices]] = directors[np.arange(realizations)[indices], muca[indices]]
			indices = np.where(directions < 0)
			reference[np.arange(realizations)[indices], muca[indices]] = (-1) * directors[np.arange(realizations)[indices], muca[indices]]
				
			# update number of assigned mesogens 
			count[remaining] += 1
			
			# update LAM
			rows = np.arange(realizations)[remaining]
			cols = muca[remaining]
			LAM[rows, cols] = count[remaining]

			# update LAM0
			LAM0 = np.where(neighbors[center] == -1, -1, LAM[np.arange(realizations)[:, None], neighbors[center]])
		
		return LAM, count, reference

	def propagate(self, positions, directors, realizations, start_indices):
		"""Parallel TADA algorithm - returns only modified director field"""
		# Maximum number of realizations are the number of mesogens
		realizations = realizations * (realizations < len(positions)) + len(positions) * (realizations >= len(positions))

		directors = np.repeat(directors[np.newaxis, :, :], realizations, axis=0)
			
		# Initialize number of mesogens assigned a vector
		count = np.zeros(realizations)
		
		# Initialize reference director field
		reference = np.zeros_like(directors)
		
		# Initialize LAM
		LAM = np.zeros([realizations, len(positions)])
		
		# Use provided start indices as centers
		center = start_indices
		
		reference[np.arange(realizations), center] = directors[np.arange(realizations), center]
		count += 1
		
		# Add mesogen to LAM
		current = 1
		LAM[np.arange(realizations), center] = current

		# Check if all realizations have been assigned a vector
		check = np.max(LAM, axis=1) == len(positions)
		
		# Loop till all mesogens have been assigned a vector
		while np.all(check) == 0:

			# Loop through neighbors of most recently added mesogen of LCM, and assign vectors accordingly
			LAM, count, reference = self.growth(LAM, count, reference, directors, self.neighbors_list, center, realizations)

			# Choose new central mesogen
			current += 1
			center = np.where(LAM == current)[1]
			
			# Check if all realizations have been assigned an appropriate vector field
			check = np.max(LAM, axis=1) == len(positions)
			
		# Return only the modified director field (reference), not positions
		return reference

	def batched(self, batch_start, batch_size, positions, directors, save_dir=None):
		"""Process a batch of realizations with TADA_parallel"""
		# Calculate how many realizations to process in this batch
		remaining_points = len(positions) - batch_start
		batch_realizations = min(batch_size, remaining_points)
		
		if batch_realizations <= 0:
			return None
		
		# Create start indices for this batch - ensure no overlap with other batches
		start_indices = np.arange(batch_start, batch_start + batch_realizations)
		
		# Validate that we don't exceed the number of positions
		if np.max(start_indices) >= len(positions):
			# Adjust batch size to fit within available positions
			max_valid_start = len(positions) - 1
			batch_realizations = max(0, max_valid_start - batch_start + 1)
			if batch_realizations <= 0:
				return None
			start_indices = np.arange(batch_start, batch_start + batch_realizations)
		
		# Use TADA_parallel with specific start indices and pre-computed neighbors
		results = self.propagate(positions, directors, batch_realizations, start_indices)
		
		# Save results if save_dir is provided
		if save_dir is not None:
			# Create save directory if it doesn't exist
			os.makedirs(save_dir, exist_ok=True)
			
			# Create filename indicating the starting points used
			start_range = f"start_{batch_start:06d}_to_{batch_start + batch_realizations - 1:06d}"
			filename = f"director_field_{start_range}.npz"
			filepath = os.path.join(save_dir, filename)
			
			# Save only essential information
			np.savez_compressed(
				filepath,
				director_field=results,
				start_indices=start_indices,
				start_seeds=start_indices  # The start seeds are the same as start indices for systematic
			)
			
			print(f"Saved batch results to: {filepath}")
			
			# Delete results to free memory
			del results
			return None
		
		return results

	def tada(self, positions, directors, n_cores=-1, batch_size=100, save_dir="director_fields"):
		"""
		Analyze defect alignments using parallel TADA algorithm
		
		Args:
			positions: (N, 3) array of atom positions
			directors: (N, 3) array of directors
			n_cores: number of CPU cores to use (-1 for all available cores)
			batch_size: number of realizations per batch
			save_dir: directory to save results
			
		Returns:
			None (results are saved to files)
		"""
		N = len(positions)
		realizations = N  # Always use all points
		
		# Calculate number of batches needed
		num_batches = (realizations + batch_size - 1) // batch_size
		
		# Create batch starting indices - ensure no overlap
		batch_starts = [i * batch_size for i in range(num_batches)]
		
		# Validate that all batch starts are within bounds
		valid_batch_starts = [start for start in batch_starts if start < N]
		
		if not valid_batch_starts:
			raise ValueError("No valid batch starts found. Check batch_size and number of positions.")
		
		# Process batches in parallel using joblib and save results
		Parallel(n_jobs=n_cores)(
			delayed(self.batched)(batch_start, batch_size, positions, directors, save_dir) 
			for batch_start in tqdm(valid_batch_starts, desc="Processing and saving batches")
		)
		
		# Create a summary file with batch information
		summary_file = os.path.join(save_dir, "batch_summary.txt")
		with open(summary_file, 'w') as f:
			f.write(f"TADA Parallel Results Summary\n")
			f.write(f"============================\n")
			f.write(f"Total particles: {N}\n")
			f.write(f"Total realizations: {realizations}\n")
			f.write(f"Batch size: {batch_size}\n")
			f.write(f"Number of batches: {len(valid_batch_starts)}\n")
			f.write(f"Radius: {self.radius}\n")
			f.write(f"Number of cores used: {n_cores if n_cores > 0 else 'all available'}\n")
			f.write(f"\nBatch information:\n")
			for i, batch_start in enumerate(valid_batch_starts):
				batch_end = min(batch_start + batch_size - 1, N - 1)
				f.write(f"Batch {i+1}: points {batch_start}-{batch_end}\n")
		
		print(f"All batch results saved to: {save_dir}")
		print(f"Summary file created: {summary_file}")
		return None

	def run_analysis(self):
		"""
		Run complete TADA analysis workflow.
		
		Returns:
			None (results are saved to files)
		"""
		try:
			self.print_parameters()
			
			# Step 1: Load data
			self.load_data()
			
			# Step 2: Find neighbors (once, outside parallel loop)
			self.find_neighbors()
			
			# Step 3: Run TADA analysis
			print(f"Running TADA analysis with {len(self.positions)} realizations...")
			self.tada(
				positions=self.positions,
				directors=self.directors,
				n_cores=self.n_cores,
				batch_size=self.batch_size,
				save_dir=self.output_dir
			)
			
			print("Analysis complete!")
			
		except Exception as e:
			print(f"Error during analysis: {e}")
			raise


def main():
	"""Main function with CLI interface for TADA analysis."""
	parser = argparse.ArgumentParser(description='TADA (Topological Analysis of Defect Alignments)')
	parser.add_argument('--long-axis', type=float, required=True,
	                   help='Long axis length of uniaxial mesogen (radius = 1.5 * long_axis)')
	parser.add_argument('--file', type=str, required=True,
	                   help='Input file containing positions and directors (format: x y z dx dy dz)')
	parser.add_argument('--output', type=str, default='tada_output',
	                   help='Output directory for results (default: tada_output)')
	parser.add_argument('--cores', type=int, default=-1,
	                   help='Number of CPU cores to use (-1 for all available, default: -1)')
	parser.add_argument('--batch-size', type=int, default=100,
	                   help='Number of realizations per batch (default: 100)')
	
	args = parser.parse_args()
	
	try:
		# Create TADA analyzer
		analyzer = TADAAnalyzer(
			long_axis=args.long_axis,
			file_path=args.file,
			output_dir=args.output,
			n_cores=args.cores,
			batch_size=args.batch_size
		)
		
		# Run analysis
		analyzer.run_analysis()
		
	except Exception as e:
		print(f"Error: {e}")
		sys.exit(1)

if __name__ == "__main__":
	main() 
