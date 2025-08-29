import numpy as np
import matplotlib.pyplot as plt
from pyapothesis.config import KMCConfig
from matplotlib.colors import ListedColormap

class LatticeNode:
    def __init__(self, x, y, z=0, species="X"):
        self.x = x
        self.y = y
        self.z = z
        self.species = species
        self.neighbors = []

class Lattice:
    default_species_colors = {
        "X": "black",
        "O*": "blue",
        "CO*": "red",
        "CO2*": "green"
    }

    def __init__(self, lattice_type="SimpleCubic", nx=None, ny=None, nz=1, initial_species="X", filename=None):
        self.lattice_type = lattice_type
        self.initial_species = initial_species
        if filename:
            self._load_from_dat(filename)
        else:
            if nx is None or ny is None:
                raise ValueError("nx and ny must be provided if filename is not set.")
            self.nx, self.ny, self.nz = nx, ny, nz
            self.nodes = np.empty((nx, ny, nz), dtype=object)
            for x in range(nx):
                for y in range(ny):
                    for z in range(nz):
                        self.nodes[x, y, z] = LatticeNode(x, y, z, initial_species)
        self.species_colors = dict(self.default_species_colors)
        

    def _load_from_dat(self, filename):
        lines = []
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or ':' in line:
                    continue
                lines.append(line.split())
        data = lines
        self.ny = len(data)
        self.nx = len(data[0])
        self.nz = 1
        self.nodes = np.empty((self.nx, self.ny, self.nz), dtype=object)
        for y, row in enumerate(data):
            if len(row) != self.nx:
                raise ValueError(f"Inconsistent row length in {filename} at line {y}")
            for x, token in enumerate(row):
                self.nodes[x, y, 0] = LatticeNode(x, y, 0, token)
        self.lattice_type = filename.split('/')[-1].split('.dat')[0]

    def set_species_color(self, species, color):
        if color in self.species_colors.values():
            raise ValueError(f"Color '{color}' is already in use.")
        self.species_colors[species] = color

    def plot(self, show_grid=True):
        # 1. Map species to integers
        unique_species = sorted(set(node.species for x in range(self.nx) for y in range(self.ny) for node in [self.nodes[x, y, 0]]))
        species_to_int = {sp: i for i, sp in enumerate(unique_species)}
        int_grid = np.zeros((self.nx, self.ny), dtype=int)
        for x in range(self.nx):
            for y in range(self.ny):
                int_grid[x, y] = species_to_int[self.nodes[x, y, 0].species]

        # 2. Build colormap
        color_list = [self.species_colors.get(sp, 'gray') for sp in unique_species]
        cmap = ListedColormap(color_list)

        # 3. Plot with imshow
        plt.imshow(int_grid.T, origin='lower', cmap=cmap, interpolation='nearest')
        if show_grid:
            plt.grid(color='gray', linestyle=':', linewidth=0.5)
        plt.title(f"{self.nx}×{self.ny} surface")
        plt.xlabel("X index")
        plt.ylabel("Y index")
        # 4. Legend
        handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, label=sp, markersize=8)
                for sp, color in zip(unique_species, color_list)]
        plt.legend(handles=handles)
        plt.show()

    def save(self, config: KMCConfig):
        # cfg = KMCConfig("input.kmc") # <-- REMOVE THIS
        config.set_lattice(self.lattice_type, self.nx, self.ny, self.nz, self.initial_species)
        # cfg.save() # <-- REMOVE THIS