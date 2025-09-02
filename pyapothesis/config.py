class KMCConfig:
    # KMCConfig class for managing Kinetic Monte Carlo simulation configurations
    def __init__(self, filename):
        self.filename = filename
        self.lines = []
        self._parse_file()

    def _parse_file(self):
        self.lines = []
        with open(self.filename, 'r') as file:
            for line in file:
                self.lines.append(line.rstrip('\n'))

    def _update_line(self, key, new_values):
        found = False
        for i, line in enumerate(self.lines):
            if line.strip().startswith(key):
                self.lines[i] = f"{key} {new_values}"
                found = True
        if not found:
            self.lines.append(f"{key} {new_values}")

    def set_temperature(self, T):
        self._update_line("temperature:", str(T))

    def set_pressure(self, P):
        self._update_line("pressure:", str(P))

    def set_lattice(self, typ, nx, ny, nz, label):
        self._update_line("lattice:", f"{typ} {nx} {ny} {nz} {label}")

    def set_time(self, start, duration):
        self._update_line("time_start:", str(start))
        self._update_line("time_duration:", str(duration))

    def set_random(self, seed):
        self._update_line("random:", str(seed))

    def speices(self, species):
        if "CO2" in species:
            for i, line in enumerate(self.lines):
                if line.strip().startswith("# growth:"):
                    raise ValueError("Cannot set species 'CO2' when growth is disabled.")
        self._update_line("report:", species)
    
    def set_growth(self, bool):
        if bool:
            # remove # the growth line if it exists
            for i, line in enumerate(self.lines):
                if line.strip().startswith("# growth:"):
                    self.lines[i] = line.lstrip("# ")
                    break
        else:
            # add # the growth line if it exists
            for i, line in enumerate(self.lines):
                if line.strip().startswith("growth:"):
                    self.lines[i] = "# " + line
                    break

    def set_adsorption_constant(self, species, new_constant):
        """
        Update the constant for the adsorption reaction involving the given species.
        Example: species="CO", new_constant=0.8
        """
        for i, line in enumerate(self.lines):
            if line.strip().startswith("#"):
                continue
            if species in line and "->" in line and "constant" in line and "*" in line:
                parts = line.split("constant")
                prefix = parts[0].strip()
                # preserve any trailing text (e.g., 'all')
                trailing = parts[1].split()
                trailing = " ".join(trailing[1:]) if len(trailing) > 1 else ""
                self.lines[i] = f"{prefix}constant {new_constant} {trailing}".strip()
                break
    
    def set_reaction_constant(self, constant):
        """
        Add or update a reaction constant in the format:
        reactants -> products constant <value>
        Example: reactants="CO O2", products="CO2", constant=0.5
        """
        reaction_line = "CO* + O* -> CO2*"
        for i, line in enumerate(self.lines):
            if line.strip().startswith(reaction_line):
                parts = line.split("constant")
                prefix = parts[0].strip()
                self.lines[i] = f"{prefix}constant {constant}".strip()
                return
    
    def set_desorbtion_constant(self, constant):
        """"
        Update the desorption constant for CO* species.
        Example: new_constant=0.1
        """
        reaction_line = "CO2* -> * + CO2:"
        for i, line in enumerate(self.lines):
            if line.strip().startswith(reaction_line):
                parts = line.split("constant")
                prefix = parts[0].strip()
                self.lines[i] = f"{prefix}constant {constant}".strip()
                return
    
    def set_refresh_time(self, rate):
        """
        Set the refresh rate for the simulation.
        Example: rate=100
        """
        self._update_line("write: log ", str(rate))
        self._update_line("write: lattice ", str(rate))

    def save(self):
        with open(self.filename, 'w') as file:
            for line in self.lines:
                file.write(line + '\n')