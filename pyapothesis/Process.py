"""
Controlling the various reaction rates

"""
from pyapothesis.config import KMCConfig

class Process:
    def __init__(self, reactants, products, rate, reaction_type=None):
        self.reactants = reactants
        self.products = products
        self.rate = rate
        if reaction_type is None:
            # Simple inference
            if "*" in reactants and "*" not in products:
                self.reaction_type = "adsorption"
            elif "*" in products and "*" not in reactants:
                self.reaction_type = "desorption"
            else:
                self.reaction_type = "reaction"
        else:
            self.reaction_type = reaction_type
    
    def __str__(self):
        reactants_str = " + ".join(self.reactants)
        products_str = " + ".join(self.products)
        return f"{reactants_str} -> {products_str}: constant {self.rate}"
    
    
    def __repr__(self):
        return f"Process({self.reaction_type}, {self.reactants}, {self.products}, {self.rate})"
    
    def save(self,config: KMCConfig):
        """
        Save the process to a KMCConfig file.
        """
        if self.reaction_type == "adsorption":
            config.set_adsorption_constant(self.reactants[0], self.rate)
        elif self.reaction_type == "desorption":
            config.set_desorbtion_constant(self.rate)
        elif self.reaction_type == "reaction":
            config.set_reaction_constant(self.rate)
        else:
            raise ValueError(f"Unknown reaction type: {self.reaction_type}")
        


    