"""
DNAfold2 Configuration Module

Handles parsing, validation, and management of folding configuration parameters.
"""

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Dict, Any
import re

logger = logging.getLogger(__name__)


@dataclass
class FoldingConfig:
    """Configuration for DNA folding simulations.
    
    Attributes:
        sampling_method: Sampling algorithm - "remc" (Replica Exchange Monte Carlo) 
                        or "sa" (Simulated Annealing)
        folding_steps: Number of steps for REMC/SA simulations during structure folding
        optimizing_steps: Number of steps in energy optimization after sampling
        na_concentration: Concentration of Na+ ions in mM
        mg_concentration: Concentration of Mg2+ ions in mM
        n_structures: Number of predicted 3D structures to output
        n_threads: Number of temperature replicas (1-20, default 15)
        conf_output_freq: Conformation output frequency (default 500)
        print_freq: Screen output frequency (default 10000)
    """
    sampling_method: str = "remc"
    folding_steps: int = 150000
    optimizing_steps: int = 50000
    na_concentration: float = 1000.0
    mg_concentration: float = 0.0
    n_structures: int = 10
    n_threads: int = 4
    conf_output_freq: int = 500
    print_freq: int = 10000
    
    def __post_init__(self):
        """Validate configuration parameters."""
        logger.debug("Validating configuration parameters...")
        if self.sampling_method not in ("remc", "sa"):
            raise ValueError(f"sampling_method must be 'remc' or 'sa', got {self.sampling_method}")
        if self.folding_steps < 1:
            raise ValueError("folding_steps must be positive")
        if self.optimizing_steps < 1:
            raise ValueError("optimizing_steps must be positive")
        if self.na_concentration < 0:
            raise ValueError("na_concentration cannot be negative")
        if self.mg_concentration < 0:
            raise ValueError("mg_concentration cannot be negative")
        if self.n_structures < 1:
            raise ValueError("n_structures must be at least 1")
        if not 1 <= self.n_threads <= 20:
            raise ValueError("n_threads must be between 1 and 20")
        if self.conf_output_freq < 1:
            raise ValueError("conf_output_freq must be positive")
        if self.print_freq < 1:
            raise ValueError("print_freq must be positive")
        logger.debug("Configuration valid: method=%s, steps=%d, threads=%d",
                     self.sampling_method, self.folding_steps, self.n_threads)
    
    @classmethod
    def from_file(cls, config_path: str | Path) -> "FoldingConfig":
        """Load configuration from a config.dat file.
        
        Args:
            config_path: Path to the configuration file
            
        Returns:
            FoldingConfig instance with parsed values
            
        Example config.dat format:
            Sampling 1                 #Sampling 1 means REMC, 2 means SA
            Folding_steps 150000       #Number of REMC/SA steps
            Optimizing_steps 50000     #Number of optimization steps
            CNa 1000                   #Na+ concentration (mM)
            CMg 0                      #Mg2+ concentration (mM)
            Ncout 10                   #Number of output structures
        """
        config_path = Path(config_path)
        logger.debug("Loading configuration from %s", config_path)
        if not config_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {config_path}")
        
        params: Dict[str, Any] = {}
        
        with open(config_path, 'r') as f:
            for line in f:
                # Remove comments and strip whitespace
                line = line.split('#')[0].strip()
                if not line:
                    continue
                
                parts = line.split()
                if len(parts) < 2:
                    continue
                
                key, value = parts[0].lower(), parts[1]
                
                if key == "sampling":
                    params["sampling_method"] = "remc" if value == "1" else "sa"
                elif key == "folding_steps":
                    params["folding_steps"] = int(value)
                elif key == "optimizing_steps":
                    params["optimizing_steps"] = int(value)
                elif key == "cna":
                    params["na_concentration"] = float(value)
                elif key == "cmg":
                    params["mg_concentration"] = float(value)
                elif key == "ncout":
                    params["n_structures"] = int(value)
        
        logger.debug("Loaded configuration: %s", params)
        return cls(**params)
    
    def to_file(self, config_path: str | Path) -> None:
        """Write configuration to a config.dat file.
        
        Args:
            config_path: Path to write the configuration file
        """
        config_path = Path(config_path)
        config_path.parent.mkdir(parents=True, exist_ok=True)
        
        sampling_value = 1 if self.sampling_method == "remc" else 2
        
        content = f"""Sampling {sampling_value}\t\t   #Sampling 1 means REMC, 2 means SA
Folding_steps {self.folding_steps}        #Number of steps for REMC/SA simulations
Optimizing_steps {self.optimizing_steps}     #Number of steps in optimizing progress
CNa {int(self.na_concentration)}   \t           #Concentration of Na (unit:mM)
CMg {int(self.mg_concentration)}\t\t\t   #Concentration of Mg (unit:mM)
Ncout {self.n_structures}                   #Number of predicted 3D structures
Nthreads {self.n_threads}                   #Number of temperature replicas (1-20)
Conf_freq {self.conf_output_freq}              #Conformation output frequency
Print_freq {self.print_freq}             #Screen output frequency
"""
        
        with open(config_path, 'w') as f:
            f.write(content)
        logger.debug("Configuration written to %s", config_path)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary."""
        return {
            "sampling_method": self.sampling_method,
            "folding_steps": self.folding_steps,
            "optimizing_steps": self.optimizing_steps,
            "na_concentration": self.na_concentration,
            "mg_concentration": self.mg_concentration,
            "n_structures": self.n_structures,
            "n_threads": self.n_threads,
            "conf_output_freq": self.conf_output_freq,
            "print_freq": self.print_freq,
        }


def validate_sequence(sequence: str) -> str:
    """Validate and normalize a DNA sequence.
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        Normalized uppercase sequence
        
    Raises:
        ValueError: If sequence contains invalid characters
    """
    logger.debug("Validating sequence (raw length=%d)", len(sequence))
    sequence = sequence.upper().strip()
    
    # Remove whitespace and newlines
    sequence = re.sub(r'\s+', '', sequence)
    
    # Check for valid DNA nucleotides only
    valid_chars = set('ATCG')
    invalid_chars = set(sequence) - valid_chars
    
    if invalid_chars:
        raise ValueError(f"Invalid characters in sequence: {invalid_chars}. "
                        "Only A, T, C, G are allowed.")
    
    if len(sequence) < 4:
        raise ValueError("Sequence must be at least 4 nucleotides long")
    
    logger.debug("Sequence validated: length=%d", len(sequence))
    return sequence
