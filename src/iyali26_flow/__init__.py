"""Likelihood-free parameter-inference experiments for the iYali26 GEM."""

__version__ = "0.1.0"

from .config import ExperimentConfig, load_experiment_config
from .core import ParameterPoint, Phase1Decision, SimulationResult, Simulator
from .experiment import ExperimentRunner
from .simulator import R4R1846CapacitySimulator

__all__ = [
    "ExperimentConfig",
    "ExperimentRunner",
    "ParameterPoint",
    "Phase1Decision",
    "R4R1846CapacitySimulator",
    "SimulationResult",
    "Simulator",
    "load_experiment_config",
]
