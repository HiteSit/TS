"""
Thompson Sampling Package

This package provides tools for molecular optimization using Thompson Sampling.
"""

# Import key classes and functions to expose at the top level
from .src.thompson_sampling import BaseReaction, ThompsonSampler
from .src.ts_utils import construct_json, cansmi, create_reagents, read_reagents
from .src.ts_main import run_ts, parse_input_dict
from .src.reagent import Reagent
from .src.evaluators import FPEvaluator, ROCSEvaluator, FredEvaluator
