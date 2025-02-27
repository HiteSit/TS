# Thomson Sampling Environment Guide

## Build/Test Commands
- Run package: `python src/ts_main.py <path-to-json-params-file>.json`
- Run example: `python src/ts_main.py examples/amide_fp_sim.json`
- Install: `pip install -e .`

## Code Style Guidelines
- **Imports**: Standard lib first, third-party second, local modules last
- **Type Hints**: Use Python type hints (`from typing import List, Optional, Union, Dict, Tuple`)
- **Docstrings**: Functions have docstrings with param explanations and return types
- **Naming**: Snake_case for functions/variables, CamelCase for classes
- **Error Handling**: Use appropriate exceptions with message context
- **Structure**: Follow abstract base classes pattern where appropriate
- **Logging**: Use the logger from ts_logger.py, not print statements

## Key Dependencies
- RDKit (>=2022.03)
- Pandas
- NumPy
- Scikit-learn
- Openeye (optional for ROCS scoring)