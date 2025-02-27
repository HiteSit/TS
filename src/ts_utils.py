from typing import List, Union, Optional, Dict
from pathlib import Path
import importlib

from .reagent import Reagent


def construct_json(reag_smiles: List[Path],
                   reaction,
                   ts_mode: Optional[str] = "maximize",
                   warmup_trials: Optional[int] = 10,
                   ts_iterations: Optional[int] = 5000,
                   evaluator_class: str = "FPEvaluator",
                   evaluator_args: Dict = None):
    """
    Constructs a JSON-like dictionary for a Thompson Sampling experiment.

    Args:
        reag_smiles (List[Path]): List of file paths containing SMILES strings for reagents.
        reaction (Union[str, Path]): Reaction SMARTS or file path to the reaction file.
        ts_mode (Optional[str]): Mode for Thompson Sampling, either "maximize" or "minimize". Default is "maximize".
        warmup_trials (Optional[int]): Number of warmup trials before Thompson Sampling. Default is 10.
        ts_iterations (Optional[int]): Number of Thompson Sampling iterations. Default is 5000.
        evaluator_class (str): Name of the evaluator class to use.
            - **FPEvaluator**
            - **ROCSEvaluator**
            - **FredEvaluator**
        evaluator_args (Dict): Arguments for the evaluator class. Default is None.
            - If "evaluator_class" is **FPEvaluator**
                - {"query_smiles": "SMILES"}
            - If "evaluator_class" is **ROCSEvaluator**
                - {"query_molfile": Path}
            - If "evaluator_class" is **FredEvaluator**
                - {"design_unit_file": Path}

    Returns:
        Dict[str, Any]: A dictionary containing the configuration for the Thompson Sampling experiment.
    """

    from .thompson_sampling import BaseReaction

    if isinstance(reaction, (str, Path)):
        json_file = {
            "reagent_file_list": [str(x) for x in reag_smiles],
            "reaction_smarts": str(reaction),
            "num_warmup_trials": warmup_trials,
            "num_ts_iterations": ts_iterations,
            "evaluator_class_name": evaluator_class,
            "evaluator_arg": evaluator_args,
            "ts_mode": ts_mode,
            "log_filename": "ts_logs.txt"
        }
        return json_file
    elif isinstance(reaction, BaseReaction):
        mymodule = importlib.import_module(reaction.__module__)
        myclass = getattr(mymodule, reaction.__class__.__name__)
        rxn_input = myclass()

        json_file = {
            "reagent_file_list": [str(x) for x in reag_smiles],
            "reaction_smarts": rxn_input,
            "num_warmup_trials": warmup_trials,
            "num_ts_iterations": ts_iterations,
            "evaluator_class_name": evaluator_class,
            "evaluator_arg": evaluator_args,
            "ts_mode": ts_mode,
            "log_filename": "ts_logs.txt"
        }
        return json_file


def cansmi(smile, isomeric=True, kekule=True):
    """
    Generate a canonical SMILES representation of a oemolecule.

    Args:
        oemol (oechem.OEMol): The oemolecule to generate SMILES for.
        isomeric (bool): Flag indicating whether to include isomeric information in the SMILES.
        kekule (bool): Flag indicating whether to kekulize the oemolecule before generating SMILES.

    Returns:
        str: The canonical SMILES representation of the oemolecule.
    """
    from openeye import oechem
    oemol = oechem.OEMol()
    oechem.OESmilesToMol(oemol, smile)

    oechem.OEFindRingAtomsAndBonds(oemol)
    oechem.OEAssignAromaticFlags(oemol, oechem.OEAroModel_OpenEye)
    smiflag = oechem.OESMILESFlag_Canonical
    if isomeric:
        smiflag |= oechem.OESMILESFlag_ISOMERIC

    if kekule:
        for bond in oemol.GetBonds(oechem.OEIsAromaticBond()):
            bond.SetIntType(5)
        oechem.OECanonicalOrderAtoms(oemol)
        oechem.OECanonicalOrderBonds(oemol)
        oechem.OEClearAromaticFlags(oemol)
        oechem.OEKekulize(oemol)

    # Strip Salt
    oechem.OEDeleteEverythingExceptTheFirstLargestComponent(oemol)

    smile = oechem.OECreateSmiString(oemol, smiflag)
    return smile

def create_reagents(filename: str, num_to_select: Optional[int] = None) -> List[Reagent]:
    """
    Creates a list of Reagents from a file
    :param filename: a smiles file containing the reagents
    :param num_to_select: For dev purposes; the number of molecules to return
    :return: List of Reagents
    """
    reagent_list = []
    with open(filename, 'r') as f:
        for line in f.readlines():
            smiles, reagent_name = line.split()
            reagent = Reagent(reagent_name=reagent_name, smiles=smiles)
            reagent_list.append(reagent)
    if num_to_select is not None and len(reagent_list) > num_to_select:
        reagent_list = reagent_list[:num_to_select]
    return reagent_list


def read_reagents(reagent_file_list, num_to_select: Optional[int]) -> List[Reagent]:
    """
    Read the reagents SMILES files
    :param reagent_file_list: a list of filenames containing reagents for the reaction. Each file list contains smiles
    strings for a single component of the reaction.
    :param num_to_select: select how many reagents to read, mostly a development function
    :return: List of Reagents
    """
    reagents = []
    for reagent_filename in reagent_file_list:
        reagent_list = create_reagents(filename=reagent_filename, num_to_select=num_to_select)
        reagents.append(reagent_list)
    return reagents
