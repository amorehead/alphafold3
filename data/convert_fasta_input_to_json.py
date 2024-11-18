import argparse
import json
import logging
import string

from collections import defaultdict
from itertools import product
from typing import Dict, List

logging.basicConfig(
    format="[%(asctime)s] {%(filename)s:%(lineno)d} %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def read_fasta(filename: str) -> Dict[str, str]:
    """Read sequences from a FASTA file and return them as a dictionary.

    :param filename: Path to the input FASTA file.
    :return: Dictionary containing sequences with identifiers as keys.
    """
    sequences = {}
    current_id = ""
    with open(filename) as file:
        for line in file:
            if line.startswith(">"):
                current_id = line.strip()[1:]
                sequences[current_id] = ""
            else:
                sequences[current_id] += line.strip()
    return sequences


def combine_sequences(sequences: Dict[str, str]) -> Dict[str, str]:
    """Combine sequences of the same protein complex by their PDB codes.

    :param sequences: Dictionary containing sequences with identifiers as keys.
    :return: Dictionary containing combined sequences with PDB codes as keys.
    """
    combined_sequences = defaultdict(list)
    for identifier, sequence in sequences.items():
        pdb_code = identifier.split("_chain")[0]
        chain_id = identifier.split("chain_")[-1]
        combined_sequences[pdb_code].append((chain_id, sequence))
    return combined_sequences


def generate_chain_ids() -> List[str]:
    """
    Generate chain IDs in the format of single letter IDs followed by double letter IDs.

    :return: List of chain IDs.
    """
    # Generate single letter IDs
    single_letter_ids = list(string.ascii_uppercase)

    # Generate double letter IDs in reverse spreadsheet order
    double_letter_ids = []
    for first, second in product(string.ascii_uppercase, repeat=2):
        double_letter_ids.append(second + first)

    # Combine both lists
    return single_letter_ids + double_letter_ids


def write_combined_json(combined_sequences: Dict[str, str], output_filename: str):
    """Write combined sequences to an output JSON file.

    :param combined_sequences: Dictionary containing combined sequences with PDB codes as keys.
    :param output_filename: Path to the output JSON file.
    """
    chain_ids = generate_chain_ids()

    for pdb_code, seq_list in combined_sequences.items():
        output_filepath = output_filename.replace(".json", f"_{pdb_code}.json")
        with open(output_filepath, "w") as file:
            output_data = {
                "name": pdb_code,
                "sequences": [
                    {
                        "protein": {
                            "id": [chain_ids[chain_index]],
                            "sequence": seq[1],
                        }
                    }
                    for chain_index, seq in enumerate(seq_list)
                ],
                "modelSeeds": [1],
                "dialect": "alphafold3",
                "version": 1,
            }

            json.dump(output_data, file)


def main():
    """Read protein chain sequences from a FASTA file, combine sequences of the same protein
    complex, and write the combined sequences to an output JSON file.
    """
    parser = argparse.ArgumentParser(
        description="Convert FASTA input file to JSON input file."
    )
    parser.add_argument(
        "--input_fasta_file",
        type=str,
        help="Path to the input FASTA file.",
        default="data/posebusters_benchmark/sequences/reference_posebusters_benchmark_sequences.fasta",
    )
    parser.add_argument(
        "--output_json_file",
        type=str,
        help="Path to the output JSON file.",
        default="data/posebusters_benchmark/sequences/fold_input.json",
    )
    args = parser.parse_args()

    sequences = read_fasta(args.input_fasta_file)
    combined_sequences = combine_sequences(sequences)
    write_combined_json(combined_sequences, args.output_json_file)

    logger.info("AlphaFold 3 JSON input file preparation complete.")


if __name__ == "__main__":
    main()
