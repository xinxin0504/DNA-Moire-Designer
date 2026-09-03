"""Small process-isolated entry point for ATHENA's original PDBGen."""

import importlib.util
import os
import sys


def main():
    cndo_path = os.path.abspath(sys.argv[1])
    output_dir = os.path.abspath(sys.argv[2])
    module_path = os.path.join(os.path.dirname(__file__), "pdbgen.py")
    spec = importlib.util.spec_from_file_location("cadnano_athena_pdbgen",
                                                   module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    stem = os.path.splitext(os.path.basename(cndo_path))[0]
    input_dir = os.path.dirname(cndo_path) + os.path.sep
    output_dir = output_dir + os.path.sep
    log_path = os.path.join(output_dir, stem + "-pdbgen.log")
    with open(log_path, "w", encoding="utf-8") as log:
        module.pdbgen(stem, "B", "DNA", input_dir, output_dir, log)
    return 0


if __name__ == "__main__":
    sys.exit(main())
