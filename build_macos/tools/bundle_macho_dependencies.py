#!/usr/bin/env python3
"""Bundle non-system dylibs for one external Mach-O command-line tool."""

from __future__ import annotations

import argparse
from collections import deque
from pathlib import Path
import shutil
import subprocess


SYSTEM_PREFIXES = ("/System/Library/", "/usr/lib/")


def dependencies(path: Path) -> list[str]:
    output = subprocess.check_output(["otool", "-L", str(path)], text=True)
    result = []
    for line in output.splitlines()[1:]:
        value = line.strip().split(" (", 1)[0]
        if value:
            result.append(value)
    return result


def resolve_dependency(owner: Path, dependency: str) -> Path | None:
    if dependency.startswith(SYSTEM_PREFIXES):
        return None
    if dependency.startswith("@rpath/"):
        return (owner.resolve().parent / dependency.removeprefix("@rpath/")).resolve()
    if dependency.startswith("@loader_path/"):
        return (owner.resolve().parent /
                dependency.removeprefix("@loader_path/")).resolve()
    if dependency.startswith("@executable_path/"):
        return (owner.resolve().parent /
                dependency.removeprefix("@executable_path/")).resolve()
    return Path(dependency).resolve()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    source = args.source.resolve()
    destination = args.destination.resolve()
    lib_dir = destination / "lib"
    lib_dir.mkdir(parents=True, exist_ok=True)
    bundled_tool = destination / source.name
    shutil.copy2(source, bundled_tool)

    queue = deque([source])
    originals: dict[str, Path] = {}
    while queue:
        current = queue.popleft()
        for dependency in dependencies(current):
            resolved = resolve_dependency(current, dependency)
            if resolved is None:
                continue
            if not resolved.is_file():
                raise FileNotFoundError(f"Mach-O dependency not found: {dependency}")
            name = resolved.name
            if name in originals:
                if originals[name] != resolved:
                    raise RuntimeError(f"Conflicting dylib basenames: {resolved} and {originals[name]}")
                continue
            originals[name] = resolved
            shutil.copy2(resolved, lib_dir / name)
            queue.append(resolved)

    targets = [(source, bundled_tool, "@loader_path/lib/")]
    targets.extend((original, lib_dir / name, "@loader_path/")
                   for name, original in originals.items())
    for original, target, prefix in targets:
        if target.suffix == ".dylib":
            subprocess.run(["install_name_tool", "-id",
                            "@loader_path/" + target.name, str(target)], check=True)
        for dependency in dependencies(original):
            resolved = resolve_dependency(original, dependency)
            if resolved is None:
                continue
            if resolved.name in originals:
                subprocess.run(["install_name_tool", "-change", dependency,
                                prefix + resolved.name, str(target)], check=True)
        subprocess.run(["codesign", "--force", "--sign", "-", str(target)],
                       check=True)


if __name__ == "__main__":
    main()
