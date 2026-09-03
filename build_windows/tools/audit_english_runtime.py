#!/usr/bin/env python3
"""Audit the shipped Designer runtime for untranslated CJK string literals."""

from __future__ import annotations

import argparse
import ast
import importlib.util
import json
from pathlib import Path
import re
import sys


CJK_RE = re.compile(r"[\u3400-\u4dbf\u4e00-\u9fff\uf900-\ufaff]")

PUBLIC_VENDOR_FILES = {
    Path("designer_vendor/cadnano2/model/io/orthogonalseq.py"),
    Path("designer_vendor/cadnano2/model/io/primer3analysis.py"),
    Path("designer_vendor/cadnano2/model/parts/part.py"),
    Path("designer_vendor/cadnano2/views/orthogonalsequences.py"),
    Path("designer_vendor/cadnano2/views/primer3analysis.py"),
}

# This is an operational parser literal, never presentation text.  The
# regular expression reads both legacy Chinese and current English logs.
INTERNAL_LITERAL_ALLOWLIST = {
    r"总长度\D*(\d+)\s*nt",
}


def _load_i18n(path: Path):
    spec = importlib.util.spec_from_file_location("moire_release_i18n", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _iter_runtime_files(root: Path):
    for path in sorted(root.rglob("*.py")):
        parts = set(path.parts)
        if "tests" in parts or "__pycache__" in parts:
            continue
        if path.name == "i18n.py":
            continue
        relative = path.relative_to(root)
        if relative.parts and relative.parts[0] == "designer_vendor" and \
                relative not in PUBLIC_VENDOR_FILES:
            # The Windows product launches a separate, clean cadnano build.
            # Legacy private-engine twist/bending and curved/frame UI modules
            # are not reachable from Designer; only the sequence dialogs and
            # model methods listed above are shipped presentation surfaces.
            continue
        yield path


def audit(root: Path):
    i18n_path = root / "moire_designer" / "i18n.py"
    i18n = _load_i18n(i18n_path)
    records = []
    for path in _iter_runtime_files(root):
        source = path.read_text(encoding="utf-8")
        try:
            tree = ast.parse(source, filename=str(path))
        except SyntaxError as error:
            records.append({
                "file": str(path), "line": error.lineno or 0,
                "value": "<syntax error>", "translated": str(error),
            })
            continue
        for node in ast.walk(tree):
            if not isinstance(node, ast.Constant) or not isinstance(node.value, str):
                continue
            if not CJK_RE.search(node.value):
                continue
            if node.value in INTERNAL_LITERAL_ALLOWLIST:
                continue
            translated = i18n.translate(node.value, "en")
            records.append({
                "file": str(path.relative_to(root)),
                "line": int(getattr(node, "lineno", 0)),
                "value": node.value,
                "translated": translated,
                "resolved": not bool(CJK_RE.search(translated)),
            })
    return records


def audit_companion(root: Path):
    records = []
    for suffix in ("*.py", "*.ui", "*.json"):
        for path in sorted(root.rglob(suffix)):
            if "__pycache__" in path.parts:
                continue
            for line_number, line in enumerate(
                    path.read_text(encoding="utf-8").splitlines(), 1):
                if CJK_RE.search(line):
                    records.append({
                        "file": str(path.relative_to(root)),
                        "line": line_number,
                        "value": line.strip(),
                        "translated": line.strip(),
                        "resolved": False,
                    })
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--json", type=Path)
    parser.add_argument("--unresolved-only", action="store_true")
    parser.add_argument("--companion-root", type=Path)
    args = parser.parse_args()
    records = audit(args.root.resolve())
    i18n = _load_i18n(
        args.root.resolve() / "moire_designer" / "i18n.py")
    catalog_records = [{
        "file": "moire_designer/i18n.py",
        "line": 0,
        "value": source,
        "translated": target,
        "resolved": False,
    } for source, target in i18n._catalogs.get("en", {}).items()
        if CJK_RE.search(str(target))]
    records.extend(catalog_records)
    if args.companion_root:
        records.extend(audit_companion(args.companion_root.resolve()))
    selected = ([item for item in records if not item.get("resolved", False)]
                if args.unresolved_only else records)
    payload = {
        "runtime_cjk_literals": len(records),
        "unresolved": sum(not item.get("resolved", False) for item in records),
        "records": selected,
    }
    rendered = json.dumps(payload, ensure_ascii=False, indent=2)
    if args.json:
        args.json.write_text(rendered + "\n", encoding="utf-8")
    else:
        print(rendered)
    return 1 if payload["unresolved"] else 0


if __name__ == "__main__":
    sys.exit(main())
