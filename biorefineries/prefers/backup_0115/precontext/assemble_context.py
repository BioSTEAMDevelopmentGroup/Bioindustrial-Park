# -*- coding: utf-8 -*-
"""
assemble_context.py - Generate Codebase Context Map

Phase 3: Create CODEBASE_CONTEXT.txt with flattened view of v1/ codebase.

This script should be placed in: prefers/v1/precontext/
Run from: prefers/v1/precontext/ directory
"""

import os
from pathlib import Path
from datetime import datetime

# Configuration
SCRIPT_DIR = Path(__file__).parent
V1_DIR = SCRIPT_DIR.parent  # v1/ is parent of precontext/

OUTPUT_FILE = SCRIPT_DIR / "CODEBASE_CONTEXT.txt"

# Exclusion patterns
EXCLUDE_DIRS = {"__pycache__", ".git", ".DS_Store", "backup_0115", "precontext"}
EXCLUDE_FILES = {".pyc", ".pyo", ".DS_Store"}


def should_exclude(path: Path) -> bool:
    """Check if path should be excluded."""
    # Check directory exclusions
    for part in path.parts:
        if part in EXCLUDE_DIRS:
            return True
    
    # Check file exclusions
    if path.suffix in EXCLUDE_FILES:
        return True
    
    return False


def get_relative_path(file_path: Path) -> str:
    """Get path relative to v1/ directory."""
    try:
        return str(file_path.relative_to(V1_DIR))
    except ValueError:
        return str(file_path)


def collect_python_files() -> list[Path]:
    """Recursively collect all .py files in v1/."""
    py_files = []
    
    for root, dirs, files in os.walk(V1_DIR):
        root_path = Path(root)
        
        # Filter out excluded directories
        dirs[:] = [d for d in dirs if d not in EXCLUDE_DIRS]
        
        for file in files:
            if file.endswith('.py'):
                file_path = root_path / file
                if not should_exclude(file_path):
                    py_files.append(file_path)
    
    # Sort for consistent ordering
    py_files.sort()
    return py_files


def generate_context():
    """Generate the CODEBASE_CONTEXT.txt file."""
    print("=" * 60)
    print("CODEBASE CONTEXT ASSEMBLY")
    print("=" * 60)
    
    print(f"\nScanning: {V1_DIR}")
    
    # Collect files
    py_files = collect_python_files()
    print(f"Found {len(py_files)} Python files")
    
    if not py_files:
        print("  [WARN] No Python files found!")
        return
    
    # Generate context
    lines = []
    
    # Header
    lines.append("=" * 70)
    lines.append("PREFERS V1 CODEBASE CONTEXT")
    lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(f"Source: {V1_DIR}")
    lines.append("=" * 70)
    lines.append("")
    
    # Table of contents
    lines.append("TABLE OF CONTENTS")
    lines.append("-" * 40)
    for i, file_path in enumerate(py_files, 1):
        rel_path = get_relative_path(file_path)
        lines.append(f"  {i:2d}. {rel_path}")
    lines.append("")
    
    # File contents
    for file_path in py_files:
        rel_path = get_relative_path(file_path)
        
        lines.append("")
        lines.append("")
        lines.append(f"================ FILE: {rel_path} ================")
        lines.append("")
        
        try:
            content = file_path.read_text(encoding='utf-8', errors='ignore')
            lines.append(content)
            print(f"  [OK] Added: {rel_path} ({len(content)} chars)")
        except Exception as e:
            lines.append(f"# ERROR: Could not read file: {e}")
            print(f"  [ERR] Error: {rel_path} - {e}")
    
    # Write output
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_FILE.write_text("\n".join(lines), encoding='utf-8')
    
    print(f"\n  [OK] Context written to: {OUTPUT_FILE}")
    print(f"  Total size: {OUTPUT_FILE.stat().st_size:,} bytes")


def main():
    generate_context()
    
    print("\n" + "=" * 60)
    print("CONTEXT ASSEMBLY COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()
