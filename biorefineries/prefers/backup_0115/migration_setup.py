# -*- coding: utf-8 -*-
"""
migration_setup.py - Safe Migration to v1 Production Environment

Phase 1: Backup & Separation
1. Create backup_0115/ with full copy of current state
2. Create v1/ production directory
3. Move refactored artifacts into v1/

Run from: prefers/ directory
"""

import os
import shutil
from pathlib import Path

# Configuration
SCRIPT_DIR = Path(__file__).parent
PREFERS_ROOT = SCRIPT_DIR  # Script should be in prefers/

BACKUP_DIR = PREFERS_ROOT / "backup_0115"
V1_DIR = PREFERS_ROOT / "v1"

# Artifacts to move to v1
FOLDERS_TO_MOVE = ["LegH", "HemeIn"]
FILES_TO_MOVE = ["units.py", "process_settings.py"]

# Exclusions for backup (don't copy these)
EXCLUDE_PATTERNS = ["__pycache__", ".git", "backup_*", "v1"]


def should_exclude(path: Path) -> bool:
    """Check if path should be excluded from backup."""
    for pattern in EXCLUDE_PATTERNS:
        if pattern in str(path):
            return True
    return False


def phase1_backup():
    """Create immutable backup of current state."""
    print("\n" + "=" * 60)
    print("PHASE 1.1: CREATING BACKUP")
    print("=" * 60)
    
    if BACKUP_DIR.exists():
        print(f"  ⚠️  Backup directory already exists: {BACKUP_DIR}")
        response = input("  Overwrite? (y/N): ").strip().lower()
        if response != 'y':
            print("  Skipping backup...")
            return False
        shutil.rmtree(BACKUP_DIR)
    
    BACKUP_DIR.mkdir(exist_ok=True)
    print(f"  Created: {BACKUP_DIR}")
    
    # Copy all files and folders (except exclusions)
    copied_count = 0
    for item in PREFERS_ROOT.iterdir():
        if should_exclude(item):
            print(f"  Skipping: {item.name}")
            continue
        
        dest = BACKUP_DIR / item.name
        if item.is_dir():
            shutil.copytree(item, dest, ignore=shutil.ignore_patterns(*EXCLUDE_PATTERNS))
            print(f"  Copied folder: {item.name}/")
        else:
            shutil.copy2(item, dest)
            print(f"  Copied file: {item.name}")
        copied_count += 1
    
    print(f"\n  ✓ Backup complete: {copied_count} items copied")
    return True


def phase2_create_v1():
    """Create v1 production directory structure."""
    print("\n" + "=" * 60)
    print("PHASE 1.2: CREATING V1 ENVIRONMENT")
    print("=" * 60)
    
    if V1_DIR.exists():
        print(f"  ⚠️  v1 directory already exists: {V1_DIR}")
        response = input("  Overwrite? (y/N): ").strip().lower()
        if response != 'y':
            print("  Skipping v1 creation...")
            return False
        shutil.rmtree(V1_DIR)
    
    V1_DIR.mkdir(exist_ok=True)
    print(f"  Created: {V1_DIR}")
    
    # Create __init__.py for v1 package
    init_content = '''# -*- coding: utf-8 -*-
"""
PREFERS v1 - Production Package

Versioned, modular BioSTEAM package for precision fermentation simulations.
"""

from . import units
from . import process_settings
from . import LegH
from . import HemeIn

__all__ = ['units', 'process_settings', 'LegH', 'HemeIn']
'''
    (V1_DIR / "__init__.py").write_text(init_content)
    print(f"  Created: v1/__init__.py")
    
    return True


def phase3_migrate_artifacts():
    """Move refactored artifacts into v1."""
    print("\n" + "=" * 60)
    print("PHASE 1.3: MIGRATING ARTIFACTS TO V1")
    print("=" * 60)
    
    moved_count = 0
    
    # Move folders
    for folder_name in FOLDERS_TO_MOVE:
        src = PREFERS_ROOT / folder_name
        dest = V1_DIR / folder_name
        
        if src.exists():
            if dest.exists():
                print(f"  ⚠️  Destination exists, removing: {dest}")
                shutil.rmtree(dest)
            shutil.move(str(src), str(dest))
            print(f"  Moved folder: {folder_name}/ → v1/{folder_name}/")
            moved_count += 1
        else:
            print(f"  ⚠️  Source not found: {src}")
    
    # Move files
    for file_name in FILES_TO_MOVE:
        src = PREFERS_ROOT / file_name
        dest = V1_DIR / file_name
        
        if src.exists():
            shutil.move(str(src), str(dest))
            print(f"  Moved file: {file_name} → v1/{file_name}")
            moved_count += 1
        else:
            print(f"  ⚠️  Source not found: {src}")
    
    print(f"\n  ✓ Migration complete: {moved_count} items moved")
    return True


def main():
    print("=" * 60)
    print("PREFERS MIGRATION SETUP")
    print("=" * 60)
    print(f"\nWorking directory: {PREFERS_ROOT}")
    print(f"Backup target:     {BACKUP_DIR}")
    print(f"V1 target:         {V1_DIR}")
    
    # Execute phases
    phase1_backup()
    phase2_create_v1()
    phase3_migrate_artifacts()
    
    print("\n" + "=" * 60)
    print("MIGRATION COMPLETE")
    print("=" * 60)
    print("\nNext steps:")
    print("  1. Run restore_tests.py to restore legacy test blocks")
    print("  2. Run assemble_context.py to generate codebase context")


if __name__ == "__main__":
    main()
