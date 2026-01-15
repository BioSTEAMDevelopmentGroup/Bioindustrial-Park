# -*- coding: utf-8 -*-
"""
restore_tests.py - Restore Legacy Test Blocks

Phase 2: Extract `if __name__ == "__main__":` blocks from backup
and inject them into v1 files that lack them.

Run from: prefers/ directory (after migration_setup.py)
"""

import re
from pathlib import Path

# Configuration
SCRIPT_DIR = Path(__file__).parent
PREFERS_ROOT = SCRIPT_DIR

BACKUP_DIR = PREFERS_ROOT / "backup_0115"
V1_DIR = PREFERS_ROOT / "v1"

# Mapping: backup file → v1 target file
FILE_MAPPINGS = {
    "_models.py": "LegH/models.py",
    "_tea.py": "LegH/tea.py",
    "_chemicals.py": "LegH/chemicals.py",
    "_units.py": "units.py",
    "_process_settings.py": "process_settings.py",
}

# Also check system files in systems/ backup
SYSTEM_MAPPINGS = {
    "systems/LegH/LegH.py": "LegH/system.py",
    "systems/HemeIn/HemeIn.py": "HemeIn/system.py",
}

LEGACY_HEADER = "\n\n# --- Legacy Test Block (Restored) ---\n"


def extract_main_block(file_path: Path) -> str | None:
    """
    Extract the 'if __name__ == "__main__":' block from a file.
    Returns the block content or None if not found.
    """
    if not file_path.exists():
        return None
    
    content = file_path.read_text(encoding='utf-8', errors='ignore')
    
    # Find the main block pattern
    pattern = r'(?:^|\n)(if\s+__name__\s*==\s*[\'"]__main__[\'"]\s*:.*)'
    match = re.search(pattern, content, re.DOTALL)
    
    if match:
        return match.group(1)
    
    return None


def has_main_block(file_path: Path) -> bool:
    """Check if file already has a main block."""
    if not file_path.exists():
        return False
    
    content = file_path.read_text(encoding='utf-8', errors='ignore')
    return 'if __name__' in content and '__main__' in content


def inject_main_block(target_path: Path, main_block: str) -> bool:
    """Inject main block into target file."""
    if not target_path.exists():
        print(f"    ⚠️  Target not found: {target_path}")
        return False
    
    content = target_path.read_text(encoding='utf-8', errors='ignore')
    
    # Append the legacy block
    new_content = content.rstrip() + LEGACY_HEADER + main_block
    
    target_path.write_text(new_content, encoding='utf-8')
    return True


def process_mappings(mappings: dict, backup_base: Path):
    """Process a set of file mappings."""
    for backup_file, target_file in mappings.items():
        backup_path = backup_base / backup_file
        target_path = V1_DIR / target_file
        
        print(f"\n  Processing: {backup_file} → {target_file}")
        
        # Extract main block from backup
        main_block = extract_main_block(backup_path)
        
        if main_block is None:
            print(f"    ℹ️  No main block found in backup")
            continue
        
        # Check if target already has main block
        if has_main_block(target_path):
            print(f"    ✓ Target already has main block, skipping")
            continue
        
        # Inject main block
        if inject_main_block(target_path, main_block):
            print(f"    ✓ Injected main block ({len(main_block)} chars)")
        else:
            print(f"    ✗ Failed to inject")


def main():
    print("=" * 60)
    print("RESTORE LEGACY TEST BLOCKS")
    print("=" * 60)
    
    if not BACKUP_DIR.exists():
        print(f"\n  ✗ Backup directory not found: {BACKUP_DIR}")
        print("  Run migration_setup.py first!")
        return
    
    if not V1_DIR.exists():
        print(f"\n  ✗ V1 directory not found: {V1_DIR}")
        print("  Run migration_setup.py first!")
        return
    
    print(f"\nBackup source: {BACKUP_DIR}")
    print(f"V1 target:     {V1_DIR}")
    
    # Process root-level file mappings
    print("\n--- Processing root-level files ---")
    process_mappings(FILE_MAPPINGS, BACKUP_DIR)
    
    # Process system file mappings
    print("\n--- Processing system files ---")
    process_mappings(SYSTEM_MAPPINGS, BACKUP_DIR)
    
    print("\n" + "=" * 60)
    print("RESTORATION COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()
