#!/usr/bin/env python3
"""
universal_rename.py
===================
A script to perform bulk renaming of identifiers across the repository.

This script:
1. Replaces target strings in all text file contents with safety checks
2. Renames directories matching target names
3. Renames files matching target names

Usage:
    python universal_rename.py [target_dir] [--dry-run]
    
    target_dir: Directory to process (default: current directory)
    --dry-run:  Preview changes without applying them
"""

import os
import re
import sys
import argparse
from pathlib import Path

# =============================================================================
# CONFIGURATION
# =============================================================================

RENAME_MAPPINGS = {
    # Old name -> New name
    "LegH": "LegHb",
    "HemeIn": "HemDx",
}

# Directories to ignore completely (glob patterns)
EXCLUDE_DIRS = {
    ".git", 
    ".github",
    "__pycache__", 
    "*.egg-info", 
    "backup*", 
    "venv", 
    "env",
    "node_modules",
    "$RECYCLE.BIN",
    "System Volume Information"
}

# File extensions to process for content replacement
FILE_EXTENSIONS = {".py", ".md", ".rst", ".txt", ".json", ".yaml", ".yml", ".html", ".css", ".js"}

# =============================================================================
# SCRIPT LOGIC
# =============================================================================

def should_exclude(path_name: str) -> bool:
    """Check if a directory or file name matches any exclusion pattern."""
    from fnmatch import fnmatch
    for pattern in EXCLUDE_DIRS:
        if fnmatch(path_name, pattern):
            return True
    return False

def create_safe_pattern(old_names: list[str]) -> re.Pattern:
    """
    Create a regex pattern that matches old names safely.
    
    Safety features:
    1. Negative lookahead (?!...) to prevent double-renaming (LegH -> LegHbb)
    2. Word boundary check (?![a-z]) to prevent substring matching in other words
       (e.g., matching 'LegH' in 'LegHemoglobin')
       
    Returns compiled regex.
    """
    # Sort by length descending to match longer strings first
    sorted_names = sorted(old_names, key=len, reverse=True)
    
    pattern_parts = []
    for old_name in sorted_names:
        new_name = RENAME_MAPPINGS[old_name]
        
        # Determine what suffix constitutes "already renamed"
        suffix_check = ""
        if new_name.startswith(old_name):
            suffix = new_name[len(old_name):]
            if suffix:
                suffix_check = f"(?!{re.escape(suffix)})"
        
        # Combine protections:
        # 1. Don't match if it's already the new name (suffix check)
        # 2. Don't match if followed by other lowercase letters (prevent partial word match)
        #    We allow matching at end of string or followed by non-alpha, upper case, numbers, etc.
        #    This allows "LegH_System" but blocks "LegHemoglobin"
        part = rf"{re.escape(old_name)}{suffix_check}(?![a-z])"
        pattern_parts.append(part)
    
    return re.compile("|".join(pattern_parts))


def replace_content(content: str, pattern: re.Pattern) -> tuple[str, int]:
    """
    Replace all occurrences of old names with new names in content.
    Returns (new_content, replacement_count).
    """
    count = 0
    
    def replacer(match):
        nonlocal count
        old_name_raw = match.group(0)
        
        # We need to find which key matched because the regex is robust
        # The raw match is the text found, which maps to one of our keys
        # giving the highest priority to longest match in our sorted logic
        
        target_key = None
        # Simple lookup strategy since our pattern ensures exact prefix match
        # We check keys again to be sure which one we matched
        for key in sorted(RENAME_MAPPINGS.keys(), key=len, reverse=True):
            if old_name_raw.startswith(key):
                target_key = key
                break
        
        if target_key:
            new_name = RENAME_MAPPINGS[target_key]
            # If the raw match was just the key, replace it
            if old_name_raw == target_key:
                if old_name_raw != new_name:
                    count += 1
                return new_name
            # If for some reason regex matched more (shouldn't with current logic), return as is
        
        return old_name_raw
    
    new_content = pattern.sub(replacer, content)
    return new_content, count


def process_file_content(file_path: Path, pattern: re.Pattern, dry_run: bool) -> int:
    """
    Process a single file, replacing content as needed.
    Returns number of replacements made.
    """
    # Skip ignored names explicitly even if walked
    if should_exclude(file_path.name):
        return 0

    try:
        # Strict UTF-8 to avoid corrupting binaries that look like text
        content = file_path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        # Skip binary files silently or with minor log if verbose
        return 0
    except Exception as e:
        print(f"  WARNING: Could not read {file_path}: {e}")
        return 0
    
    new_content, count = replace_content(content, pattern)
    
    if count > 0:
        if dry_run:
            print(f"  [DRY-RUN] Modify content: {file_path} ({count} replacements)")
        else:
            try:
                file_path.write_text(new_content, encoding="utf-8")
                print(f"  Modified: {file_path} ({count} replacements)")
            except Exception as e:
                print(f"  ERROR: Could not write {file_path}: {e}")
                return 0
    
    return count


def rename_path(old_path: Path, dry_run: bool) -> tuple[Path | None, str | None]:
    """
    Rename a file or directory if its name matches any old name.
    Returns (new_path, old_name_matched) or (None, None).
    """
    name = old_path.name
    
    # Skip if should exclude
    if should_exclude(name):
        return None, None

    new_name = name
    matched_old = None
    
    for old_name in sorted(RENAME_MAPPINGS.keys(), key=len, reverse=True):
        new_name_target = RENAME_MAPPINGS[old_name]
        
        # Check if this old_name exists in the path name
        if old_name in name:
            # SAFETY: Ensure we don't partial-match middle of words recklessly
            # e.g. "LegHemoglobin" should NOT become "LegHbemoglobin"
            
            # We use regex replacement on the name string for safety
            # Pattern: old_name followed by NOT a lowercase letter
            # But here we are dealing with filenames, which might use different conventions
            # Let's use a similar safety pattern: old_name not followed by [a-z]
            
            # Simple python replacement is risky: name.replace("LegH", "LegHb")
            # Let's use re for the name string
            
            name_pattern = rf"{re.escape(old_name)}(?![a-z])"
            
            # Check if we already have the new name preventing replacement?
            # e.g. LegHb_Process. LegH matches, but it is followed by 'b'. 'b' is [a-z]. 
            # So (?![a-z]) will FAIL for "LegHb..." which properly prevents re-renaming.
            # However, logic: "LegH" in "LegHb" -> index 0. Next char 'b'. 'b' is [a-z].
            # So regex won't match. Great.
            
            if re.search(name_pattern, new_name):
                 new_name_regex = re.sub(name_pattern, new_name_target, new_name)
                 if new_name_regex != new_name:
                     new_name = new_name_regex
                     matched_old = old_name
    
    if new_name != name:
        new_path = old_path.parent / new_name
        return new_path, matched_old
    
    return None, None


def collect_paths_to_rename(root_dir: Path) -> list[tuple[Path, Path]]:
    """
    Collect all paths (files and directories) that need renaming.
    Returns paths sorted by depth (deepest first).
    """
    paths_to_rename = []
    
    # Better approach: 2-pass or rely on os.walk default (topdown=True) to prune, 
    # then reverse list for execution.
    
    all_paths_depth_first = []
    
    for dirpath, dirnames, filenames in os.walk(root_dir, topdown=True):
        # Prune excluded directories in-place
        dirnames[:] = [d for d in dirnames if not should_exclude(d)]
        
        current_dir = Path(dirpath)
        
        # Check files
        for filename in filenames:
            if not should_exclude(filename):
                all_paths_depth_first.append(current_dir / filename)
        
        # We add directory itself processing after its children, 
        # but in this loop we can just track it.
        if current_dir != root_dir and not should_exclude(current_dir.name):
            all_paths_depth_first.append(current_dir)

    # Now we have all valid paths. We need to process them deepest first.
    # So we sort by path length descending (or component count).
    all_paths_depth_first.sort(key=lambda p: len(p.parts), reverse=True)
    
    for path in all_paths_depth_first:
        new_path, _ = rename_path(path, dry_run=True)
        if new_path:
            paths_to_rename.append((path, new_path))

    return paths_to_rename


def main():
    parser = argparse.ArgumentParser(description="Universal Rename Script")
    parser.add_argument("target_dir", nargs="?", default=".", help="Directory to process")
    parser.add_argument("--dry-run", "-n", action="store_true", help="Preview changes without applying")
    args = parser.parse_args()

    script_dir = Path(__file__).parent.resolve()
    target_path = Path(args.target_dir).resolve()
    
    if not target_path.exists():
        print(f"ERROR: Target directory not found: {target_path}")
        sys.exit(1)
    
    print("=" * 70)
    print("UNIVERSAL RENAME SCRIPT (IMPROVED)")
    print("=" * 70)
    print(f"Target directory: {target_path}")
    print(f"Dry run: {args.dry_run}")
    print()
    print("Rename mappings:")
    for old, new in RENAME_MAPPINGS.items():
        print(f"  {old} -> {new}")
    print("\nExcluding directories:")
    for pat in sorted(EXCLUDE_DIRS):
        print(f"  {pat}")
    print()
    
    if args.dry_run:
        print("*** DRY RUN MODE - No changes will be made ***")
        print()
    
    pattern = create_safe_pattern(list(RENAME_MAPPINGS.keys()))
    
    # ==========================================================================
    # STEP 1: Content Replacement
    # ==========================================================================
    print("-" * 70)
    print("STEP 1: Content Replacement")
    print("-" * 70)
    
    files_modified = 0
    total_replacements = 0
    
    # We use os.walk with prune logic again for content
    for dirpath, dirnames, filenames in os.walk(target_path, topdown=True):
        dirnames[:] = [d for d in dirnames if not should_exclude(d)]
        
        for filename in filenames:
            file_path = Path(dirpath) / filename
            
            if should_exclude(filename):
                continue
                
            if file_path.suffix.lower() not in FILE_EXTENSIONS:
                continue
            
            # Don't process script itself
            if file_path.resolve() == Path(__file__).resolve():
                continue

            count = process_file_content(file_path, pattern, args.dry_run)
            if count > 0:
                files_modified += 1
                total_replacements += count
    
    print()
    print(f"Content replacement summary:")
    print(f"  Files modified: {files_modified}")
    print(f"  Total string replacements: {total_replacements}")
    print()
    
    # ==========================================================================
    # STEP 2: Path Renaming
    # ==========================================================================
    print("-" * 70)
    print("STEP 2: Path Renaming")
    print("-" * 70)
    
    paths_to_rename = collect_paths_to_rename(target_path)
    
    dirs_renamed = 0
    files_renamed = 0
    
    for old_path, new_path in paths_to_rename:
        is_dir = old_path.is_dir()
        path_type = "directory" if is_dir else "file"
        
        # Don't rename script inputs!
        if old_path.resolve() == Path(__file__).resolve():
            continue

        if args.dry_run:
            print(f"  [DRY-RUN] Rename {path_type}:")
            print(f"            {old_path.name} -> {new_path.name}")
            if is_dir: dirs_renamed += 1
            else: files_renamed += 1
        else:
            try:
                # Re-check existence, parent might have moved if we did this naively
                # But since we sorted by depth (deepest first), children move before parents
                # So renaming parent last is correct.
                # However, pathlib rename handles full path. 
                # If we rename child: /a/b/child -> /a/b/new
                # Then parent: /a/b -> /a/new_b
                # This works.
                
                old_path.rename(new_path)
                print(f"  Renamed {path_type}: {old_path.name} -> {new_path.name}")
                if is_dir: dirs_renamed += 1
                else: files_renamed += 1
            except Exception as e:
                print(f"  ERROR: Could not rename {old_path}: {e}")
    
    print()
    print(f"Path renaming summary:")
    print(f"  Directories renamed: {dirs_renamed}")
    print(f"  Files renamed: {files_renamed}")
    print()
    
    # ==========================================================================
    # FINAL SUMMARY
    # ==========================================================================
    print("=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)
    print(f"  Files with content changes: {files_modified}")
    print(f"  Total string replacements:  {total_replacements}")
    print(f"  Directories renamed:        {dirs_renamed}")
    print(f"  Files renamed:              {files_renamed}")
    print()
    
    if args.dry_run:
        print("This was a DRY RUN. To apply changes, run without --dry-run")
        print(f"  python {Path(__file__).name} {args.target_dir if args.target_dir != '.' else ''}")
    else:
        print("Renaming complete!")


if __name__ == "__main__":
    main()
