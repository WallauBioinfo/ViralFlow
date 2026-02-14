#!/usr/bin/env python3
import os
import sys

def configure_singularity():
    """
    Configures singularity.conf to work with ViralFlow
    Enables sif fuse to allow SIF container mounting
    """
    # Detect singularity.conf path
    conda_prefix = os.environ.get('CONDA_PREFIX', '')
    
    if not conda_prefix:
        print("WARNING: CONDA_PREFIX not found. Skipping Singularity configuration.")
        print("If you are using micromamba/conda, activate the environment first.")
        return
    
    singularity_conf = os.path.join(conda_prefix, 'etc', 'singularity', 'singularity.conf')
    
    if not os.path.exists(singularity_conf):
        print(f"WARNING: {singularity_conf} not found. Skipping configuration.")
        return
    
    print(f"Configuring Singularity at: {singularity_conf}")
    
    # Create backup of original file
    backup_file = singularity_conf + '.backup'
    if not os.path.exists(backup_file):
        with open(singularity_conf, 'r') as original:
            with open(backup_file, 'w') as backup:
                backup.write(original.read())
        print(f"Backup created at: {backup_file}")
    
    # Read current file
    with open(singularity_conf, 'r') as f:
        lines = f.readlines()
    
    # Remove ALL existing sif fuse lines and check if already configured
    new_lines = []
    found_lines = []
    already_correct = False
    
    for line in lines:
        # Check if this line contains sif fuse
        if 'sif fuse' in line.lower():
            found_lines.append(line.strip())
            # Check if it's already set to yes (and not commented)
            if not line.strip().startswith('#') and '= yes' in line:
                already_correct = True
            # Skip this line (we'll add it back correctly later)
            continue
        else:
            # Keep all other lines
            new_lines.append(line)
    
    # If already correctly configured, just inform and keep one correct line
    if already_correct and len(found_lines) == 1:
        print(f"  ✓ sif fuse is already configured correctly")
        new_lines.append('sif fuse = yes\n')
        # No need to write, file is already correct
        return
    
    # Add the correct configuration at the end
    if found_lines:
        print(f"  ✓ Found {len(found_lines)} existing sif fuse line(s), replacing with correct configuration")
    else:
        print(f"  ✓ No existing sif fuse configuration found, adding new line")
    
    # Add a comment and the correct configuration
    new_lines.append('sif fuse = yes\n')
    
    # Write back the corrected file
    with open(singularity_conf, 'w') as f:
        f.writelines(new_lines)
    
    print("✓ Singularity configuration completed successfully!")
    print(f"  Original file backed up at: {backup_file}")

if __name__ == '__main__':
    configure_singularity()