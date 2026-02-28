#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate YHapLZ - Integrate A-T complete haplogroup trees

Usage:
    python generate_YHapLZ.py

Required files in current directory:
    ATREE.txt, BTREE.txt, CTREE.txt, DTREE.txt, ETREE.txt, 
    FTREE.txt, GTREE.txt, HTREE.txt, ITREE.txt, JTREE.txt,
    KTREE.txt, LTREE.txt, MTREE.txt, NTREE.txt, OTREE.txt,
    PTREE.txt, QTREE.txt, RTREE.txt, STREE.txt, TTREE.txt
    
    and proYHapLZ.py (base classifier code)

Output:
    YHapLZ.py
"""

import sys
import os
import re

def parse_isogg_tree(file_path):
    """
    Parse ISOGG format tree file
    Returns: [(indent_level, haplogroup, main_SNP), ...]
    """
    nodes = []
    
    if not os.path.exists(file_path):
        print(f"  Warning: {file_path} does not exist, skipping")
        return nodes
    
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            if not line.strip():
                continue
            
            # Count tab indentation
            tabs = 0
            for ch in line:
                if ch == '\t':
                    tabs += 1
                else:
                    break
            
            content = line[tabs:].rstrip('\n\r')
            parts = content.split('\t', 1)
            haplo = parts[0].strip()
            
            # Skip uncertain nodes
            if '~' in haplo or '*' in haplo or not haplo:
                continue
            
            # Extract main SNP
            if len(parts) > 1 and parts[1].strip():
                snp_str = parts[1].strip()
                first_snp = snp_str.split(',')[0].strip()
                first_snp = first_snp.split('/')[0].strip()
                first_snp = first_snp.replace('^^', '').replace('^', '').replace('~', '').strip()
            else:
                first_snp = haplo
            
            nodes.append((tabs, haplo, first_snp))
    
    return nodes

def format_branch(nodes, base_indent):
    """Convert node list to tree string"""
    if not nodes:
        return ""
    lines = []
    for tabs, haplo, snp in nodes:
        indent = '\t' * (base_indent + tabs)
        lines.append(f"{indent}{haplo}\t{snp}")
    return '\n'.join(lines)

def generate_yhaplz():
    """Generate YHapLZ.py"""
    
    print("="*60)
    print("YHapLZ Generator")
    print("Integrating A-T complete haplogroup trees")
    print("="*60)
    
    # Define all haplogroups and their base indentation in tree
    haplogroups = {
        # Haplogroup: (filename, base_indent)
        'A': ('ATREE.txt', 1),
        'B': ('BTREE.txt', 7),
        'C': ('CTREE.txt', 9),
        'D': ('DTREE.txt', 9),
        'E': ('ETREE.txt', 9),
        'F': ('FTREE.txt', 9),
        'G': ('GTREE.txt', 11),
        'H': ('HTREE.txt', 12),
        'I': ('ITREE.txt', 14),
        'J': ('JTREE.txt', 14),
        'K': ('KTREE.txt', 13),
        'L': ('LTREE.txt', 15),
        'M': ('MTREE.txt', 14),
        'N': ('NTREE.txt', 16),
        'O': ('OTREE.txt', 16),
        'P': ('PTREE.txt', 15),
        'Q': ('QTREE.txt', 17),
        'R': ('RTREE.txt', 17),
        'S': ('STREE.txt', 14),
        'T': ('TTREE.txt', 15),
    }
    
    # Parse all tree files
    print("\nParsing tree files...")
    parsed_trees = {}
    for hg, (filename, _) in haplogroups.items():
        nodes = parse_isogg_tree(filename)
        if nodes:
            parsed_trees[hg] = nodes
            print(f"  {hg}: {len(nodes)} nodes")
        else:
            print(f"  {hg}: skipped (file not found or empty)")
    
    # Build complete OFFICIAL_TREE
    tree_parts = []
    
    # === Root and A branch ===
    tree_parts.append('"""\nY\tRoot')
    
    if 'A' in parsed_trees:
        tree_parts.append(format_branch(parsed_trees['A'], 1))
    else:
        tree_parts.append('''\tA0000\tA8864
\tA000-T\tA8835
\t\tA000\tA10805
\t\tA00-T\tPR2921
\t\t\tA00\tAF6/L1284
\t\t\tA0-T\tL1085
\t\t\t\tA0\tCTS2809.1
\t\t\t\tA1\tP305
\t\t\t\t\tA1a\tV161
\t\t\t\t\tA1b\tP108
\t\t\t\t\t\tA1b1\tL419''')
    
    # === BT branch ===
    tree_parts.append('\t\t\t\t\t\tBT\tM91')
    
    if 'B' in parsed_trees:
        tree_parts.append(format_branch(parsed_trees['B'], 7))
    else:
        tree_parts.append('''\t\t\t\t\t\t\tB\tM60
\t\t\t\t\t\t\t\tB1\tM236
\t\t\t\t\t\t\t\tB2\tM182''')
    
    # === CT branch ===
    tree_parts.append('\t\t\t\t\t\t\tCT\tM168')
    
    # === DE branch ===
    tree_parts.append('\t\t\t\t\t\t\t\tDE\tM145')
    
    if 'D' in parsed_trees:
        tree_parts.append(format_branch(parsed_trees['D'], 9))
    else:
        tree_parts.append('\t\t\t\t\t\t\t\t\tD\tM174')
    
    if 'E' in parsed_trees:
        tree_parts.append(format_branch(parsed_trees['E'], 9))
    else:
        tree_parts.append('''\t\t\t\t\t\t\t\t\tE\tM96
\t\t\t\t\t\t\t\t\t\tE1\tP147
\t\t\t\t\t\t\t\t\t\t\tE1a\tM33
\t\t\t\t\t\t\t\t\t\t\tE1b\tP179
\t\t\t\t\t\t\t\t\t\t\t\tE1b1\tP2
\t\t\t\t\t\t\t\t\t\t\t\t\tE1b1a\tM2
\t\t\t\t\t\t\t\t\t\t\t\t\tE1b1b\tM35''')
    
    # === CF branch ===
    tree_parts.append('\t\t\t\t\t\t\t\tCF\tP143')
    
    if 'C' in parsed_trees:
        tree_parts.append(format_branch(parsed_trees['C'], 9))
    else:
        tree_parts.append('''\t\t\t\t\t\t\t\t\tC\tM130
\t\t\t\t\t\t\t\t\t\tC1\tF3393
\t\t\t\t\t\t\t\t\t\tC2\tM217''')
    
    tree_parts.append('\t\t\t\t\t\t\t\t\tF\tM89')
    
    if 'F' in parsed_trees:
        f_nodes = [(t+1, h, s) for t, h, s in parsed_trees['F'] if not h.startswith('F\t')]
        if f_nodes:
            tree_parts.append(format_branch(f_nodes, 9))
    
    # === GHIJK branch ===
    tree_parts.append('\t\t\t\t\t\t\t\t\t\tGHIJK\tF1329')
    
    if 'G' in parsed_trees:
        tree_parts.append(format_branch(parsed_trees['G'], 11))
    else:
        tree_parts.append('''\t\t\t\t\t\t\t\t\t\t\tG\tM201
\t\t\t\t\t\t\t\t\t\t\t\tG1\tM285
\t\t\t\t\t\t\t\t\t\t\t\tG2\tP287''')
    
    # === HIJK branch ===
    tree_parts.append('\t\t\t\t\t\t\t\t\t\t\tHIJK\tF929')
    
    if 'H' in parsed_trees:
        tree_parts.append(format_branch(parsed_trees['H'], 12))
    else:
        tree_parts.append('''\t\t\t\t\t\t\t\t\t\t\t\tH\tL901
\t\t\t\t\t\t\t\t\t\t\t\t\tH1\tM52
\t\t\t\t\t\t\t\t\t\t\t\t\tH2\tP96''')
    
    # === IJK branch ===
    tree_parts.append('\t\t\t\t\t\t\t\t\t\t\t\tIJK\tL15')
    
    # === IJ branch ===
    tree_parts.append('\t\t\t\t\t\t\t\t\t\t\t\t\tIJ\tM429')
    
    if 'I' in parsed_trees:
        tree_parts.append(format_branch(parsed_trees['I'], 14))
    else:
        tree_parts.append('''\t\t\t\t\t\t\t\t\t\t\t\t\t\tI\tM170
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tI1\tM253
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tI2\tM438''')
    
    if 'J' in parsed_trees:
        tree_parts.append(format_branch(parsed_trees['J'], 14))
    else:
        tree_parts.append('''\t\t\t\t\t\t\t\t\t\t\t\t\t\tJ\tM304
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tJ1\tM267
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tJ2\tM172''')
    
    # === K branch ===
    tree_parts.append('\t\t\t\t\t\t\t\t\t\t\t\t\tK\tM9')
    
    # === LT branch ===
    tree_parts.append('\t\t\t\t\t\t\t\t\t\t\t\t\t\tLT\tL298')
    
    if 'L' in parsed_trees:
        tree_parts.append(format_branch(parsed_trees['L'], 15))
    else:
        tree_parts.append('''\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tL\tM20
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tL1\tM22''')
    
    if 'T' in parsed_trees:
        tree_parts.append(format_branch(parsed_trees['T'], 15))
    else:
        tree_parts.append('''\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tT\tM184
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tT1\tL162''')
    
    # === NO branch ===
    tree_parts.append('\t\t\t\t\t\t\t\t\t\t\t\t\t\tNO\tF549')
    tree_parts.append('\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tNO1\tM214')
    
    if 'N' in parsed_trees:
        tree_parts.append(format_branch(parsed_trees['N'], 16))
    else:
        tree_parts.append('''\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tN\tM231
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tN1\tCTS3750
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tN1a\tM2013''')
    
    if 'O' in parsed_trees:
        tree_parts.append(format_branch(parsed_trees['O'], 16))
    else:
        tree_parts.append('\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tO\tM175')
    
    # === K2 branch ===
    tree_parts.append('\t\t\t\t\t\t\t\t\t\t\t\t\t\tK2\tM526')
    tree_parts.append('\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tK2a\tM2308')
    
    if 'M' in parsed_trees:
        tree_parts.append(format_branch(parsed_trees['M'], 16))
    
    if 'S' in parsed_trees:
        tree_parts.append(format_branch(parsed_trees['S'], 16))
    
    # === K2b -> P branch ===
    tree_parts.append('\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tK2b\tP331')
    
    if 'P' in parsed_trees:
        tree_parts.append(format_branch(parsed_trees['P'], 16))
    else:
        tree_parts.append('\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tP\tP295')
    
    if 'Q' in parsed_trees:
        tree_parts.append(format_branch(parsed_trees['Q'], 17))
    else:
        tree_parts.append('''\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tQ\tM242
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tQ1\tL472''')
    
    if 'R' in parsed_trees:
        tree_parts.append(format_branch(parsed_trees['R'], 17))
    else:
        tree_parts.append('''\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tR\tM207
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tR1\tM173
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tR1a\tM420
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tR1b\tM343
\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tR2\tM479''')
    
    tree_parts.append('"""')
    
    # Merge tree parts
    new_tree = '\n'.join(tree_parts)
    
    # Read base classifier
    base_file = 'proYHapLZ.py'
    if not os.path.exists(base_file):
        # Try alternative path
        alt_paths = [
            'github/proYHapLZ.py',
            '../proYHapLZ.py',
            'y_haplogroup_classifier_v4.0.py'
        ]
        for alt in alt_paths:
            if os.path.exists(alt):
                base_file = alt
                break
        else:
            print(f"\nError: Cannot find base file (proYHapLZ.py)")
            print("Please ensure proYHapLZ.py exists in current directory")
            sys.exit(1)
    
    print(f"\nReading {base_file}...")
    with open(base_file, 'r', encoding='utf-8') as f:
        base_code = f.read()
    
    # Replace OFFICIAL_TREE
    pattern = r'OFFICIAL_TREE = """.*?"""'
    new_tree_assignment = f'OFFICIAL_TREE = {new_tree}'
    new_code = re.sub(pattern, new_tree_assignment, base_code, flags=re.DOTALL)
    
    # Update version and naming
    new_code = new_code.replace('VERSION = "4.0"', 'VERSION = "1.0"')
    new_code = new_code.replace("VERSION = '4.0'", "VERSION = '1.0'")
    new_code = new_code.replace('v4.0', 'v1.0')
    new_code = new_code.replace('v4.1', 'v1.0')
    
    # Update tool name in docstrings and prints
    new_code = new_code.replace('YHapLZ v4.1 (YHapLZ)', 
                                 'YHapLZ - YHapLZ v1.0')
    new_code = new_code.replace('Y染色体单倍群分型工具 v4.1', 
                                 'YHapLZ - YHapLZ v1.0')
    new_code = new_code.replace('Y染色体单倍群分型工具 v4.0', 
                                 'YHapLZ - YHapLZ v1.0')
    new_code = new_code.replace('y_haplogroup_classifier_v4.0', 'YHapLZ')
    new_code = new_code.replace('y_haplogroup_classifier_v4.1', 'YHapLZ')
    new_code = new_code.replace('Y_haplogroup_result_v4.1.txt', 'YHapLZ_result.txt')
    new_code = new_code.replace('Y_haplogroup_result_v4.0.txt', 'YHapLZ_result.txt')
    new_code = new_code.replace('Y_haplogroup_summary_v4.1.txt', 'YHapLZ_summary.txt')
    new_code = new_code.replace('Y_haplogroup_summary_v4.0.txt', 'YHapLZ_summary.txt')
    
    # Write YHapLZ.py
    output_file = 'YHapLZ.py'
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(new_code)
    
    # Statistics
    total_nodes = sum(len(nodes) for nodes in parsed_trees.values())
    
    print(f"\n{'='*60}")
    print(f"Generation complete!♡(*´∀｀*)人(*´∀｀*)♡")
    print(f"{'='*60}")
    print(f"  Output: {output_file}")
    print(f"  Total haplogroup branches: {len(parsed_trees)}")
    print(f"  Total tree nodes: {total_nodes}")
    print(f"\nUsage:")
    print(f"  python YHapLZ.py -v <input.vcf> -i <indexdata.csv> -o <output_dir>")
    print(f"\nExample:")
    print(f"  python YHapLZ.py -v sample.vcf -i indexdata.csv -o results/")


if __name__ == "__main__":
    generate_yhaplz()
