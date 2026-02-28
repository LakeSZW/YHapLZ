#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Y-chromosome Haplogroup Classifier v4.0
Optimized for low-density SNP data with PLINK compatibility

v4.0Core improvements:
1. Combines v2.4 ancestor chain validation + v3.4 clear framework
2. Strict classification logic: traverse tree downward, validate each node
3. Supports downstream inference: downstream derived + upstream missing → can infer upstream
4. Conflict detection: downstream derived + upstream ancestral → mark warning
5. Detailed diagnostic output for review and debugging

Classification logic (simulates manual classification):
1. Start from root node, check defining SNPs at each level
2. Derived → enter node, continue downward
3. Ancestral → stop, do not enter branch
4. Missing → check downstream for derived evidence (can infer)
5. Final result = deepest validated node
"""

import os
import sys
import argparse
import gzip
from collections import defaultdict
from datetime import datetime

VERSION = "4.0"

# ============================================================
# Phylogenetic tree definition (based on ISOGG 2019-2020)
# Using indented text format for easy maintenance and review
# ============================================================
OFFICIAL_TREE = """
Y	Root
	A0000	A8864
	A000-T	A8835
		A000	A10805
		A00-T	PR2921
			A00	AF6/L1284
			A0-T	L1085
				A0	CTS2809.1/L991.1
				A1	P305
					A1a	V161
					A1b	P108
						A1b1	L419/PF712
						BT	M91
							B	M60
							CT	M168/PF1416
								DE	M145/P205/PF1444
									D	M174/Page3
										D1	CTS3946/Z3707
											D1a	CTS9290/F2285
												D1a1	F6251/Z27276
													D1a1a	F2889
														D1a1a1	F901
															D1a1a1a	Z27269
																D1a1a1a1	Z29448
																	D1a1a1a1a	PH344
																		D1a1a1a1a1	Z44028
																			D1a1a1a1a1a	Z44056
																			D1a1a1a1a1b	F738
																	D1a1a1a1b	F5498
															D1a1a1b	Z41260
																D1a1a1b1	Z41261
																	D1a1a1b1a	Z41262
																		D1a1a1b1a1	Z27277
																		D1a1a1b1a2	Z41250
														D1a1a2	F5439
													D1a1b	Z34171
														D1a1b1	Z34166
															D1a1b1a	Y14750
																D1a1b1a1	SK554
																	D1a1b1a1a	MF77628
												D1a2	Z27276
													D1a2a	P37.1
														D1a2a1	CTS2728
															D1a2a1a	CTS7187
															D1a2a1b	Z17320
															D1a2a1c	Z14850
																D1a2a1c1	Z14879
																	D1a2a1c1a	CTS321
																		D1a2a1c1a1	CTS11048
																			D1a2a1c1a1a	Z39387
																			D1a2a1c1a1b	Z40671
																				D1a2a1c1a1b1	Z40686
														D1a2a2	FGC14584
															D1a2a2a	FGC14590
																D1a2a2a1	FGC14591
																D1a2a2a2	Z43926
																	D1a2a2a2a	Z43932
																	D1a2a2a2b	CTS5960
										D2	L1378
									E	M96/PF182
										E1	P147
											E1a	M33
											E1b	P179
												E1b1	P2
													E1b1a	M2
													E1b1b	M35
								CF	P143/PF2587
									C	M130/Page51/RPS4Y711
										C1	F3393
											C1a	CTS11043
											C1b	K281
												C1b1	F1370
													C1b1a	M356
														C1b1a1	P92
															C1b1a1a	M94
																C1b1a1a1	F948
										C2	M217
											C2a	F1144
												C2a1	F3553
													C2a1a	M407
														C2a1a1	P39
														C2a1a2	BY66307
									F	M89/PF2746
										GHIJK	F1329/M3658/PF2622
											G	M201/PF2957
												G1	M285
												G2	P287
													G2a	P15
											HIJK	F929/M578/PF3494
												H	L901/M2939
												IJK	L15/M523/PF3492
													IJ	M429/P125/PF3535
														I	M170/PF3715
															I1	M253
															I2	M438
														J	M304/Page16/PF4609
															J1	M267
															J2	M172
													K	M9
														LT	L298/P326
															L	M20/PF5570
															T	M184/Page34
														NO	F549/M2335
															NO1	M214/Page39
																N	M231/Page91
																	N1	CTS3750
																		N1a	M2013
																			N1a1	M2291
																				N1a1a	CTS5858
																					N1a1a1	M178
																O	M175
																	O1	M119
																		O1a	M307
																		O1b	P31
																			O1b1	M268
																			O1b2	M176
																	O2	M122
																		O2a	P201
																			O2a1	L465
																			O2a2	M324
																				O2a2a	F1336
																				O2a2b	P201.1
																					O2a2b1	M134
																						O2a2b1a	M117
																							O2a2b1a1	F5
																							O2a2b1a2	CTS7634
																						O2a2b1b	F48
																					O2a2b2	F742
														P	P295/PF5866
															P1	M45
																Q	M242
																	Q1	L232
																R	M207/Page37/UTY2
																	R1	M173
																		R1a	M420
																		R1b	M343
																	R2	M479
"""


class GenotypeStatus:
    """Genotype status enumeration"""
    DERIVED = "derived"      # Derived
    ANCESTRAL = "ancestral"  # Ancestral
    HETEROZYGOUS = "het"     # Heterozygous
    MISSING = "missing"      # Missing
    MISMATCH = "mismatch"    # Mismatch


class HaplogroupClassifier:
    """
    Y-chromosome Haplogroup Classifier
    
    Core logic: simulates manual classification
    1. Traverse phylogenetic tree from root downward
    2. Check defining SNP status at each node
    3. Derived → enter, Ancestral → stop, Missing → check downstream
    4. Validate ancestor chain consistency
    """
    
    # Terminal main branches (includes all major global branches)
    # Note: F is ancestor of GHIJK, most cases should be subdivided to G/H/I/J/K downstream
    # True basal F is rare, handled by special logic
    MAIN_BRANCHES = ['A', 'B', 'C', 'D', 'E', 'G', 'H', 'I', 'J', 'L', 'N', 'O', 'Q', 'R', 'T']
    
    # Intermediate nodes (not used as final result unless cannot subdivide)
    INTERMEDIATE_NODES = {
        'K', 'NO', 'NO1', 'P', 'P1', 'LT', 'F', 'CF', 'DE', 'CT', 
        'BT', 'GHIJK', 'HIJK', 'IJK', 'IJ', 'A1b', 'A1b1',
        'A0-T', 'A00-T', 'A000-T'
    }
    
    # Main branch defining SNPs
    MAIN_BRANCH_DEFINING_SNPS = {
        'A': 'L1085',  # A0-TDefining SNP，Covers all A lineages
        'B': 'M60',    # BHaplogroup defining SNP
        'C': 'M130', 'D': 'M174', 'E': 'M96',
        'G': 'M201', 'H': 'L901', 'I': 'M170',
        'J': 'M304', 'L': 'M20', 'N': 'M231',
        'O': 'M175',  # indel，May not be detectable
        'Q': 'M242', 'R': 'M207', 'T': 'M184'
    }
    
    # FHaplogroupSpecial handling（basal FDetection）
    F_DEFINING_SNP = 'M89'
    
    def __init__(self, isogg_file, het_mode="moderate"):
        """
        Initialize classifier
        
        het_mode:
            strict:   Heterozygous treated as missing (most conservative)
            moderate: Heterozygous treated as derived but marked (default)
            lenient:  Heterozygous fully treated as derived
        """
        self.het_mode = het_mode
        
        print(f"[1] Parsing phylogenetic tree...")
        self.phylo_tree, self.node_snps = self._parse_tree(OFFICIAL_TREE)
        print(f"    Node count: {len(self.phylo_tree)}")
        
        print(f"\n[2] Reading ISOGG index: {isogg_file}")
        self.snp_to_info = {}      # snp_name -> (pos, ref, alt)
        self.pos_to_haplo = defaultdict(list)  # pos -> [(haplo, snp, ref, alt), ...]
        self._load_isogg(isogg_file)
        
        # Build node to SNP position mapping
        self._build_node_snp_positions()
        
        # Check main branch SNPs
        self._check_main_branch_snps()
    
    def _parse_tree(self, tree_text):
        """Parse indented format phylogenetic tree"""
        phylo_tree = {}  # child -> parent
        node_snps = {}   # node -> [snp1, snp2, ...]
        
        lines = tree_text.strip().split('\n')
        path_stack = []  # [(indent, node_name), ...]
        
        for line in lines:
            if not line.strip():
                continue
            
            # Calculate indentation
            indent = 0
            for char in line:
                if char == '\t':
                    indent += 1
                else:
                    break
            
            content = line.strip()
            parts = content.split('\t')
            
            if len(parts) < 1:
                continue
            
            node_name = parts[0].strip()
            
            # Skip root node and invalid lines
            if node_name in ['Y', 'Root', ''] or 'see ' in node_name.lower():
                continue
            
            # Extract SNPs
            snps = []
            if len(parts) >= 2:
                snp_text = parts[1].strip()
                for snp_part in snp_text.split(','):
                    snp_part = snp_part.strip()
                    if snp_part and not snp_part.startswith('('):
                        for alias in snp_part.split('/'):
                            alias = alias.strip().rstrip('~^*')
                            if alias:
                                snps.append(alias)
            
            node_snps[node_name] = snps
            
            # Maintain path stack
            while path_stack and path_stack[-1][0] >= indent:
                path_stack.pop()
            
            if path_stack:
                parent = path_stack[-1][1]
                phylo_tree[node_name] = parent
            
            path_stack.append((indent, node_name))
        
        return phylo_tree, node_snps
    
    def _load_isogg(self, isogg_file):
        """Load ISOGG index file"""
        count = 0
        
        with open(isogg_file, 'r', encoding='utf-8', errors='ignore') as f:
            header = f.readline()
            sep = '\t' if '\t' in header else ','
            
            for line in f:
                parts = line.strip().split(sep)
                if len(parts) < 7:
                    continue
                
                snp_name = parts[0].strip().rstrip('~^*')
                haplo = parts[1].strip().rstrip('~^*')
                pos_str = parts[4].strip()
                mutation = parts[6].strip() if len(parts) > 6 else ""
                
                if not pos_str or '..' in pos_str:
                    continue
                
                try:
                    pos = int(pos_str)
                except:
                    continue
                
                ref, alt = '', ''
                if '->' in mutation:
                    try:
                        ref, alt = mutation.split('->')[:2]
                        ref, alt = ref.strip().upper(), alt.strip().upper()
                    except:
                        pass
                
                self.snp_to_info[snp_name] = (pos, ref, alt)
                self.pos_to_haplo[pos].append((haplo, snp_name, ref, alt))
                count += 1
                
                # Process aliases
                if len(parts) > 2:
                    for alt_name in parts[2].split(';'):
                        alt_name = alt_name.strip().rstrip('~^*')
                        if alt_name and alt_name not in self.snp_to_info:
                            self.snp_to_info[alt_name] = (pos, ref, alt)
        
        print(f"    Loaded SNPs: {count}")
    
    def _build_node_snp_positions(self):
        """Build SNP position mapping for tree nodes"""
        self.node_snp_positions = {}
        
        for node, snps in self.node_snps.items():
            positions = []
            for snp in snps:
                if snp in self.snp_to_info:
                    pos, ref, alt = self.snp_to_info[snp]
                    positions.append((snp, pos, ref, alt))
            self.node_snp_positions[node] = positions
        
        nodes_with_pos = sum(1 for p in self.node_snp_positions.values() if p)
        print(f"    {nodes_with_pos}/{len(self.node_snps)} tree nodes have SNP positions")
    
    def _check_main_branch_snps(self):
        """Check availability of main branch defining SNPs"""
        print(f"\n    Main branch SNP check:")
        
        for branch in self.MAIN_BRANCHES:
            snp = self.MAIN_BRANCH_DEFINING_SNPS.get(branch, '')
            if snp and snp in self.snp_to_info:
                pos = self.snp_to_info[snp][0]
                print(f"      {branch}: ✓ {snp} @ {pos}")
            else:
                print(f"      {branch}: ✗ {snp} (Needs inference from downstream)")
    
    def get_ancestors(self, node):
        """Get all ancestors of node (nearest to farthest)"""
        ancestors = []
        current = node
        seen = set()
        
        while current in self.phylo_tree and current not in seen:
            seen.add(current)
            parent = self.phylo_tree[current]
            ancestors.append(parent)
            current = parent
        
        return ancestors
    
    def get_depth(self, node):
        """Get node depth in tree"""
        return len(self.get_ancestors(node))
    
    def get_main_branch(self, haplo):
        """Extract main branch from haplogroup name"""
        if not haplo:
            return None
        
        # Direct match
        for branch in self.MAIN_BRANCHES:
            if haplo == branch:
                return branch
            if haplo.startswith(branch) and len(haplo) > 1:
                next_char = haplo[len(branch)] if len(haplo) > len(branch) else ''
                if next_char.isdigit() or next_char == '':
                    return branch
        
        # Check if intermediate node
        if haplo in self.INTERMEDIATE_NODES:
            return None
        
        # First letter match
        first = haplo[0] if haplo else None
        if first in self.MAIN_BRANCHES:
            return first
        
        return None
    
    def check_genotype(self, gt, vcf_ref, vcf_alt, exp_ref, exp_alt):
        """
        Check genotype status
        
        Returns: GenotypeStatus
        """
        if gt in ['./.', '.', '']:
            return GenotypeStatus.MISSING
        
        vcf_ref = vcf_ref.upper()
        vcf_alt = vcf_alt.upper() if vcf_alt != '.' else ''
        
        is_het = gt in ['0/1', '0|1', '1/0', '1|0']
        is_hom_alt = gt in ['1/1', '1|1', '1']
        is_hom_ref = gt in ['0/0', '0|0', '0']
        
        # Has explicit ref/alt information
        if exp_ref and exp_alt:
            # Forward match: VCF REF=ancestral, ALT=derived
            if vcf_ref == exp_ref.upper() and vcf_alt == exp_alt.upper():
                if is_hom_alt:
                    return GenotypeStatus.DERIVED
                elif is_het:
                    return GenotypeStatus.HETEROZYGOUS
                elif is_hom_ref:
                    return GenotypeStatus.ANCESTRAL
            # Reverse match: VCF REF=derived, ALT=ancestral
            elif vcf_ref == exp_alt.upper() and vcf_alt == exp_ref.upper():
                if is_hom_ref:
                    return GenotypeStatus.DERIVED
                elif is_het:
                    return GenotypeStatus.HETEROZYGOUS
                elif is_hom_alt:
                    return GenotypeStatus.ANCESTRAL
        
        # When no explicit info, infer from VCF
        if is_hom_alt:
            return GenotypeStatus.DERIVED
        elif is_het:
            return GenotypeStatus.HETEROZYGOUS
        elif is_hom_ref:
            return GenotypeStatus.ANCESTRAL
        
        return GenotypeStatus.MISMATCH
    
    def is_derived(self, status):
        """
        Determine if treated as derived based on het_mode
        
        Returns: True/False/None (NoneIndicates missing)
        """
        if status == GenotypeStatus.DERIVED:
            return True
        elif status == GenotypeStatus.ANCESTRAL:
            return False
        elif status == GenotypeStatus.HETEROZYGOUS:
            if self.het_mode == "strict":
                return None
            else:  # moderate or lenient
                return True
        else:  # MISSING or MISMATCH
            return None
    
    def check_node_status(self, node, sample_geno):
        """
        Check status of single node
        
        Returns: {
            'status': 'derived'/'ancestral'/'missing'/'mixed',
            'derived_snps': [...],
            'ancestral_snps': [...],
            'het_snps': [...],
            'missing_snps': [...]
        }
        """
        positions = self.node_snp_positions.get(node, [])
        
        derived_snps = []
        ancestral_snps = []
        het_snps = []
        missing_snps = []
        
        for snp, pos, ref, alt in positions:
            if pos not in sample_geno:
                missing_snps.append(snp)
                continue
            
            gt, vcf_ref, vcf_alt = sample_geno[pos]
            status = self.check_genotype(gt, vcf_ref, vcf_alt, ref, alt)
            
            if status == GenotypeStatus.DERIVED:
                derived_snps.append(snp)
            elif status == GenotypeStatus.ANCESTRAL:
                ancestral_snps.append(snp)
            elif status == GenotypeStatus.HETEROZYGOUS:
                het_snps.append(snp)
            else:
                missing_snps.append(snp)
        
        # Overall determination
        if derived_snps:
            overall = 'derived'
        elif ancestral_snps and not derived_snps:
            overall = 'ancestral'
        elif het_snps and not derived_snps and not ancestral_snps:
            overall = 'het'
        else:
            overall = 'missing'
        
        return {
            'status': overall,
            'derived_snps': derived_snps,
            'ancestral_snps': ancestral_snps,
            'het_snps': het_snps,
            'missing_snps': missing_snps
        }
    
    def classify_sample(self, sample, sample_geno):
        """
        Classify single sample
        
        Core logic：
        1. Collect all derived haplogroups
        2. Group by main branch
        3. Validate ancestor chain for each candidate
        4. Select deepest validated node
        """
        het_count = 0
        
        # ========================================
        # Step 1: Collect all derived SNPs/haplogroups
        # ========================================
        derived_haplos = defaultdict(list)  # haplo -> [(snp, pos), ...]
        all_node_status = {}  # Record status of all checked nodes
        
        for pos, entries in self.pos_to_haplo.items():
            if pos not in sample_geno:
                continue
            
            gt, vcf_ref, vcf_alt = sample_geno[pos]
            
            # Count heterozygous
            if gt in ['0/1', '0|1', '1/0', '1|0']:
                het_count += 1
            
            for haplo, snp_name, ref, alt in entries:
                status = self.check_genotype(gt, vcf_ref, vcf_alt, ref, alt)
                
                # Record node status
                if haplo not in all_node_status:
                    all_node_status[haplo] = {
                        'derived': [], 'ancestral': [], 'het': [], 'missing': []
                    }
                
                if status == GenotypeStatus.DERIVED:
                    all_node_status[haplo]['derived'].append(snp_name)
                elif status == GenotypeStatus.ANCESTRAL:
                    all_node_status[haplo]['ancestral'].append(snp_name)
                elif status == GenotypeStatus.HETEROZYGOUS:
                    all_node_status[haplo]['het'].append(snp_name)
                
                # Determine if derived
                if self.is_derived(status):
                    derived_haplos[haplo].append((snp_name, pos))
        
        if not derived_haplos:
            return self._make_result(
                sample, 'Unknown', 'Unknown', 0, 0, het_count, 
                [], 'No derived SNPs', {}
            )
        
        # ========================================
        # Step 2: Group by main branch
        # ========================================
        branch_candidates = defaultdict(list)  # branch -> [(haplo, n_snps, snps), ...]
        
        for haplo, snps in derived_haplos.items():
            branch = self.get_main_branch(haplo)
            if branch:
                branch_candidates[branch].append({
                    'haplo': haplo,
                    'n_snps': len(snps),
                    'snps': [s[0] for s in snps],
                    'depth': self.get_depth(haplo)
                })
        
        if not branch_candidates:
            return self._make_result(
                sample, 'Unknown', 'Unknown', 0, 0, het_count,
                [], 'Cannot determine main branch', {}
            )
        
        # ========================================
        # Step 3: Determine main branch
        # ========================================
        # Check main branch defining SNPs
        branch_direct_status = {}
        for branch in self.MAIN_BRANCHES:
            snp = self.MAIN_BRANCH_DEFINING_SNPS.get(branch)
            if snp and snp in self.snp_to_info:
                pos, ref, alt = self.snp_to_info[snp]
                if pos in sample_geno:
                    gt, vcf_ref, vcf_alt = sample_geno[pos]
                    status = self.check_genotype(gt, vcf_ref, vcf_alt, ref, alt)
                    branch_direct_status[branch] = status
                else:
                    branch_direct_status[branch] = GenotypeStatus.MISSING
            else:
                branch_direct_status[branch] = GenotypeStatus.MISSING
        
        # Select main branch: prioritize defining SNPs, then downstream evidence
        main_branch = None
        main_branch_note = ""
        
        # First find those with derived defining SNPs
        for branch in self.MAIN_BRANCHES:
            status = branch_direct_status.get(branch)
            if self.is_derived(status) and branch in branch_candidates:
                main_branch = branch
                main_branch_note = f"{self.MAIN_BRANCH_DEFINING_SNPS[branch]} derived"
                break
        
        # Not found, infer from downstream evidence
        if not main_branch:
            # Exclude branches with clear ancestral defining SNPs
            excluded = set()
            for branch, status in branch_direct_status.items():
                if status == GenotypeStatus.ANCESTRAL:
                    excluded.add(branch)
            
            # Select one with most downstream evidence
            best_branch = None
            best_count = 0
            
            for branch, candidates in branch_candidates.items():
                if branch in excluded:
                    continue
                total_snps = sum(c['n_snps'] for c in candidates)
                if total_snps > best_count:
                    best_count = total_snps
                    best_branch = branch
            
            if best_branch:
                main_branch = best_branch
                main_branch_note = f"inferred from {best_count} downstream SNPs"
        
        # Special handling: detect basal F
        # If no main branch found but M89 derived, may be rare basal F
        if not main_branch:
            if self.F_DEFINING_SNP in self.snp_to_info:
                f_pos = self.snp_to_info[self.F_DEFINING_SNP][0]
                if f_pos in sample_geno:
                    gt, vcf_ref, vcf_alt = sample_geno[f_pos]
                    f_ref, f_alt = self.snp_to_info[self.F_DEFINING_SNP][1], self.snp_to_info[self.F_DEFINING_SNP][2]
                    f_status = self.check_genotype(gt, vcf_ref, vcf_alt, f_ref, f_alt)
                    if self.is_derived(f_status):
                        # M89Derived but no downstream branch evidence = basal F
                        main_branch = 'F'
                        main_branch_note = "basal F (M89 derived, no downstream)"
                        # Add F to branch_candidates for subsequent processing
                        if 'F' not in branch_candidates:
                            branch_candidates['F'] = [{'haplo': 'F', 'n_snps': 1, 'snps': ['M89'], 'depth': 0}]
        
        if not main_branch:
            return self._make_result(
                sample, 'Unknown', 'Unknown', 0, 0, het_count,
                [], 'Cannot determine main branch（Insufficient or conflicting evidence）', {}
            )
        
        # ========================================
        # Step 4: Subdivide within main branch with ancestor chain validation
        # ========================================
        candidates = branch_candidates[main_branch]
        
        # Sort by depth and SNP count
        candidates.sort(key=lambda x: (x['depth'], x['n_snps']), reverse=True)
        
        # Validate ancestor chain for each candidate
        best_haplo = main_branch
        best_n_snps = 0
        best_evidence = []
        ancestor_conflicts = []
        
        for cand in candidates:
            haplo = cand['haplo']
            ancestors = self.get_ancestors(haplo)
            
            # Only check ancestors within main branch
            relevant_ancestors = []
            for anc in ancestors:
                if anc == main_branch:
                    break
                if self.get_main_branch(anc) == main_branch or anc == main_branch:
                    relevant_ancestors.append(anc)
            
            # Check ancestor chain
            chain_valid = True
            conflicts = []
            
            for anc in relevant_ancestors:
                if anc in all_node_status:
                    anc_info = all_node_status[anc]
                    has_derived = len(anc_info['derived']) > 0
                    has_ancestral = len(anc_info['ancestral']) > 0
                    has_het = len(anc_info['het']) > 0
                    
                    # Key: if ancestor is clearly ancestral with no derived evidence
                    if has_ancestral and not has_derived:
                        # Can heterozygous save it?
                        if has_het and self.het_mode in ['moderate', 'lenient']:
                            pass  # Acceptable but marked
                        else:
                            chain_valid = False
                            conflicts.append(f"{anc}(ancestral)")
            
            if chain_valid:
                best_haplo = haplo
                best_n_snps = cand['n_snps']
                best_evidence = cand['snps'][:5]
                break  # Found first deepest validated node
            else:
                ancestor_conflicts.extend(conflicts)
        
        # If all candidates have conflicts，Select one with fewest conflicts
        if best_haplo == main_branch and candidates:
            # Re-evaluate, select one with fewest conflicts
            min_conflicts = float('inf')
            for cand in candidates:
                haplo = cand['haplo']
                ancestors = self.get_ancestors(haplo)
                
                conflict_count = 0
                for anc in ancestors:
                    if anc in all_node_status:
                        anc_info = all_node_status[anc]
                        if anc_info['ancestral'] and not anc_info['derived']:
                            conflict_count += 1
                
                if conflict_count < min_conflicts:
                    min_conflicts = conflict_count
                    best_haplo = haplo
                    best_n_snps = cand['n_snps']
                    best_evidence = cand['snps'][:5]
        
        # ========================================
        # Step 5: Calculate confidence
        # ========================================
        direct_status = branch_direct_status.get(main_branch, GenotypeStatus.MISSING)
        
        # Base confidence
        if self.is_derived(direct_status):
            # Main branch defining SNP directly derived
            if best_n_snps >= 5:
                confidence = 0.95
            elif best_n_snps >= 3:
                confidence = 0.90
            else:
                confidence = 0.85
        else:
            # Inferred from downstream
            if best_n_snps >= 10:
                confidence = 0.85
            elif best_n_snps >= 5:
                confidence = 0.75
            elif best_n_snps >= 3:
                confidence = 0.65
            else:
                confidence = 0.50
        
        # Penalty factors
        if ancestor_conflicts:
            confidence -= 0.15
        if het_count > 10:
            confidence -= 0.05
        
        confidence = max(0.0, min(1.0, confidence))
        
        # Build diagnostic information
        diagnostics = {
            'main_branch_status': str(direct_status),
            'main_branch_note': main_branch_note,
            'ancestor_conflicts': ancestor_conflicts[:3] if ancestor_conflicts else [],
            'total_derived_haplos': len(derived_haplos),
            'candidates_in_branch': len(candidates)
        }
        
        return self._make_result(
            sample, main_branch, best_haplo, best_n_snps, 
            confidence, het_count, best_evidence, '', diagnostics
        )
    
    def _make_result(self, sample, main_branch, haplo, n_snps, 
                     confidence, het_count, evidence, note, diagnostics):
        """Build result dictionary"""
        return {
            'sample': sample,
            'main_branch': main_branch,
            'haplogroup': haplo,
            'n_snps': n_snps,
            'confidence': confidence,
            'het_count': het_count,
            'evidence': evidence,
            'note': note,
            'diagnostics': diagnostics
        }


def load_vcf(vcf_file):
    """Load VCF file"""
    opener = open
    mode = 'r'
    if vcf_file.endswith('.gz'):
        opener = gzip.open
        mode = 'rt'
    
    samples = []
    sample_geno = defaultdict(dict)
    y_count = 0
    
    with opener(vcf_file, mode) as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                parts = line.strip().split('\t')
                samples = parts[9:]
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 10:
                continue
            
            chrom = parts[0]
            if chrom not in ['Y', 'chrY', '24', 'chr24']:
                continue
            
            pos = int(parts[1])
            ref = parts[3]
            alt = parts[4]
            
            for i, sample in enumerate(samples):
                gt = parts[9 + i].split(':')[0]
                sample_geno[sample][pos] = (gt, ref, alt)
            
            y_count += 1
    
    return samples, sample_geno, y_count


def print_summary(results, het_mode):
    """Print summary statistics"""
    print("\n" + "=" * 70)
    print("Classification Results")
    print("=" * 70)
    
    print(f"\nHeterozygosity handling modes: {het_mode}")
    
    # Confidence distribution
    high = sum(1 for r in results if r['confidence'] >= 0.8)
    med = sum(1 for r in results if 0.5 <= r['confidence'] < 0.8)
    low = sum(1 for r in results if r['confidence'] < 0.5)
    
    print(f"\n【Confidence】")
    print(f"   High(>=0.8): {high} ({high/len(results)*100:.1f}%)")
    print(f"   Medium(0.5-0.8): {med} ({med/len(results)*100:.1f}%)")
    print(f"   Low(<0.5): {low} ({low/len(results)*100:.1f}%)")
    
    # Main branch distribution
    main_counts = defaultdict(int)
    for r in results:
        main_counts[r['main_branch']] += 1
    
    print(f"\n【Main branch】")
    for branch in sorted(main_counts.keys(), key=lambda x: -main_counts[x]):
        count = main_counts[branch]
        print(f"   {branch}: {count} ({count/len(results)*100:.1f}%)")
    
    # Detailed haplogroups
    sub_counts = defaultdict(int)
    for r in results:
        sub_counts[r['haplogroup']] += 1
    
    print(f"\n【Detailed haplogroups (top 15)】")
    for haplo in sorted(sub_counts.keys(), key=lambda x: -sub_counts[x])[:15]:
        count = sub_counts[haplo]
        print(f"   {haplo:30s}: {count:4d} ({count/len(results)*100:.1f}%)")
    
    # Ancestor chain conflict statistics
    conflict_count = sum(1 for r in results 
                         if r.get('diagnostics', {}).get('ancestor_conflicts'))
    if conflict_count > 0:
        print(f"\n【Warning】")
        print(f"   {conflict_count} samples have ancestor chain conflicts (confidence reduced)")


def save_results(results, output_dir, het_mode):
    """Save results"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Main result file
    result_file = os.path.join(output_dir, 'Y_haplogroup_result_v4.0.txt')
    with open(result_file, 'w') as f:
        f.write("Sample\tMain_Branch\tHaplogroup\tConfidence\tN_SNPs\t")
        f.write("Het_Count\tMain_Status\tEvidence\tConflicts\n")
        
        for r in results:
            evidence = ','.join(r['evidence'][:5]) if r['evidence'] else '-'
            diag = r.get('diagnostics', {})
            main_status = diag.get('main_branch_status', '-')
            conflicts = ','.join(diag.get('ancestor_conflicts', [])) or '-'
            
            f.write(f"{r['sample']}\t{r['main_branch']}\t{r['haplogroup']}\t")
            f.write(f"{r['confidence']:.3f}\t{r['n_snps']}\t{r['het_count']}\t")
            f.write(f"{main_status}\t{evidence}\t{conflicts}\n")
    
    # Summary file
    summary_file = os.path.join(output_dir, 'Y_haplogroup_summary_v4.0.txt')
    with open(summary_file, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write(f"Y-chromosomeHaplogroup Classification Report v{VERSION}\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"Sample count: {len(results)}\n")
        f.write(f"Heterozygosity handling modes: {het_mode}\n\n")
        
        # Confidence
        high = sum(1 for r in results if r['confidence'] >= 0.8)
        med = sum(1 for r in results if 0.5 <= r['confidence'] < 0.8)
        low = sum(1 for r in results if r['confidence'] < 0.5)
        
        f.write("Confidence distribution:\n")
        f.write(f"  High(>=0.8): {high} ({high/len(results)*100:.1f}%)\n")
        f.write(f"  Medium(0.5-0.8): {med} ({med/len(results)*100:.1f}%)\n")
        f.write(f"  Low(<0.5): {low} ({low/len(results)*100:.1f}%)\n\n")
        
        # Main branch
        main_counts = defaultdict(int)
        for r in results:
            main_counts[r['main_branch']] += 1
        
        f.write("Main branch distribution:\n")
        for branch in sorted(main_counts.keys(), key=lambda x: -main_counts[x]):
            count = main_counts[branch]
            f.write(f"  {branch}: {count} ({count/len(results)*100:.1f}%)\n")
        
        # Subdivide
        sub_counts = defaultdict(int)
        for r in results:
            sub_counts[r['haplogroup']] += 1
        
        f.write("\nDetailed haplogroups:\n")
        for haplo in sorted(sub_counts.keys(), key=lambda x: -sub_counts[x])[:30]:
            count = sub_counts[haplo]
            f.write(f"  {haplo}: {count} ({count/len(results)*100:.1f}%)\n")
    
    print(f"\n[5] Results saved:")
    print(f"    {result_file}")
    print(f"    {summary_file}")
    
    return result_file, summary_file


def main():
    parser = argparse.ArgumentParser(
        description=f'Y-chromosome Haplogroup Classifier v{VERSION}',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Classification logic description:
  This tool simulates manual classification, traversing phylogenetic tree from root:
  1. Check defining SNP status at each node
  2. Derived → enter node, continue downward
  3. Ancestral → stop, do not enter branch  
  4. Missing → Check downstream forDerivedEvidence（Can be inferred）
  5. Validate ancestor chain consistency, exclude conflicting paths

Heterozygosity handling modes:
  strict:   Heterozygous treated as missing (most conservative)
  moderate: Heterozygous treated as derived but marked (default)
  lenient:  Heterozygous fully treated as derived
        '''
    )
    
    parser.add_argument('-v', '--vcf', required=True, help='VCF file')
    parser.add_argument('-i', '--isogg', required=True, help='ISOGG indexdata.csv')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('--het-mode', choices=['strict', 'moderate', 'lenient'],
                        default='moderate', help='Heterozygosity handling mode (Default: moderate)')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.vcf):
        print(f"Error: {args.vcf} Does not exist")
        sys.exit(1)
    if not os.path.exists(args.isogg):
        print(f"Error: {args.isogg} Does not exist")
        sys.exit(1)
    
    print("=" * 70)
    print(f"Y-chromosome Haplogroup Classifier v{VERSION}")
    print("=" * 70)
    
    # Initialize classifier
    classifier = HaplogroupClassifier(args.isogg, args.het_mode)
    
    # LoadVCF
    print(f"\n[3] Reading VCF: {args.vcf}")
    samples, sample_geno, y_count = load_vcf(args.vcf)
    print(f"    Sample count: {len(samples)}")
    print(f"    Y-chromosome sites: {y_count}")
    
    if not samples:
        print("Error: No samples")
        sys.exit(1)
    
    # Classification
    print(f"\n[4] Classifying... (Heterozygosity mode: {args.het_mode})")
    results = []
    
    for i, sample in enumerate(samples):
        if (i + 1) % 50 == 0 or i == len(samples) - 1:
            print(f"    {i + 1}/{len(samples)}", end='\r')
        
        result = classifier.classify_sample(sample, sample_geno[sample])
        results.append(result)
    
    print("\n    Complete!♡(*´∀｀*)人(*´∀｀*)♡")
    
    # Output
    print_summary(results, args.het_mode)
    save_results(results, args.output, args.het_mode)
    
    print("\n" + "=" * 70)
    print("Complete!♡(*´∀｀*)人(*´∀｀*)♡")
    print("=" * 70)


if __name__ == "__main__":
    main()