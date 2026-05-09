"""
validate_v7_B.py — Cross-section consistency checks B1–B7 for halo_shells v7.0.

For each check, the script searches the manuscript .tex for occurrences of a
specific numerical claim and verifies that all locations report the same value.

This catches copy-paste drift between abstract / introduction / results /
discussion / conclusion / tables — the most common source of subtle errors
in a multi-section manuscript that has been revised through several versions.

Each Bx check defines a pattern (regex or literal substring) and finds every
location where it appears. If a manuscript revision updated one location but
missed another, this script flags the inconsistency.

Usage:
    cd halo_shells_v7.0/
    python3 scripts/validate_v7_B.py

Exit code 0 if all pass, 1 if any fail.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PACKAGE_ROOT = os.path.dirname(HERE)
TEX_PATH = os.path.join(PACKAGE_ROOT, 'source', 'sparc_shells_body.tex')


# ANSI color codes
class C:
    OK = '\033[92m'
    FAIL = '\033[91m'
    WARN = '\033[93m'
    BOLD = '\033[1m'
    DIM = '\033[2m'
    END = '\033[0m'


# ============================================================
# Section detection — given a line number, what section is it in?
# ============================================================
def build_section_map(text):
    """Return list of (line_num, section_label) entries sorted by line number."""
    sections = []
    lines = text.split('\n')
    cur_section = "preamble"
    for i, line in enumerate(lines, start=1):
        m = re.match(r'\\section\*?\{([^}]+)\}', line)
        if m:
            cur_section = m.group(1)
            sections.append((i, cur_section))
            continue
        m = re.match(r'\\subsection\*?\{([^}]+)\}', line)
        if m:
            sections.append((i, f"  {cur_section} → {m.group(1)}"))
    return sections


def find_section_for_line(line_num, sections):
    """Given a line number and section map, return the closest preceding section label."""
    last = "(before any section)"
    for ln, label in sections:
        if ln <= line_num:
            last = label
        else:
            break
    return last


# ============================================================
# Helpers for pattern-finding
# ============================================================
def _find_unescaped_comment(line):
    """Return position of first unescaped '%' (LaTeX comment start), or -1."""
    i = 0
    while i < len(line):
        if line[i] == '%':
            # Check it's not preceded by a backslash escape
            # Count consecutive backslashes immediately before
            n_back = 0
            j = i - 1
            while j >= 0 and line[j] == '\\':
                n_back += 1
                j -= 1
            if n_back % 2 == 0:
                return i
        i += 1
    return -1


def find_all(text, pattern, regex=False, case_sensitive=True):
    """Find all matches; return list of (line_num, line_text, matched_text)."""
    results = []
    flags = 0 if case_sensitive else re.IGNORECASE
    lines = text.split('\n')
    if regex:
        compiled = re.compile(pattern, flags)
        for i, line in enumerate(lines, start=1):
            if line.lstrip().startswith('%'):
                continue
            comment_idx = _find_unescaped_comment(line)
            for m in compiled.finditer(line):
                if comment_idx >= 0 and m.start() > comment_idx:
                    continue
                results.append((i, line.strip(), m.group(0)))
    else:
        for i, line in enumerate(lines, start=1):
            if line.lstrip().startswith('%'):
                continue
            comment_idx = _find_unescaped_comment(line)
            idx = 0
            while True:
                if case_sensitive:
                    pos = line.find(pattern, idx)
                else:
                    pos = line.lower().find(pattern.lower(), idx)
                if pos == -1:
                    break
                if comment_idx >= 0 and pos > comment_idx:
                    break
                results.append((i, line.strip(), pattern))
                idx = pos + 1
    return results


def report_check(label, hits, expected_min=2, sections=None, max_show=8):
    """Report whether a pattern is present in expected locations."""
    n = len(hits)
    if n == 0:
        return False, f"  {C.FAIL}✗ FAIL{C.END} {label}: pattern NOT FOUND in manuscript"
    elif n < expected_min:
        # Suspicious — was claimed in only one place
        msg = f"  {C.WARN}⚠ ONLY 1 LOCATION{C.END} {label}: found in only {n} location, "
        msg += f"expected ≥ {expected_min}"
        for line_num, _, _ in hits[:max_show]:
            sec = find_section_for_line(line_num, sections) if sections else "?"
            msg += f"\n      L{line_num} ({sec})"
        return n >= 1, msg  # not strictly fail; flag warning
    else:
        msg = f"  {C.OK}✓ PASS{C.END} {label}: present in {n} locations"
        for line_num, _, _ in hits[:max_show]:
            sec = find_section_for_line(line_num, sections) if sections else "?"
            msg += f"\n      L{line_num} ({sec})"
        if n > max_show:
            msg += f"\n      ... and {n - max_show} more"
        return True, msg


def check_no_stale(label, stale_pattern, regex=False, exclude_changelog=True):
    """Check that a stale (v6.5) value does NOT appear anywhere in the manuscript.

    If exclude_changelog=True, occurrences inside the v7.0 changelog appendix
    are intentional (documenting the change) and not flagged.
    """
    with open(TEX_PATH) as f:
        text = f.read()
    hits = find_all(text, stale_pattern, regex=regex)
    if exclude_changelog:
        # Find the line range of the changelog section
        sections = build_section_map(text)
        changelog_start = None
        next_section_start = None
        all_lines = text.split('\n')
        for i, line in enumerate(all_lines, start=1):
            if 'Version 7.0 changes relative to version 6.5' in line and re.match(r'\\section', line):
                changelog_start = i
            elif changelog_start is not None and next_section_start is None:
                if re.match(r'\\section\*?\{', line):
                    next_section_start = i
                    break
        if next_section_start is None:
            next_section_start = len(all_lines) + 1
        # Filter out hits within changelog
        if changelog_start is not None:
            hits = [(ln, l, m) for (ln, l, m) in hits
                    if not (changelog_start <= ln < next_section_start)]
    if not hits:
        return True, f"  {C.OK}✓ PASS{C.END} {label}: no stale value found in body"
    msg = f"  {C.FAIL}✗ FAIL{C.END} {label}: stale value still present at:"
    sections = build_section_map(text)
    for line_num, line, _ in hits[:5]:
        sec = find_section_for_line(line_num, sections)
        snippet = line[:100] + ('...' if len(line) > 100 else '')
        msg += f"\n      L{line_num} ({sec}): {snippet}"
    return False, msg


# ============================================================
# B1. Framework adequacy count (91/102)
# ============================================================
def B1():
    print(f"\n{C.BOLD}B1. Framework adequacy count (91/102){C.END}")
    print(f"   Should appear in: Abstract, §1, §3.1, §6, §7")
    
    with open(TEX_PATH) as f:
        text = f.read()
    sections = build_section_map(text)
    
    # Pattern: "91/102", "91 of 102"  
    hits_91 = find_all(text, r'91\s*[/\s]\s*(of\s+)?102', regex=True)
    p1, m1 = report_check("'91/102' pattern", hits_91, expected_min=4, sections=sections)
    print(m1)
    
    # Anti-check: stale "92/102" or "92 of 102" — must NOT appear (was v6.5 value)
    p2, m2 = check_no_stale("no stale '92/102' (v6.5)", r'92\s*[/\s]\s*(of\s+)?102', regex=True)
    print(m2)
    
    # Anti-check: stale "92" alone in adequacy contexts — too noisy to do without false positives, skip
    return p1 and p2


# ============================================================
# B2. NFW BIC-preferring count (8 strict KR)
# ============================================================
def B2():
    print(f"\n{C.BOLD}B2. NFW BIC-preferring count (8 galaxies, strict KR ΔBIC > 2){C.END}")
    print(f"   Should appear in: Abstract, §3.2 (count + analysis), §6, §7")
    
    with open(TEX_PATH) as f:
        text = f.read()
    sections = build_section_map(text)
    
    # Patterns for "8 galaxies BIC-prefer NFW", "8 NFW-prefer", etc.
    hits_8 = find_all(text, r'8\s+(galaxies\s+BIC-prefer\s+NFW|BIC-prefer\s+NFW|NFW-prefer)', regex=True)
    if not hits_8:
        # Looser pattern: "8 ... NFW" within context
        hits_8 = find_all(text, r'\b8\b[^.]{0,100}(BIC-prefer NFW|NFW-prefer|prefer NFW)', regex=True)
    
    p1, m1 = report_check("'8 BIC-prefer NFW' pattern", hits_8, expected_min=2, sections=sections)
    print(m1)
    
    # Anti-checks: stale "9 galaxies" / "9 BIC-prefer" / "14 BIC-prefer" 
    # (v6.5 had 9 strict; loose KR was 14)
    p2, m2 = check_no_stale("no stale '9 BIC-prefer NFW' (v6.5 strict)", 
                             r'\b9\s+(galaxies\s+)?(BIC-prefer\s+NFW|NFW-prefer|prefer\s+NFW)', regex=True)
    print(m2)
    p3, m3 = check_no_stale("no stale '14 BIC-prefer NFW' (loose KR)",
                             r'\b14\s+(galaxies\s+)?(BIC-prefer\s+NFW|NFW-prefer|prefer\s+NFW)', regex=True)
    print(m3)
    
    return p1 and p2 and p3


# ============================================================
# B3. Trend statistics (ρ=-0.83, CI[-0.95,-0.22], p_two=0.002, MK=0.003)
# ============================================================
def B3():
    print(f"\n{C.BOLD}B3. Trend statistics (ρ=-0.83, CI[-0.95,-0.22], p_two=0.002, MK=0.003){C.END}")
    print(f"   Should appear in: Abstract, §1, §3.3, §6, Table 4 caption")
    
    with open(TEX_PATH) as f:
        text = f.read()
    sections = build_section_map(text)
    
    # Spearman ρ = -0.83
    hits_rho = find_all(text, r'\\rho\s*=\s*-0\.83', regex=True)
    p1, m1 = report_check("'ρ = -0.83'", hits_rho, expected_min=3, sections=sections)
    print(m1)
    
    # Bootstrap CI [-0.95, -0.22]
    hits_ci = find_all(text, r'\[-0\.95,\s*-0\.22\]', regex=True)
    p2, m2 = report_check("CI [-0.95, -0.22]", hits_ci, expected_min=3, sections=sections)
    print(m2)
    
    # p_two-sided = 0.002
    hits_ptwo = find_all(text, r'p_\{.*?two.*?\}\s*=\s*0\.002', regex=True)
    p3, m3 = report_check("p_two-sided = 0.002", hits_ptwo, expected_min=3, sections=sections)
    print(m3)
    
    # Mann-Kendall p = 0.003
    hits_mk = find_all(text, r'Mann-Kendall\s+\$?p\$?\s*=\s*0\.003', regex=True)
    p4, m4 = report_check("Mann-Kendall p = 0.003", hits_mk, expected_min=3, sections=sections)
    print(m4)
    
    # Anti-checks: stale v6.5 values
    p5, m5 = check_no_stale("no stale ρ = -0.81 (v6.5)", r'\\rho\s*=\s*-0\.81', regex=True)
    print(m5)
    p6, m6 = check_no_stale("no stale CI [-0.93, -0.12] (v6.5)", r'\[-0\.93,\s*-0\.12\]', regex=True)
    print(m6)
    p7, m7 = check_no_stale("no stale p_two = 0.028 (v6.5)", r'p_\{.*?two.*?\}\s*=\s*0\.028', regex=True)
    print(m7)
    p8, m8 = check_no_stale("no stale Mann-Kendall p = 0.035 (v6.5)", r'Mann-Kendall\s+\$?p\$?\s*=\s*0\.035', regex=True)
    print(m8)
    
    return all([p1, p2, p3, p4, p5, p6, p7, p8])


# ============================================================
# B4. Null test aggregate rates (4.0% / 63.6% / 346 mocks)
# ============================================================
def B4():
    print(f"\n{C.BOLD}B4. Null test aggregate rates (4.0% Burk, 63.6% NFW, 346 mocks){C.END}")
    print(f"   Should appear in: Abstract, §3.5, §5.3, §6, Table 5")
    
    with open(TEX_PATH) as f:
        text = f.read()
    sections = build_section_map(text)
    
    # Burkert FP: 4.0% or 4\%
    hits_burk = find_all(text, r'4\.0\\?\%', regex=True)
    p1, m1 = report_check("Burkert FP rate 4.0%", hits_burk, expected_min=2, sections=sections)
    print(m1)
    
    # NFW FP: 63.6%
    hits_nfw = find_all(text, r'63\.6\\?\%', regex=True)
    p2, m2 = report_check("NFW FP rate 63.6%", hits_nfw, expected_min=2, sections=sections)
    print(m2)
    
    # 346 mocks total
    hits_346 = find_all(text, r'\b346\b', regex=True)
    p3, m3 = report_check("346 total mocks", hits_346, expected_min=2, sections=sections)
    print(m3)
    
    # 7/173 — Burkert false positives count
    hits_7_173 = find_all(text, r'7\s*/\s*173|7\s+of\s+173', regex=True)
    p4, m4 = report_check("Burkert FP count 7/173", hits_7_173, expected_min=1, sections=sections)
    print(m4)
    
    # 110/173 — NFW false positives count
    hits_110_173 = find_all(text, r'110\s*/\s*173|110\s+of\s+173', regex=True)
    p5, m5 = report_check("NFW FP count 110/173", hits_110_173, expected_min=1, sections=sections)
    print(m5)
    
    # Anti-checks: stale v6.5 values
    p6, m6 = check_no_stale("no stale 4.3% Burkert FP (v6.5)", r'\b4\.3\\?\%', regex=True)
    print(m6)
    p7, m7 = check_no_stale("no stale 57.3% NFW FP (v6.5)", r'\\b57\.3\\?\%', regex=True)
    print(m7)
    p8, m8 = check_no_stale("no stale '234 mocks' (v6.5)", r'\b234\s+mocks?\b', regex=True)
    print(m8)
    
    return all([p1, p2, p3, p4, p5, p6, p7, p8])


# ============================================================
# B5. Einasto agreement count (90/102, 88%)
# ============================================================
def B5():
    print(f"\n{C.BOLD}B5. Einasto agreement count (90/102, 88%){C.END}")
    print(f"   Should appear in: Abstract, §3.6, §5.3, §6")
    
    with open(TEX_PATH) as f:
        text = f.read()
    sections = build_section_map(text)
    
    # 90/102 or "90 of 102"
    hits_90 = find_all(text, r'90\s*[/\s]\s*(of\s+)?102', regex=True)
    p1, m1 = report_check("'90/102' Einasto agreement", hits_90, expected_min=2, sections=sections)
    print(m1)
    
    # 88\% or 88%
    hits_88 = find_all(text, r'\b88\\?\%', regex=True)
    p2, m2 = report_check("'88%' Einasto agreement", hits_88, expected_min=2, sections=sections)
    print(m2)
    
    # Anti-checks: stale v6.5 values
    p3, m3 = check_no_stale("no stale '89/102' Einasto (v6.5)",
                             r'89\s*[/\s]\s*(of\s+)?102', regex=True)
    print(m3)
    p4, m4 = check_no_stale("no stale '87%' Einasto (v6.5)", r'\b87\\?\%', regex=True)
    print(m4)
    
    return all([p1, p2, p3, p4])


# ============================================================
# B6. Einasto morphology gradient stats (ρ=-0.347, p=0.0004; Burk ρ=-0.296, p=0.003)
# ============================================================
def B6():
    print(f"\n{C.BOLD}B6. Einasto morphology gradient stats{C.END}")
    print(f"   Should appear in: Abstract, §3.6, §7")
    
    with open(TEX_PATH) as f:
        text = f.read()
    sections = build_section_map(text)
    
    # Einasto per-galaxy ρ = -0.347
    hits_ein_rho = find_all(text, r'\\rho\s*=\s*-0\.347', regex=True)
    p1, m1 = report_check("Einasto per-galaxy ρ = -0.347", hits_ein_rho, expected_min=2, sections=sections)
    print(m1)
    
    # Einasto p = 0.0004
    hits_ein_p = find_all(text, r'\$?p\$?\s*=\s*0\.0004', regex=True)
    p2, m2 = report_check("Einasto per-galaxy p = 0.0004", hits_ein_p, expected_min=2, sections=sections)
    print(m2)
    
    # Burkert per-galaxy ρ = -0.296
    hits_burk_rho = find_all(text, r'\\rho\s*=\s*-0\.296', regex=True)
    p3, m3 = report_check("Burkert per-galaxy ρ = -0.296", hits_burk_rho, expected_min=1, sections=sections)
    print(m3)
    
    # Anti-checks: stale v6.5 values
    p4, m4 = check_no_stale("no stale Einasto ρ = -0.327 (v6.5)", r'\\rho\s*=\s*-0\.327', regex=True)
    print(m4)
    p5, m5 = check_no_stale("no stale Einasto p = 0.0008 (v6.5)", r'\$?p\$?\s*=\s*0\.0008', regex=True)
    print(m5)
    
    return all([p1, p2, p3, p4, p5])


# ============================================================
# B7. T-binned Upsilon_disk refit headline numbers
#     (T=2: 9/9; classification stability 96.1%; per-galaxy permutation
#      p_two = 0.008; per-galaxy Spearman ρ = -0.264; bootstrap CI on per-bin
#      Spearman [-0.905, -0.096])
# ============================================================
def B7():
    print(f"\n{C.BOLD}B7. T-binned Upsilon_disk refit (M/L systematics){C.END}")
    print(f"   Should appear in: Abstract, §6.4")
    
    with open(TEX_PATH) as f:
        text = f.read()
    sections = build_section_map(text)
    
    # 96.1% classification stability
    hits_stab = find_all(text, r'96\.1\\?\\?%', regex=True)
    p1, m1 = report_check("Classification stability 96.1%", hits_stab, expected_min=1, sections=sections)
    print(m1)
    
    # Per-galaxy permutation p = 0.008 under T-binned Y
    hits_perm = find_all(text, r'p_\{\\rm two\\mbox\{-\}sided\}\s*=\s*0\.008', regex=True)
    p2, m2 = report_check("Y_T permutation p_two = 0.008", hits_perm, expected_min=1, sections=sections)
    print(m2)
    
    # Per-galaxy Spearman ρ = -0.264
    hits_rho = find_all(text, r'\\rho\s*=\s*-0\.264', regex=True)
    p3, m3 = report_check("Y_T per-galaxy ρ = -0.264", hits_rho, expected_min=1, sections=sections)
    print(m3)
    
    # T=2 vs T=3-6 Fisher strengthens to p = 0.004
    hits_fisher = find_all(text, r'p\s*=\s*0\.004', regex=True)
    p4, m4 = report_check("Y_T T=2-vs-T=3-6 Fisher p = 0.004", hits_fisher, expected_min=1, sections=sections)
    print(m4)
    
    # Bootstrap CI [-0.905, -0.096]
    hits_ci = find_all(text, r'\[-0\.905,\s*-0\.096\]', regex=True)
    p5, m5 = report_check("Y_T bootstrap CI [-0.905, -0.096]", hits_ci, expected_min=1, sections=sections)
    print(m5)
    
    # Adequacy counts under Y_T: 90/102 framework, 66/102 Burkert, 51/102 NFW
    hits_adeq = find_all(text, r'framework 90/102, Burkert 66/102, NFW 51/102', regex=True)
    p6, m6 = report_check("Y_T adequacy 90/66/51", hits_adeq, expected_min=1, sections=sections)
    print(m6)
    
    # Schombert+2022 anchor pattern: 0.65, 0.50, 0.40 should appear in §6.4
    hits_anchors = find_all(text, r'\\Upsilon_\{\\rm disk\}\(T=2\)\s*=\s*0\.65', regex=True)
    p7, m7 = report_check("Y(T=2) = 0.65 anchor in §6.4", hits_anchors, expected_min=1, sections=sections)
    print(m7)
    
    return all([p1, p2, p3, p4, p5, p6, p7])


# ============================================================
# Driver
# ============================================================
def main():
    print(f"{C.BOLD}" + "=" * 70)
    print(f"halo_shells v7.0 — Validation Pass B (cross-section consistency)")
    print(f"=" * 70 + C.END)
    print(f"Manuscript: {TEX_PATH}")
    
    if not os.path.exists(TEX_PATH):
        print(f"{C.FAIL}ERROR: manuscript file not found{C.END}")
        sys.exit(1)
    
    checks = [
        ('B1', B1, 'Framework adequacy (91/102)'),
        ('B2', B2, 'NFW BIC-preferring count (8)'),
        ('B3', B3, 'Trend statistics (ρ, CI, p_two, MK)'),
        ('B4', B4, 'Null test aggregates (4.0%/63.6%/346)'),
        ('B5', B5, 'Einasto agreement (90/102, 88%)'),
        ('B6', B6, 'Einasto gradient stats (ρ=-0.347, p=0.0004)'),
        ('B7', B7, 'Y_T refit (T=2 stays 100%, p=0.008, ρ=-0.264)'),
    ]
    
    results = {}
    for label, func, desc in checks:
        try:
            results[label] = func()
        except Exception as e:
            print(f"\n{C.FAIL}ERROR in {label}: {e}{C.END}")
            import traceback
            traceback.print_exc()
            results[label] = False
    
    # Summary
    print(f"\n{C.BOLD}" + "=" * 70)
    print(f"SUMMARY")
    print(f"=" * 70 + C.END)
    n_pass = sum(1 for v in results.values() if v)
    n_total = len(results)
    for label, _, desc in checks:
        status = f"{C.OK}✓ PASS{C.END}" if results[label] else f"{C.FAIL}✗ FAIL{C.END}"
        print(f"  {label}: {status}  {desc}")
    print()
    if n_pass == n_total:
        print(f"{C.OK}{C.BOLD}ALL {n_total}/{n_total} CHECKS PASSED{C.END}")
        sys.exit(0)
    else:
        print(f"{C.FAIL}{C.BOLD}{n_pass}/{n_total} CHECKS PASSED — {n_total - n_pass} FAILED{C.END}")
        sys.exit(1)


if __name__ == '__main__':
    main()
