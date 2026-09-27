import os, sys, subprocess

base_dir = os.path.dirname(os.path.abspath(__file__))

tasks = [
    ("Task 1: Official Model Authenticity Audit", "task1_audit_official_models.py"),
    ("Task 2: Strict Scaffold x Lineage OOD", "task2_strict_scaffold_lineage.py"),
    ("Task 3: TOST Equivalence Testing", "task3_tost_equivalence.py"),
    ("Task 4: Training-Only Stress Decomposition", "task4_training_only_stress.py"),
    ("Task 5: Multi-Database Pathway Replication", "task5_multidatabase_pathway.py"),
    ("Task 6: ProteinTalks Matching Audit", "task6_proteintalks_audit.py"),
    ("Task 7: PRISM Directionality Unit Test", "prism_direction_unit_test.py"),
    ("Task 8: GDSC / CTRP External Replication", "task8_gdsc_ctrp_replication.py"),
    ("Task 9: d_joint Independent Validation", "task9_djoint_validation.py")
]

print("=" * 80)
print("RUNNING MASTER FINAL-ROUND DEEP-DIVE VALIDATION SUITE")
print("=" * 80)

for name, script in tasks:
    script_path = os.path.normpath(os.path.join(base_dir, script))
    print(f"\n>>> Running {name} ({os.path.basename(script_path)})...", flush=True)
    res = subprocess.run([sys.executable, script_path], capture_output=True, text=True, cwd=base_dir)
    if res.returncode != 0:
        print(f"FAILED {name}:")
        print(res.stderr)
        sys.exit(1)
    else:
        print(res.stdout.strip())

print("\n" + "=" * 80)
print("SUCCESS: ALL 9 FINAL-ROUND VALIDATION TASKS COMPLETED AND VERIFIED!")
print("=" * 80)
