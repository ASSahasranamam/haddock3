# Issues to Fix Before Resubmitting PR

## 1. CI FAILURE: `test_alascan_with_ligand_topar` fails on Python 3.14 + ARM

Root cause: **NOT duplicate DIHEdral** (those 3 lines at 1232-1234 are original upstream
lines, just shifted down by 3 due to our earlier insertions — they are NOT duplicates).

Still investigating. The alascan test calls CNS with the modified param file. Something
in our param additions is causing CNS to fail silently on Python 3.14/ARM, leaving no
output TSV. Candidates:

### a) IMPRoper at line 813 (added by us, unique):
```
IMPRoper  C    CH2E NH2  O   500.00 {sd=0.031}  0  0.0124  ! AHX->CTN planar
```
Upstream has NO equivalent. Could this be conflicting with existing C/CH2E interactions?

### b) IMPRoper at line 855 (added by us, unique):
```
IMPRoper  CH0  C    NH1  CH2E  500.00 {sd=0.031}  0  -180.0067  ! AHX after IVA
```
- `-180.0067` equilibrium value — likely a numeric artifact, should probably be `-180.0000`
- Verify this entry is actually needed and correct

### c) The diff says these lines were added near the top of the DIHEdral section
(before lines 1232-1234). Check if inserting parameters mid-block causes CNS
to mis-parse the rest of the block on some platforms.

---

## 2. `covalions.cns` — Prof. Bonvin request ✓ DONE
Added IVA, LZH, ADI, AHX, VHP to:
`src/haddock/modules/refinement/flexref/cns/covalions.cns`  
Committed: 7e4813cdd

---

## 3. Codacy coverage upload failure
Not fixable from PR — maintainers need their Codacy secret.

---

## Next steps
1. Identify actual cause of alascan failure (run test with verbose CNS output)
2. Fix param file issue
3. Push and re-run CI
4. If clean → close #1556 and open new PR
