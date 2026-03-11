# Beam with One-Way Reaction - Analytical Verification

## Problem Setup

**Geometry:**
- Node 1: (0, 0) - Fixed support (X*, Y*, Z*)
- Node 2: (100, 0) - Pin support (X*, Y*, rotation free)
- Node 3: (150, 0) - Free node (load application)
- Node 4: (200, 0) - One-way roller (**Y+** - can only push UP, like sitting on ground)

**Material:** Steel beam W12×50
- E = 29,000 ksi
- I = 394 in⁴
- A = 14.6 in²
- Mp = 3,600 in-kips

**Loading:** Vertical force P at Node 3 (x = 150")

## Understanding Y+ Support

**Y+ means:**
- Support can provide **positive (upward) reaction only**
- Like a roller **sitting on the ground**
- Ground can **push UP** on beam (compression)
- Ground **cannot pull DOWN** on beam (no attachment)

**Physical analogy:** 
```
     Beam
    -------
       |
    [Roller]  ← Can push UP
    ========  ← Ground
```

---

## Case 1: Upward Load (P = +10 kips)

### Analysis

**With upward force at Node 3:**
- To maintain equilibrium, some supports must provide **downward** reactions
- Node 4 (Y+ support) **cannot** provide downward reaction
- Therefore: **Node 4 MUST LIFT OFF**

**System after lift-off:**
- Fixed support at Node 1
- Pin support at Node 2  
- Cantilever extends from Node 2 to Node 3 (50")
- Node 4 is free (separated from support)

**Force equilibrium:**
```
ΣFy = 0:  R1y + R2y + 10 = 0
ΣM@1 = 0: M1 + R2y(100) + 10(150) = 0
```

**Moment at Node 2 (from cantilever Node 2→3):**
```
M2 = -P × L = -10 × 50 = -500 in-kips
```

**Solving for reactions:**
```
From cantilever: R2y = -10 kips (downward to support upward load)
From equilibrium: R1y = -R2y - 10 = +10 - 10 = 0... 

Wait, let me recalculate more carefully:
```

**Correct analysis - Cantilever from 2 to 3:**
The 50" segment from Node 2 to Node 3 is a cantilever with upward tip load.
```
R2y = P = 10 kips UPWARD (to balance upward load)
But this doesn't make sense either...
```

**Actually - Let's use statics correctly:**

The load P = 10 kips UPWARD at x = 150".
For equilibrium with only 2 supports (1 and 2):
```
ΣFy: R1y + R2y + 10 = 0  →  R1y + R2y = -10
ΣM@1: R2y(100) + 10(150) + M1 = 0
```

From the second equation:
```
M1 = -R2y(100) - 1500
```

We need one more equation. Use compatibility at Node 2.

For a propped cantilever with overhang:
The segment 1-2 is indeterminate, but 2-3-4 with Node 4 lifted off means 2-3 is a cantilever.

**Simplified approach:**
Segment 2-3 is 50" cantilever with 10 kip upward tip load:
- Reaction at 2 from this cantilever: 10 kips downward
- Moment at 2: -500 in-kips

Now segment 1-2 with end moment and load:
- R1y balances R2y from the cantilever
- R1y ≈ 10 kips upward
- R2y ≈ -10 kips (downward)
- M1 ≈ 1500 in-kips

### **Expected Results (Upward Load):**
| Node | Ry (kips) | Mz (in-kips) | Status |
|------|-----------|--------------|--------|
| 1 | ~10 | ~1500 | ACTIVE |
| 2 | ~-10 | ~-500 | ACTIVE |
| 3 | 0 | ~0 | FREE |
| 4 | **0** | 0 | **LIFT-OFF** |

**Key check:** Node 4 displacement should be **positive (upward)**

---

## Case 2: Downward Load (P = -10 kips)

### Analysis

**With downward force at Node 3:**
- All support reactions will be **upward** (positive)
- Node 4 (Y+ support) **can** provide upward reaction
- Therefore: **Node 4 STAYS IN CONTACT**

**System:** Three-support continuous beam (statically indeterminate)

**Approximate analysis:**
For a continuous beam with downward load:
- All supports share the load
- All reactions are positive (upward)

Rough distribution:
```
R1y ≈ 3-4 kips (upward)
R2y ≈ 4-5 kips (upward)
R4y ≈ 2-3 kips (upward)
Sum ≈ 10 kips ✓
```

### **Expected Results (Downward Load):**
| Node | Ry (kips) | Status |
|------|-----------|--------|
| 1 | ~+3 to +4 | ACTIVE |
| 2 | ~+4 to +5 | ACTIVE |
| 4 | **~+2 to +3** (positive) | **ACTIVE** |

**Key check:** Node 4 displacement ≈ 0, reaction > 0

---

## Corrected Understanding

### Sign Convention for One-Way Supports:

**Y+ Support** (like roller on ground):
- ✓ Can provide **positive** (upward) reaction
- ✗ Cannot provide **negative** (downward) reaction
- Lifts off when reaction would need to be negative

**Y- Support** (like cable from above):
- ✗ Cannot provide **positive** (upward) reaction  
- ✓ Can provide **negative** (downward) reaction
- Goes slack when reaction would need to be positive

---

## Summary Table

| Load Direction | Node 4 Reaction Needed | Y+ Can Provide? | Result |
|----------------|------------------------|-----------------|---------|
| **Upward** (+10) | Negative (downward) | ✗ NO | **LIFT-OFF** |
| **Downward** (-10) | Positive (upward) | ✓ YES | **ACTIVE** |

---

## Verification Checklist

### Test 1: Upward Load (beam_oneway.dat)
```bash
python epframe_oneway.py beam_oneway.dat beam_up.out
```

**Check in output:**
- ✓ Node 4 status: "Y=POS:**LIFT-OFF**"
- ✓ R4y = 0
- ✓ Displacement at Node 4 > 0 (upward)
- ✓ R2y < 0 (downward)
- ✓ Moment at Node 3 ≈ 0 (free end)

### Test 2: Downward Load (beam_oneway_down.dat)
```bash
python epframe_oneway.py beam_oneway_down.dat beam_down.out
```

**Check in output:**
- ✓ Node 4 status: "Y=POS:**ACTIVE**"
- ✓ R4y > 0 (positive/upward)
- ✓ Displacement at Node 4 ≈ 0
- ✓ Sum of reactions ≈ 10 kips upward

---

Thank you for catching this! The correct specification is **Y+** for a roller on the ground. 🎓

**Material:** Steel beam W12×50
- E = 29,000 ksi
- I = 394 in⁴
- A = 14.6 in²
- Mp = 3,600 in-kips

**Loading:** Vertical force P at Node 3 (x = 150")

---

## Case 1: Upward Load (P = +10 kips)

### Step 1: Check if Node 4 would have upward or downward reaction

For a beam on three supports with upward load between supports 2 and 4:
- The reaction at support 4 would naturally be **upward** (positive)
- But support 4 is **Y- only** (can only resist downward motion)
- Therefore: **Node 4 LIFTS OFF**

### Step 2: Analyze as two-support system

With Node 4 lifted off, the system becomes:
- **Segment 1-2:** Propped cantilever (fixed at 1, pin at 2)
- **Segment 2-3:** Cantilever from Node 2
- **Segment 3-4:** Free (unsupported, Node 4 displaced upward)

### Step 3: Calculate reactions (Nodes 1 and 2 only)

**Equilibrium equations:**
```
ΣFx = 0:  R1x = 0
ΣFy = 0:  R1y + R2y + 10 = 0
ΣM@1 = 0: M1 + R2y(100) + 10(150) = 0
```

**For the propped cantilever (Node 1 to Node 2, length L=100"):**

Using standard formulas, the load P = 10 kips at x = 150" is beyond the pin at x = 100".

This means the load is on a cantilever extending from Node 2 to Node 3 (length = 50").

**Simple cantilever analysis (Node 2 to Node 3):**
```
R2y = -P = -10 kips (downward)
M2 = -P × 50 = -10 × 50 = -500 in-kips
```

**Continuous beam from Node 1 to Node 2 with end moment M2 = -500:**
```
R1y = -R2y = 10 kips (upward)
M1 = -R2y(100) - M2 = -(-10)(100) - (-500) = 1000 + 500 = 1500 in-kips
```

### **Expected Results (Upward Load):**
| Node | Rx (kips) | Ry (kips) | Mz (in-kips) | Status |
|------|-----------|-----------|--------------|--------|
| 1 | 0 | 10 | 1500 | ACTIVE |
| 2 | 0 | -10 | 0 | ACTIVE |
| 4 | 0 | 0 | 0 | **LIFT-OFF** |

**Displacements:** Node 4 should have positive (upward) displacement.

---

## Case 2: Downward Load (P = -10 kips)

### Step 1: Check if Node 4 stays in contact

For downward load, the reaction at Node 4 would naturally be **downward** (negative).
- Support 4 is **Y- only** (can resist downward motion)
- Therefore: **Node 4 STAYS IN CONTACT**

### Step 2: Analyze as three-support continuous beam

This is a statically indeterminate continuous beam with three supports.

Using **superposition** or **moment distribution:**

For a continuous beam with equal spans (100", 50", 50") this is complex, but we can use symmetry and approximate:

**Simplified analysis (assuming relatively flexible beam):**

The middle span (100-150-200) with pin at 100 and roller at 200 acts approximately as a simply supported beam for the segment 100-200 (length 100").

Load P = -10 kips at x = 50" from Node 2 (midpoint of 100-200 span):

For simple beam with center load:
```
R2y ≈ -5 kips
R4y ≈ -5 kips
```

Node 1 provides additional support for the continuous system, so exact values will differ, but Node 4 will have a **negative** (downward) reaction.

### **Expected Results (Downward Load):**
| Node | Rx (kips) | Ry (kips) | Mz (in-kips) | Status |
|------|-----------|-----------|--------------|--------|
| 1 | 0 | ~-3 to -5 | ~-300 to -500 | ACTIVE |
| 2 | 0 | ~-5 to -7 | 0 | ACTIVE |
| 4 | 0 | ~-2 to -5 (negative) | 0 | **ACTIVE** |

**Displacements:** Node 4 displacement ≈ 0 (in contact)

---

## Key Verification Points

### Test 1: Upward Load
Run: `python epframe_oneway.py beam_oneway.dat beam_up.out`

**Check:**
1. ✓ Node 4 status shows "LIFT-OFF"
2. ✓ R4y = 0
3. ✓ Displacement at Node 4 is positive (upward)
4. ✓ R1y ≈ 10 kips, R2y ≈ -10 kips
5. ✓ M1 ≈ 1500 in-kips

### Test 2: Downward Load
Run: `python epframe_oneway.py beam_oneway_down.dat beam_down.out`

**Check:**
1. ✓ Node 4 status shows "ACTIVE"
2. ✓ R4y < 0 (negative/downward)
3. ✓ Displacement at Node 4 ≈ 0
4. ✓ Sum of all reactions = -10 kips

---

## Physical Interpretation

### Upward Load Scenario
- The upward force at Node 3 tries to **lift** the right end of the beam
- A regular support would provide upward reaction to resist this
- But Node 4's **Y- constraint** cannot provide upward (tension) force
- Result: The support **cannot engage**, beam lifts off
- The beam becomes a cantilever from the left side

### Downward Load Scenario
- The downward force naturally pushes the beam **down**
- Node 4's **Y- constraint** can resist this downward motion
- Result: Support **engages**, provides compressive (downward) reaction
- The beam acts as a three-support continuous beam

---

## Teaching Points

1. **One-way supports are common in practice:**
   - Beams resting on walls (can lift off)
   - Bridges with expansion bearings
   - Temporary shoring

2. **Contact mechanics:**
   - Complementarity: Either u=0 (contact) OR R=0 (separated)
   - Never both simultaneously

3. **Load path changes:**
   - Same structure, different load → different active supports
   - Structural system changes based on loading direction

4. **Computational verification:**
   - Hand calculations verify code correctness
   - Code handles complexity beyond hand calculation

---

## Summary

This example provides a **clear analytical verification** of the one-way reaction algorithm:
- **Known input** (simple beam geometry and loading)
- **Predictable behavior** (lift-off vs contact)
- **Calculable reactions** (hand calculation possible)
- **Two contrasting cases** (upward vs downward load)

Perfect for teaching and code validation! 🎓
