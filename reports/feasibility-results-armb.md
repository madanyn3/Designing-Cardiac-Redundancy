# M8 feasibility results — Arm B, fat-capable metabolism

2026-10-22. Protocol version 4. The structural alternative to
[Arm A](feasibility-results.md), which charges every joule of oxidative demand
to glucose. Arm B asks what changes when tissue may oxidise fatty acid.

Same conditions as Arm A unless stated: compression kappa 0.5, reference
physiology, oracle 0.5 s activation, 40 s settled horizon, and the nine
unresolved requirements attached to every verdict.

---

## 1. What Arm B changes, and what it does not

**The circulation.** Visceral adipose is **disaggregated** out of the residual
branch rather than added to the network, so total cardiac output, the baseline
fraction sum and total oxygen demand are all unchanged. The split is
hydraulically neutral by construction — two parallel paths at a shared driving
pressure have the same conductance as their parent — and the model confirms
0.98644 mL/s/mmHg for both arms. Arm B is therefore not a different
circulation; it is the same circulation with one territory named.

**The fuel.** Fatty-acid supply is limited by two declared terms, with the
binding one reported rather than assumed: lipolytic capacity, and adipose
perfusion multiplied by the arterio-venous concentration difference. **At every
flow examined, perfusion binds.** Fat availability is therefore a circulatory
question, not a metabolic one.

**The brain is excluded on physiological grounds.** Long-chain fatty acid does
not cross the blood-brain barrier in useful quantity. Ketone bodies would
change that, but ketogenesis takes days of fasting against a horizon of minutes
to hours, so ketones are a separate long-horizon axis and contribute nothing
here. The brain's demand stays wholly on glucose, which sets a floor on the
glucose fraction at its share of total demand.

**What Arm B does not change:** the pressure screen, the self-supply
requirement, the horizon, the startup allowance, or any design control. The two
arms differ by the network variant, the glucose share array and the adipose
policy, and by nothing else.

## 2. The joint optimum, against Arm A

Same search, same objective, same constraints, on the fat-capable model:
32 Sobol points plus 3 compass iterations, 50 evaluations, 23 feasible, 103 min.

| | Arm A | **Arm B** |
|---|---|---|
| Cavity scale | 0.566 | **0.443** |
| Rate | 49.2 /min | 67.8 /min |
| Coronary conductance | 1.027 | 1.591 |
| Auxiliary O₂ demand | 7.600 mL/min | **7.609** |
| MAP | — | 66.69 mmHg |
| Glucose fraction | 1.000 | **0.751** |
| **Supported duration** | 33.7 min | **47.1 min** |
| At a 3.0 mmol/L floor | 13.5 min | **18.9 min** |

The cavity and rate differences in this table are **not** attributable to the
fuel model; see below. Only the glucose fraction and the duration rows are.

**Fat buys 40 per cent more supported duration. It does not buy a smaller pump.**

An earlier draft of this section claimed the 22 per cent smaller cavity as a
fat effect. A direct cross-check refutes that. Arm B's optimum design — scale
0.443, 67.8 beats/min, conductance 1.591 — was evaluated under **Arm A's
glucose-only metabolism** and is feasible there: pressure satisfied,
self-supply satisfied, sustained unmet 0.0, at 7.627 mL O₂/min against Arm A's
reported 7.600. The design needs no fat at all.

So the two arms' optima differ by **what their searches sampled**, not by the
fuel. The objective surface is a flat valley, the two searches used different
seeds, and they landed in different parts of it. Arm A's joint search found its
smallest feasible design at scale 0.514 and sampled nine points near 0.443,
none at the combination of rate and conductance that works there.

What survives as a genuine fuel effect is the substrate budget: the glucose
fraction falls from 1.000 to 0.751, and supported duration rises from 33.7 to
47.1 minutes, 13.5 to 18.9 against a 3.0 mmol/L floor. Those follow directly
from fat displacing glucose demand and are not sampling artefacts.

**Arm B changes the fuel budget, not the machine.**

As in Arm A the optimum is a **shallow valley**: the five leading designs share
scale 0.443 and 67.8 beats/min while conductance spans 1.231 to 1.591 and the
objective varies 0.12 per cent.

## 3. Fat is nowhere near saturating

At the optimum the depot receives **52.5 mL/min**, delivering 27.0 mL O₂-eq/min
and covering only **12.6 per cent** of non-obligate demand. Lipolytic capacity
never binds at any flow examined; adipose perfusion does.

| Adipose policy | Depot flow | Fat supply | Non-obligate covered |
|---|---|---|---|
| Preserved (0.25) | 13.5 mL/min | 7.0 mL O₂-eq/min | 3.2% |
| Perfused (1.0) | 52.5 | 27.0 | 12.6% |

Perfusing the depot costs **0.91 mmHg** of mean pressure (67.70 to 66.79 at a
fixed design), consistent with the 2.1 per cent rise in effective systemic
conductance it produces. That is the price of the fuel, and it is modest.

The ceiling, if fat covered all non-obligate demand, would be a glucose
fraction of 0.266 — the brain's share — and roughly 216 min. **The gap between
12.6 per cent and that ceiling is entirely a perfusion limit**, so the result is
bounded by how much of the circulation the support policy routes to fat, not by
how much fat there is.

## 4. The adipose policy is a capability claim, not a design choice

The depot's vascular resistance is the subject's anatomy and is identical in
both policies; what the support policy elects is the fraction of that
conductance the circulation opens. The choice is categorical — adipose joins
the mandatory set or stays with the preserved territories — because any value
above full conductance would mean dilating a vessel beyond its physiological
maximum.

The policy is therefore **not searched**. The joint search runs once per policy
and the two are compared, rather than letting an optimizer select the
assumption that flatters the result.

## 5. What is not established

- Everything in Arm A's limitations section still applies: latency, kappa,
  scenario and horizon dependence are unchanged by the fuel model.
- **The arterio-venous fatty-acid difference is assumed, not sourced.** It is
  set at 1.0 mmol/L, a mid value against reported plasma concentrations of
  0.5-1.0 fasting and 1.5-2 in prolonged fasting. Fat supply scales linearly
  with it, so it deserves the evidence treatment the pressure criterion and the
  glucose floor received.
- **The visceral depot mass is assumed** at 2.0 kg, and its perfusion rate is a
  declared axis over the reported 3-5 mL/100 g/min fasting range.
- Fatty acid is a **supply controller, not a conserved species**. No
  fatty-acid concentration is tracked per compartment, in the same way the
  hepatic glucose export law is bookkeeping rather than hormone kinetics.
- Ketone contribution is zero here, which is correct for this horizon and
  wrong for a multi-day one.
- Red cells and renal medulla are obligate glucose consumers physiologically
  but are not represented as separate organs, so the brain is the only tissue
  the model holds to that constraint.
