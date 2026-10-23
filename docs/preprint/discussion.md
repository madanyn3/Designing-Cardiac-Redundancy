# Discussion

Draft for the preprint. The study asks whether an auxiliary biological pump
could maintain vital-organ perfusion after loss of native cardiac pumping, and
answers it conditionally: **within a declared envelope, for tens of minutes,
limited by fuel logistics rather than by pumping capacity, and requiring a
detection speed no demonstrated system achieves.** Each clause carries a
condition, and the conditions are the result.

## 1. What the answer depends on

Three quantities decide the outcome, and the study's contribution is largely in
identifying which.

**Cavity size, which was the unexamined variable.** An earlier negative result
from the same model — no configuration of 84 reaching the pressure screen —
described a design space in which cavity scale was held at a quarter-scale test
fixture. The parameter was always legal and the inventory never constrained it.
A negative result about a design space is only as wide as the space actually
varied, and this one was narrower than it appeared.

**Coronary self-compression, which rests on no measurement.** The fraction of
auxiliary chamber pressure applied to its own coronary bed as closing pressure
is the single most influential quantity in the self-supply constraint. At the
maximum, self-supply fails in every design tested; at the minimum it is benign.
Feasibility holds below roughly the midpoint and fails above it. This is a
one-sided wall rather than a narrow band, which matters: the question is not
whether the true value sits in a window but whether it sits below a threshold.
No measurement addresses it, and obtaining one is the highest-ranked evidence
priority carried out of this work.

**Hepatic glucose export rate, which sets duration and is not a design
variable.** Duration is limited not by stored glycogen — which would outlast
the circulating pool by two orders of magnitude — but by the rate at which
glucose can be released into circulation. The circulation delivers five times
more hepatic flow than export can use. Supported duration consequently moves
by only a few minutes across the entire admissible cavity range, while moving
by a factor of three on the declared glucose floor, a reporting choice the
model does not set.

## 2. The latency result, against existing support

The design requirement the study hands back is a detection-plus-recruitment
budget of **0.8 to 1.6 seconds**, depending on which design is built, against
an implemented controller in this same model that releases at 2.55 s. Nothing
here demonstrates that sub-second detection is achievable.

This requirement looks severe until it is compared with the alternative.
Extracorporeal resuscitation, the nearest existing approach to restoring
circulation after refractory arrest, reports times from call or arrest to
established flow of **59 minutes on average** in one randomised trial and a
**74-minute median** in another. Those intervals are not no-flow times and
their start definitions differ, so they are not directly comparable with a
detection latency — but they are three orders of magnitude away from the
requirement here, and they bound a different problem.

The comparison clarifies what an implanted dormant pump would actually trade:
**a logistics chain measured in tens of minutes for a detection problem
measured in seconds.** That is a favourable exchange if and only if the
detection problem is solvable, which this study cannot determine. It does
establish that capacity buys latency tolerance, and that the design which
avoids the study's least-supported assumption also tolerates roughly twice the
delay — so the two robustness properties are purchased together rather than
traded against one another.

## 3. Why the fuel result is not a statement about fasting

Supported duration of 33.8 minutes on glucose alone, or 47.1 with fat
available, sits oddly beside the observation that humans tolerate days of
fasting. The two are not in conflict, and the reason is specific. In fasting,
demand is met by a slow, hormonally regulated shift toward fat oxidation at
rest. Here, demand is elevated, the horizon is minutes, and the model charges
oxidative demand to substrate that must arrive by concurrent delivery, with no
tissue buffer. What binds is a *rate* of mobilisation, not a *quantity* of
reserve.

This is also where the structural alternative earns its place. Permitting
fatty-acid oxidation raises duration by 40 per cent and lowers the glucose
fraction to 0.751 — but it does **not** change the machine: under the declared
pressure margin the two fuel arms cost the same to within 0.3 per cent. What
changes is the binding resource, from hepatic export to **adipose perfusion**.
Lipolytic capacity never binds at any flow examined. Fat covers 12.6 per cent
of non-obligate demand against a ceiling of roughly 216 minutes if it covered
all of it, and the entire gap is circulatory.

**Fat availability under arrest is therefore a perfusion question, not a
metabolic one** — a result that would be invisible to an analysis treating fuel
as a budget rather than as something the circulation has to deliver.

## 4. A methodological result worth separating from the physiology

Four of this study's own conclusions were corrected during the work, three of
them from a single cause, and the cause generalises past this model.

When the objective is flat and a constraint is a hard threshold, minimizing
cost returns **the most fragile design that technically qualifies**. The
optimum located here satisfied every screen while holding 0.005 mmHg of
margin; a 0.04 per cent change in its activation period moved it to violated.
Disconfirming that design then measures the specification rather than the
physiology: a design sitting on a threshold moves off it under any
perturbation, and designs within 0.2 per cent of one another reorder under any
change. Two conclusions drawn that way — that no design survives minimum
compression, and that the selection criterion was unsound — both reversed when
the stressed sample was drawn from across the feasible set instead of from its
cheapest corner.

The fix is cheap and should be routine: declare a minimum constraint margin
before searching, report its cost, and require robustness tests to span the
feasible set. Here the margin cost 7.1 per cent of the objective. Reporting a
design that meets a threshold by an arbitrarily small amount is reporting the
threshold, not the design.

A separate instance is worth stating for the same reason. Caching keyed on
less than the full question silently answers a different one: a key missing
the screen horizon paired 40-second pressure verdicts with 10-second
self-supply verdicts, and a key missing the design returned a superseded
design's verdicts byte-identical. Both were found by cross-checks rather than
by the runs that produced them.

## 5. What this does not establish

The configuration studied does not exist, so nothing here is validated against
measurement from an arrested subject, and no result is a statement about
survival, neurological outcome, injury avoidance or clinical applicability.
Specifically:

- **No claim of constructibility.** Cavity scale is a capacity proxy, not a
  determined anatomical mass, and no developmental, immunological or surgical
  pathway is modelled.
- **No injury account.** Ischaemic and reperfusion injury, thrombosis, blood
  compatibility, CO₂ clearance and pH are not represented. A cumulative unmet
  demand is an exposure proxy and never an injury dose.
- **No regional adequacy.** Organ-level pressure minima are unresolved, so a
  satisfied central-aortic screen is not evidence of adequate perfusion in any
  particular territory, including the brain.
- **Mechanisms declared absent rather than assumed benign:** supply-limited
  contractility, finite tissue glucose transport and refill, dynamic
  autoregulation and intracranial coupling, tissue injury and reperfusion,
  thrombosis, and dormant biological viability. Contraction is prescribed, so
  unmet demand never slows the pump — a conservative choice for pressure and an
  optimistic one for metabolism.
- **Two named reproduction limitations**, R-M3-01 and L-M4-F01, stated in the
  Verification section.
- **Fuel-arm assumptions**: the arterio-venous fatty-acid difference and the
  visceral depot mass are assumed rather than sourced, and fat supply scales
  linearly with the former. Ketone contribution is zero, correct for minutes
  and wrong for days.

Every boundary reported is a boundary of what was tested. Where the study says
a design fails, it means at the values examined, and where it says one
succeeds, it means best found within the declared bounds and budget.

## 6. What would change the answer

In priority order, and with the decision each would change:

1. **Joint auxiliary pressure, flow, work and coronary supply measurements**,
   which would settle the compression fraction and so decide whether the
   feasible region exists at all.
2. **Regional pressure–flow and closing-pressure relationships**, which would
   replace the central-aortic screen with organ-level adequacy and could
   invalidate or confirm the single most load-bearing simplification.
3. **Tissue extraction and mobilisation rates, especially in the auxiliary
   myocardium**, which set the self-supply constraint that bounds the envelope
   from above.
4. **Sustained hepatic output and accessible glycogen**, which set the duration
   endpoint directly.
5. **A demonstrated detection latency** for loss of native pumping, which is
   the one requirement this study converts from assumption into specification
   and cannot itself evaluate.

The study's useful output is not the located design. It is the identification
of which three quantities decide the answer, the demonstration that two of them
are physiological rather than engineering choices, and a bounded statement of
how long the arrangement could work if they fell favourably.
