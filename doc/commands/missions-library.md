# Missions — To Mars & Back

Pedagogical Earth⇄Mars round-trip planners, built on the Astronautics transfer
machinery (universal-variable Lambert, gravity assist) and the heliocentric
ephemerides. Dates are Julian Days; ΔV in km/s. Two mission classes are covered:
**conjunction** (long stay, lowest ΔV) with RTPlan, and **opposition** (short stay,
Venus gravity-assist return) with OppRoute and RTPlanOpp.

## RTPlan

Conjunction-class Earth⇄Mars round-trip optimiser (heliocentric ΔV). It optimises the
outbound Earth→Mars leg, then the return Mars→Earth leg departing after a minimum Mars
surface stay, and returns the full itinerary. Slow, so the example sets a fast search
precision first. It reproduces the real Perseverance outbound window.

`t1lo t1hi stay` → `{ t1 t2 stay t3 t4 ΔVout ΔVret ΔVtot }`

```rpl
1 'AstronTXPrecision' STO
2459030 2459075 500 RTPlan
@ Expecting { 2459055.6 2459259.8 517.6 2459777.4 2460049.4 6.318 6.505 12.823 }
@ Earth 2020-07-25 → Mars 2021-02-14, stay 517 d, home 2023-04-14; ΔVtot 12.82 km/s
```

## RTPlanOpp

Best Venus-flyby Mars→Earth return over a Mars-departure window (opposition class). A
two-stage search (coarse grid then local refine) that minimises a combined cost
favouring a ballistic and turn-feasible Venus swing-by. Keep the window narrow (~60 d);
it is slow (~2-4 min). The result is feasible when matchErr is near zero and turnReq
does not exceed turnMax.

`t3lo t3hi` → `{ t3 tv t4 ΔVTEI matchErr turnReq turnMax tof }`

```rpl
2459920 2459980 RTPlanOpp
@ Expecting { 2459930 2460150 2460220 5.90 0.00084 30.93 49.72 290 }
@ Mars 2022-12-16 → Venus 2023-07-24 → Earth 2023-10-02; Earth entry only ~4.9 km/s
```

## OppRoute

Evaluate one Mars→Venus→Earth gravity-assist return: two Lambert arcs patched by a Venus
swing-by. It returns the Mars-departure ΔV, the v∞ mismatch at Venus (near zero for a
valid ballistic flyby), the required versus achievable turn angle, and the return time.
Instant, no search. A route is feasible when matchErr is near zero and turnReq does not
exceed turnMax.

`t3 tv t4` → `{ ΔVTEI matchErr turnReq turnMax tof }`

```rpl
2459950 2460160 2460240 OppRoute
@ Expecting { 6.857 0.268 29.61 37.70 290 } — feasible since turnReq 29.6 ≤ turnMax 37.7
```

## ♀Pf

Venus heliocentric position vector (au) from a Julian Day — the Venus position function,
sibling of the Earth and Mars position functions used by the transfer tools.

`JD` → `[x y z]`

```rpl
2459950 ♀Pf     @ Expecting [ 0.6140036 -0.3890710 -0.0407721 ]  (au, 2023-01-05)
```

## MarsRoundTrip

Compare a conjunction versus an opposition Earth⇄Mars round trip. It runs RTPlan
(conjunction) and RTPlanOpp (opposition Venus-flyby return), computes the opposition
Earth-entry speed, and leaves five tagged values on the stack for a side-by-side read.
Slow (~90 s) — it chains both optimisers. Then explore your own windows with RTPlan and
RTPlanOpp directly.

(no input) → `Conj_years Conj_dVtot Opp_return_d Opp_TEI Opp_entry`

```rpl
MarsRoundTrip
@ Expecting Conj_years:2.72 Conj_dVtot:12.82 Opp_return_d:290 Opp_TEI:5.90 Opp_entry:5.44
@ Conjunction 2.72 yr / 12.82 km/s; opposition Venus return, Earth entry only 5.4 km/s.
```
