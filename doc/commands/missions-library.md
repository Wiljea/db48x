# Missions — Earth To Mars & Back

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
valid ballistic flyby), the required versus achievable turn angle, the return time, and the
Earth entry speed. Instant, no search. A route is feasible when matchErr is near zero and
turnReq does not exceed turnMax; a low EarthEntry marks a strong Venus assist (soft entry).

`t3 tv t4` → `{ ΔVTEI matchErr turnReq turnMax tof EarthEntry }`

```rpl
2459950 2460160 2460240 OppRoute
@ Expecting { 6.857 0.268 29.61 37.70 290 4.87 } — feasible (turnReq 29.6 ≤ turnMax 37.7);
@ EarthEntry only 4.87 km/s = a soft, Venus-assisted entry
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

## MissionΔV

End-to-end propulsive ΔV budget of a conjunction round trip. Feed it an RTPlan itinerary;
it recomputes the two transfer legs (TrCost) for their departure/arrival v∞, converts the
space burns from parking orbits with TrToOrbi, and adds representative atmospheric phases.
It returns eight tagged ΔV values (seven phases plus the total), in km/s. Chains after
RTPlan; the example uses a fixed conjunction itinerary so it runs instantly.

The seven propulsive phases, in chronological order:

* **AscentEarth** — Earth surface to Low Earth Orbit (the launch to a parking orbit).
* **TMI** — Trans-Mars Injection: from LEO onto the interplanetary transfer to Mars.
* **MOI** — Mars Orbit Insertion: braking capture from the arrival hyperbola into Mars orbit.
* **EDL_Mars** — Entry, Descent and Landing: from Mars orbit down to the surface.
* **AscentMars** — Mars surface back up to Mars orbit.
* **TEI** — Trans-Earth Injection: from Mars orbit onto the return transfer to Earth.
* **EntryEarth** — Earth atmospheric entry on return (aerobraked, so ~0).

The space phases (TMI, MOI, TEI) are computed rigorously; the atmospheric phases
(ascents, EDL, entry) are representative values (Earth ascent 9.4, Mars EDL 0.6, Mars
ascent 4.1, Earth entry 0 km/s), with parking orbits LEO 6678 km and Mars 3689 km.

`{ itinerary }` → `1_AscentEarth … 8_TOTAL`

```rpl
{ 2459055.617 2459259.801 517.558 2459777.359 2460049.376 6.31786 6.50523 12.82309 } MissionΔV
@ Expecting 8_TOTAL:22.47 km/s — AscentEarth 9.4, TMI 3.81, MOI 2.07, EDL 0.6,
@ AscentMars 4.1, TEI 2.48, Entry 0. A typical chemical Mars round-trip budget.
```

## ETMB1

Earth-To-Mars-&-Back mission 2022 (launch), Venus-flyby return — a worked example. It runs
the validated opposition-class return (Mars→Venus→Earth) at fixed dates via OppRoute and
returns the return-leg itinerary as tagged values. The original 2022 mission. Venus turns ~31°, giving a gentle Earth entry of ~5.4 km/s.

![ETMB1 trajectory](img/etmb1.bmp)

`(no input)` → `MarsDeparture VenusFlyby EarthArrival ReturnDays dVTEI turnReq EarthEntry`

```rpl
ETMB1
@ Expecting MarsDeparture:2459930 VenusFlyby:2460150 EarthArrival:2460220 ReturnDays:290
@ dVTEI_kms:5.90 turnReq_deg:30.9° EarthEntry_kms:5.44
```

## ETMB2

Earth-To-Mars-&-Back mission 2035 (launch), Venus-flyby return — a worked example. It runs
the validated opposition-class return (Mars→Venus→Earth) at fixed dates via OppRoute and
returns the return-leg itinerary as tagged values. **★ Featured mission.** The strongest twin found: Venus works even harder than in ETMB1 (~39° turn) for the *softest* entry — ~4.0 km/s, gentler than ETMB1 itself. It uniquely combines the Venus synodic resonance with a **perihelic Mars opposition** (Mars near the Sun), which is what makes the swing-by so strong — the best future crewed-return candidate in the ephemeris range.

![ETMB2 trajectory](img/etmb2.bmp)

`(no input)` → `MarsDeparture VenusFlyby EarthArrival ReturnDays dVTEI turnReq EarthEntry`

```rpl
ETMB2
@ Expecting MarsDeparture:2464635 VenusFlyby:2464825 EarthArrival:2464915 ReturnDays:280
@ dVTEI_kms:7.14 turnReq_deg:39.3° EarthEntry_kms:3.97
```

## ETMB3

Earth-To-Mars-&-Back mission 2042 (launch), Venus-flyby return — a worked example. It runs
the validated opposition-class return (Mars→Venus→Earth) at fixed dates via OppRoute and
returns the return-leg itinerary as tagged values. Venus-light: the 6.4-year resonance yields a feasible passage but a weak assist (~7° turn), leaving a brutal ~14.7 km/s entry.

![ETMB3 trajectory](img/etmb3.bmp)

`(no input)` → `MarsDeparture VenusFlyby EarthArrival ReturnDays dVTEI turnReq EarthEntry`

```rpl
ETMB3
@ Expecting MarsDeparture:2466965 VenusFlyby:2467105 EarthArrival:2467170 ReturnDays:205
@ dVTEI_kms:9.19 turnReq_deg:7.5° EarthEntry_kms:14.69
```

## ETMB4

Earth-To-Mars-&-Back mission 2048 (launch), Venus-flyby return — a worked example. It runs
the validated opposition-class return (Mars→Venus→Earth) at fixed dates via OppRoute and
returns the return-leg itinerary as tagged values. The cheapest full mission (outbound v∞ only ~3.6 km/s), with a moderate ~9.6 km/s entry.

![ETMB4 trajectory](img/etmb4.bmp)

`(no input)` → `MarsDeparture VenusFlyby EarthArrival ReturnDays dVTEI turnReq EarthEntry`

```rpl
ETMB4
@ Expecting MarsDeparture:2469315 VenusFlyby:2469470 EarthArrival:2469535 ReturnDays:220
@ dVTEI_kms:6.24 turnReq_deg:9.7° EarthEntry_kms:9.58
```
