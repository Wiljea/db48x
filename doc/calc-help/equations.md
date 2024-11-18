# Equations library

The DB48X calculator features a library of equations covering mathematics,
physics, chemistry and computer science. The built-in equations can be extended
using the `config/equations.csv` configuration file.


## Columns and Beams

The 21 variables in the Columns and Beams section are:

* `o`: Eccentricity (offset) of load (dim.: length)
* `σcr`: Critical stress (dim.: pressure=force/area, in SI: pascal, Pa)
* `σmax`: Maximum stress (dim.: pressure=force/area, in SI: pascal, Pa)
* `θ`: Slope at `x` (dim.: angle)
* `A`: Cross-sectional area
* `a`: Distance to point load
* `ε`: Eccentricity (dim.: length)
* `c`: Distance to edge fiber ([Eccentric Columns](#Eccentric Columns)), or Distance to applied moment ([Beams](#Beams))
* `E`: Modulus of elasticity (dim.: pressure=force/area, in SI: pascal, Pa)
* `I`: Moment of inertia (dim.: length^4, in SI: m^4)
* `K`: Effective length factor of column
* `L`: Length of column or beam
* `M`: Applied moment (dim.: length·force, in SI: N·m)
* `Mx`: Internal bending moment at x (dim.: length·force, in SI: N·m)
* `P`: Load (Eccentric Columns), or Point load (beams) (dim.: force)
* `Pcr`: Critical load (dim.: force)
* `r`: Radius of gyration
* `V`: Shear force at x
* `w`: Distributed load (dim.: force)
* `x`: Distance along beam
* `y`: Deflection at x (dim.: length)

For simply supported beams and cantilever beams ([Simple Deflection](#Simple Deflection) through [Cantilever Shear](#Cantilever Shear)), the calculations differ depending upon the location of `x` relative to the loads.

* Applied loads are positive downward.
* The applied moment is positive counterclockwise.
* Deflection is positive upward.
* Slope is positive counterclockwise
* Internal bending moment is positive counterclockwise on the left-hand part.
* Shear force is positive downward on the left-hand part.

### Elastic Buckling

These equations apply to a slender column (`K·L/r>100`) with length factor `K`.

![Elastic Buckling](img/ElasticBuckling.bmp)

* To calculate `[Pcr_kN;σcr_kPa]` (Critical load; Critical stress) from 6 known variables:
```rpl
L=7.3152_m  r=4.1148_cm  E=199947961.502_kPa  A=53.0967_cm^2  K=0.7  I=8990598.7930_mm^4
@ Expecting [ Pcr=676.60192 6324 kN σcr=127 428.24437 8 kPa ]
'ROOT(ⒺElastic Buckling;[Pcr;σcr];[1_kN;1_kPa])'
```

### Eccentric Columns

These equations apply to a slender column (`K·L/r>100`) with length factor `K`.

![Eccentric Columns](img/EccentricColumns.bmp)

* To calculate `[σmax_kPa;I_mm^4]` (Maximum stress; Moment of inertia) from 8 known variables:
```rpl
L=6.6542_m  A=187.9351_cm^2  r=8.4836_cm  E=206842718.795_kPa  K=1  P=1908.2571_kN  c=15.24_cm  ε=1.1806_cm
@ Expecting [ σmax=140 853.09700 6 kPa I=135 259 652.161 mm↑4 ]
'ROOT(ⒺEccentric Columns;[σmax;I];[1_kPa;1_mm^4])'
```

### Simple Deflection

![Simple Deflection](img/SimpleDeflection.bmp)

* To calculate `[y_in]` (Deflection at x) from 9 known variables:
```rpl
L=20_ft  E=29000000_psi  I=40_in^4  a=10_ft  P=674.427_lbf  c=17_ft  M=3687.81_ft*lbf  w=102.783_lbf/ft  x=9_ft
@ Expecting [ y=-0.60048 54094 96 in ]
'ROOT(ⒺSimple Deflection;[y];[1_in])'
```

### Simple Slope

![Simple Slope](img/SimpleSlope.bmp)

* To calculate `[Θ_°]` (Slope at `x`) from 9 known variables:
```rpl
L=20_ft  E=29000000_psi  I=40_in^4  a=10_ft  P=674.427_lbf  c=17_ft  M=3687.81_ft*lbf
 w=102.783_lbf/ft  x=9_ft
@ Expecting [ θ=-0.08763 17825 27 ° ]
'ROOT(ⒺSimple Slope;[θ];[0_°])'
```

### Simple Moment

![Simple Moment](img/SimpleMoment.bmp)

* To calculate `[Mx_ft*lbf]` (Internal bending moment at `x`) from 7 known variables:
```rpl
L=20_ft  a=10_ft  P=674.427_lbf  c=17_ft  M=3687.81_ft*lbf  w=102.783_lbf/ft  x=9_ft
@ Expecting [ Mx=9 782.1945 lbf·ft ]
'ROOT(ⒺSimple Moment;[Mx];[1_ft*lbf])'
```

### Simple Shear

![Simple Shear](img/SimpleShear.bmp)

* To calculate `[V_lbf]` (Shear force at `x`) from 6 known variables:
```rpl
L=20_ft  a=10_ft  P=674.427_lbf  M=3687.81_ft*lbf  w=102.783_lbf/ft  x=9_ft
@ Expecting [ V=624.387 lbf ]
'ROOT(ⒺSimple Shear;[V];[1_lbf])'
```

### Cantilever Deflection

![Cantilever Deflection](img/CantileverDeflection.bmp)

* To calculate `[y_in]` (Deflection at `x`) from 9 known variables:
```rpl
L=10_ft  E=29000000_psi  I=15_in^4  P=500_lbf  M=800_ft*lbf  a=3_ft  c=6_ft  w=100_lbf/ft  x=8_ft
@ Expecting [ y=-0.33163 03448 28 in ]
'ROOT(ⒺCantilever Deflection;[y];[0_in])'
```

### Cantilever Slope

![Cantilever Slope](img/CantileverSlope.bmp)

* To calculate `[Θ_°]` (Slope at `x`) from 9 known variables:
```rpl
L=10_ft  E=29000000_psi  I=15_in^4  P=500_lbf  M=800_ft*lbf  a=3_ft  c=6_ft  w=100_lbf/ft  x=8_ft
@ Expecting [ θ=-0.26522 01876 49 ° ]
'ROOT(ⒺCantilever Slope;[θ];[0_°])'
```

### Cantilever Moment

![Cantilever Moment](img/CantileverMoment.bmp)

* To calculate `[Mx_ft*lbf]` (Internal bending moment at `x`) from 7 known variables:
```rpl
L=10_ft  P=500_lbf  M=800_ft*lbf  a=3_ft  c=6_ft  w=100_lbf/ft  x=8_ft
@ Expecting [ Mx=-200. lbf·ft ]
'ROOT(ⒺCantilever Moment;[Mx];[1_ft*lbf])'
```

### Cantilever Shear

![Cantilever Shear](img/CantileverShear.bmp)

* To calculate `[V_lbf]` (Shear force at `x`) from 5 known variables:
```rpl
L=10_ft  P=500_lbf  a=3_ft  x=8_ft  w=100_lbf/ft
@ Expecting [ V=200. lbf ]
'ROOT(ⒺCantilever Shear;[V];[1_lbf])'
```


## Electricity

The 77 variables in the Electricity section are:

* `∈r`: Relative permittivity
* `μr`: Relative permeability
* `ω`: Angular frequency (dim.: angle/time)
* `ω0`: Resonant angular frequency (dim.: angle/time)
* `φ`: Phase angle
* `φp`: Parallel phase angle
* `φs`: Series phase angle
* `θ1`: First subtended angle relative to the left end of the wire
* `θ2`: Second subtended angle relative to the right end of the wire
* `ρ`: Resistivity (dim.: resistance·length; in SI: ohm·meter, Ω·m), or Volumic charge density [Drift Speed & Current Density](#Drift Speed & Current Density) (dim.: charge/volume, in SI: C/m^3)
* `ρ0`: Resistivity at the reference temperature T0 (dim.: resistance·length; in SI: ohm·meter, Ω·m)
* `ΔI`: Current Change (dim.: charge/time; in SI: ampere, A)
* `Δt`: Time Change
* `ΔV`: Voltage Change (dim.: energy/charge; in SI: volt, V)
* `A`: Wire cross-section area ([Wire Resistance](#Wire Resistance)), or Solenoid cross-section area ([Solenoid Inductance](#Solenoid Inductance)), or Plate area ([Plate Capacitor](#Plate Capacitor))
* `C`: Capacitance (dim.: charge^2/energy; in SI: farad, F)
* `C1`: First capacitance (dim.: charge^2/energy; in SI: farad, F)
* `C2`: Second capacitance (dim.: charge^2/energy; in SI: farad, F)
* `Cp`: Parallel capacitance (dim.: charge^2/energy; in SI: farad, F)
* `Cs`: Series capacitance (dim.: charge^2/energy; in SI: farad, F)
* `d`: Plate separation (dim.: length)
* `E`: Energy
* `Er`: Electric field at distance r from the source (dim.: force/charge; in SI: N/C = V/m)
* `Ep`: Electric field over an infinite plate (dim.: force/charge; in SI: N/C = V/m)
* `Ein`: Electric field inside a plate capacitor (dim.: force/charge; in SI: N/C = V/m)
* `F`: Force between charges
* `f`: Frequency (dim.: time^-1)
* `f0`: Resonant frequency (dim.: time^-1)
* `I`: Current, or Total current ([Current Divider](#Current Divider)) (dim.: charge/time; in SI: ampere, A)
* `I1`: Current in R1 (dim.: charge/time; in SI: ampere, A)
* `Imax`: Maximum current (dim.: charge/time; in SI: ampere, A)
* `L`: Inductance (dim.: energy/current^2; in SI: henry, H), or Length ([Wire Resistance](#Wire Resistance), [Cylindrical Capacitor](#Cylindrical Capacitor))
* `Lp`: Parallel inductance (dim.: energy/current^2; in SI: henry, H)
* `Ls`: Series inductance (dim.: energy/current^2; in SI: henry, H)
* `N`: Number of turns
* `P`: Power (dim.: energy/time)
* `q`: Charge (in SI: coulomb, C)
* `q1`: First point charge (in SI: coulomb, C)
* `q2`: Second point charge (in SI: coulomb, C)
* `qtest`: Test point charge (in SI: coulomb, C)
* `Qp`: Parallel quality factor
* `Qs`: Series quality factor
* `r`: Charge distance
* `R`: Resistance (dim.: energy·time/charge^2; in SI: ohm, Ω)
* `R1`: First resistance (dim.: energy·time/charge^2; in SI: ohm, Ω)
* `R2`: Second resistance (dim.: energy·time/charge^2; in SI: ohm, Ω)
* `ri`: Inside radius
* `ro`: Outside radius
* `Rp`: Parallel resistance (dim.: energy·time/charge^2; in SI: ohm, Ω)
* `Rs`: Series resistance (dim.: energy·time/charge^2; in SI: ohm, Ω)
* `t`: Time
* `ti`: Initial time
* `tf`: Final time
* `V`: Voltage, or Total voltage (dim.: energy/charge; in SI: volt, V)
* `V1`: Voltage accross R1 (dim.: energy/charge; in SI: volt, V)
* `Vi`: Initial voltage (dim.: energy/charge; in SI: volt, V)
* `Vf`: Final voltage (dim.: energy/charge; in SI: volt, V)
* `Vmax`: Maximum voltage (dim.: energy/charge; in SI: volt, V)
* `XC`: Reactance of capacitor (dim.: energy·time/charge^2; in SI: ohm, Ω)
* `XL`: Reactance of inductor (dim.: energy·time/charge^2; in SI: ohm, Ω)
* `σ`: Conductivity (dim.: (resistance·length)^-1; in SI: siemens/meter, S/m), or Surface charge density ([E field infinite plate](#E field infinite plate), [Plate capacitor](#Plate capacitor)) (dim.: charge/area, in SI: C/m^2)
* `T`: Temperature (in SI: K)
* `T0`: Fixed reference temperature for the measurement of resistivity (In SI: K)
* `uE`: Volumic density of electric energy (dim.: energy/volume)
* `αT`: Temperature coefficient of resistivity (dim.: temperature^-1, in SI: K^-1)
* `λ`: Linear charge density (dim.: charge/length, in SI: C/m)
* `J`: Current density (dim.: current/area, in SI: A/m^2)
* `Je`: Electron current density (dim.: current/area, in SI: A/m^2)
* `Jh`: Hole current density (dim.: current/area, in SI: A/m^2)
* `vd`: Drift speed
* `n`: Charge-carrier number density (dim.: volume^-1)
* `ne`: Electron number density (dim.: volume^-1)
* `nh`: Hole number density (dim.: volume^-1)
* `μe`: Electron mobility (dim.: speed/(electric field), in SI: m^2/(V·s))
* `μh`: Hole mobility (dim.: speed/(electric field), in SI: m^2/(V·s))
* `τc`: Mean free time between collisions
* `meeff`: Electron effective mass
* `mheff`: Hole effective mass

### Coulomb’s Law & E Field

These equations describe the electrostatic force between two point charged particles and the electric field observed at the position of a test charge which replaces one of the two charges `q1` or `q2` in the expression of the electric force. A finite object carrying a net charge `q1` can be considered as a point charge if the distance to the position of the point charge `q2` is much greater than the object dimension, see example 2, for the approximate calculations of the electric force and electric field far away from a charged plate.

* **Example 1**. To calculate `[F_N;Er_N/C]` (Electric force; Electric Field) from 5 known variables:
```rpl
q1=1.6E-19_C  q2=1.6E-19_C  r=4.00E-13_cm  εr=1  qtest=1.6E-19_C
@ Expecting [ F=14.38008 28598 N Er=8.98755 17873 8⁳¹⁹ N/C ]
'ROOT(ⒺCoulomb’s Law & E Field;[F;Er];[1_N;1_N/C])'
```

* **Example 2**. A square metal plate `L = 8_cm` on a side carries a charge of `q1 = 6_μC`. Approximate values of the electric force & electric field for a point charge `q2 = 1_μC` located at `r = 3_m` can be calculated with Coulomb’s law if the separation distance is much greater than the plate dimension `r >> L`. The whole plate is indeed considered as being a point charge providing that `r > 10 · L`. Therefore, to calculate `[F_N;Er_N/C]`:
```rpl
L=8_cm r=3_m q1=6E-6_C  q2=1E-6_C  r=3_m  εr=1  qtest=1E-6_C
@ Expecting [ F=5.99170 11915 9⁳⁻³ N Er=5 991.70119 159 N/C ]
if 'r > 10*L' then
 'ROOT(ⒺCoulomb’s Law & E Field;[F;Er];[1_N;1_N/C])'
end
```

### E Field Infinite Line

The expression for the radial electric field at the distance `r` is approximately valid if this distance is such that `r << L` and therefore also applies to a wire of finite length `L`.

* To calculate `[λ_C/m;Er_N/C]` (Linear charge density; Electric Field) from 4 known variables:
```rpl
Q=5E-6_C  L=3_m  r=0.05_m  εr=1
@ Expecting [ λ=1.66666 66666 7⁳⁻⁶ C/m Er=599 170.11915 9 N/C ]
'ROOT(ⒺE Field Infinite Line;[λ;Er];[1_C/m;1_N/C])'
```
The code below saves the reference value for comparison with the example 2 in [E Field Finite Line](#E Field Finite Line):
```rpl
@ Save the reference value for comparison below
Er0=Er
```

### E Field Finite Line

The expression of the radial electric field at the distance `r` depends on the subtended angles `θ1` and `θ2` relative to the ends of the wire of finite length `L`.

![E field finite line](img/EFieldFiniteLine.bmp)

* **Example 1.** To calculate `[λ_C/m;Er_N/C]` (Electric Field; Linear charge density) from 6 known variables and also with the distance `r=(L/2)/tanθ1` and angle `θ2=360°-θ1` (see figure):
```rpl
r='(3_m)/(2*tan 30_°)' θ2='360_°-30_°'
Q=5E-6_C  L=3_m  r=2.5981_m  εr=1  θ1=30_°
@ Expecting [ λ=1.66666 66666 7⁳⁻⁶ C/m Er=5 765.46436 972 N/C ]
'ROOT(ⒺE Field Finite Line;[λ;Er];[1_C/m;1_N/C])'
```

* **Example 2.** To show that the infinite line of the previous section can approximate the finite case (if `r << L` realised when `r < L/10`), we calculate `[λ_C/m;Er_N/C]` with the angles `θ1=ATAN((L/2)/r)` and `θ2=360°-θ1` (see figure):
```rpl
Q=5E-6_C  L=3_m  r=5_cm  εr=1  θ1='atan(L/2/r)' θ2='360_°-θ1'
if 'r < L/10' then
@ Expecting [ λ=1.66666 66666 7⁳⁻⁶ C/m Er=598 837.52400 7 N/C ]
'ROOT(ⒺE Field Finite Line;[λ;Er];[1_C/m;1_N/C])'
end
```
Verify relative difference under condition `5_cm << 3_m` with the example of ([E Field Infinite Line](#E Field Infinite Line))
```rpl
Er0 Er %Ch
@ Expecting [ -5.55093 02084 6⁳⁻² ]
@ % of relative difference 
```

### E Field Infinite Plate

The expression of the perpendicular electric field is constant over an infinite plate and can approximate the field at a distance `d` from a finite plate if it is very small compare to the dimensions (length or width `L`) of the plate. On the contrary if `d >> L`, `Ep` can be approximated if we consider the whole plate as being a point charge with `q = σ·A` (where `σ` is the surface charge density), see example 2 of [Coulomb’s Law & E Field](#Coulomb’s Law & E Field).

* To calculate `[Ep_N/C;σ_C/m^2]` (Electric Field; Linear charge density) at position `[d=5_mm]` above a square plate of width `[L=8_cm]` and surface `A=L^2` where `d << L` when `d < L/10` is verified:
```rpl
L=8_cm A='L^2' d=5_mm Q=6E-6_C  A=64_cm^2  εr=1
@ Expecting [ σ=9.375⁳⁻⁸ C/cm↑2 Ep=52 941 050.0044 N/C ]
if 'd < L/10' then
 'ROOT(ⒺE Field Infinite Plate;[σ;Ep];[1_C/cm^2;1_N/C])'
end
```

### Ohm’s Law & Power

* To calculate `[R_Ω;P_W]` (Resistance; Power) from 2 known variables:
```rpl
V=24_V  I=16_A
@ Expecting [ R=1.5 Ω P=384. W ]
'ROOT(ⒺOhm’s Law & Power;[R;P];[1_Ω;1_W])'
```

### Volt Divider

* To calculate `[V1_V]` (Voltage) from 3 known variables:
```rpl
R1=40_Ω  R2=10_Ω  V=100_V
@ Expecting [ V1=80. V ]
'ROOT(ⒺVolt Divider;[V1];[1_V])'
```

### Current Divider

* To calculate `[I1_A]` (Current) from 3 known variables:
```rpl
R1=10_Ω  R2=6_Ω  I=15_A
@ Expecting [ I1=5.625 A ]
'ROOT(ⒺCurrent Divider;[I1];[1_A])'
```

### Wire Resistance

* To calculate `[R_Ω]` (Resistance) from 3 known variables:
```rpl
ρ=0.0035_Ω*cm  L=50_cm  A=1_cm^2
@ Expecting [ R=0.175 Ω ]
'ROOT(ⒺWire Resistance;[R];[1_Ω])'
```

### Resistivity & Conductivity

The electrical resistivity `ρ` of most materials changes with temperature. If the temperature `T` does not vary too much, a linear approximation can be used around the reference point (`ρ0`; `T0`).

* To calculate `[ρ_(Ω*m);σ_(S/m)]` (Resistance) from 4 known variables:
```rpl
ρ0=1.68E-8_Ω*m  αT=4.04E-3_K^-1  T0=293,15_K  T=373,15_K
@ Expecting [ ρ=2.22297 6⁳⁻⁸ Ω·m σ=44 984 741.1758 S/m ]
'ROOT(ⒺResistivity & Conductivity;[ρ;σ];[1_(Ω*m);1_(S/m)])'
```

### Series & Parallel R

![Series & Parallel R](img/Series&ParallelR.bmp)

* To calculate `[Rs_Ω;Rp_Ω]` (Resistance) from 2 known variables:
```rpl
R1=2_Ω  R2=3_Ω
@ Expecting [ Rs=5. Ω Rp=1.2 Ω ]
'ROOT(ⒺSeries & Parallel R;[Rs;Rp];[1_Ω;1_Ω])'
```

### Series & Parallel C

![Series & Parallel C](img/Series&ParallelC.bmp)

* To calculate `[Cs_μF;Cp_μF]` (Capacitance) from 2 known variables:
```rpl
C1=2_μF  C2=3_μF
@ Expecting [ Cs=1.2 μF Cp=5. μF ]
'ROOT(ⒺSeries & Parallel C;[Cs;Cp];[1_μF;1_μF])'
```

### Series & Parallel L

![Series & Parallel L](img/Series&ParallelL.bmp)

* To calculate `[Ls_mH;Lp_mH]` (Inductance) from 2 known variables:
```rpl
L1=17_mH  L2=16.5_mH
@ Expecting [ Ls=33.5 mH Lp=8.37313 43283 6 mH ]
'ROOT(ⒺSeries & Parallel L;[Ls;Lp];[1_mH;1_mH])'
```

### Capacitive Energy

* To calculate `[V_V;q_μC]` (Potential; Charge) from 2 known variables:
```rpl
E=0.025_J  C=20_μF
@ Expecting [ V=50. V q=1 000. μC ]
'ROOT(ⒺCapacitive Energy;[V;q];[1_V;1_μC])'
```

### Volumic Density Electric Energy

* To calculate `[uE_(J/m^3)]` (Volumic Density Electric Energy) from 2 known variables:
```rpl
E=5_V/m  εr=1
@ Expecting [ uE=1.10677 34772⁳⁻¹⁰ J/m↑3 ]
'ROOT(ⒺVolumic Density Electric Energy;[uE];[1_(J/m^3)])'
```

### Inductive Energy

* To calculate `[I_A]` (Current) from 2 known variables:
```rpl
E=4_J L=15_mH
@ Expecting [ I=23.09401 07676 A ]
'ROOT(ⒺInductive Energy;[I];[1_A])'
```

### RLC Current Delay

![RLC Current Delay](img/RLCCurrentDelay.bmp)

* To calculate `[ω_r/s;φs_°;φp_°;XC_Ω;XL_Ω]` (Phases; Impedances) from 4 known variables:
```rpl
f=107_Hz  C=80_μF  L=20_mH  R=5_Ω
@ Expecting [ ω=672.30082 7868 r/s φs=-45.82915 71488 ° φp=-5.87715 65317 1 ° XC=18.59286 71836 Ω XL=13.44601 65574 Ω ]
'ROOT(ⒺRLC Current Delay;[ω;φs;φp;XC;XL];[1_r/s;1_°;1_°;1_Ω;1_Ω])'
```

### DC Capacitor Current

These equations approximate the dc current required to charge the voltage on a capacitor in a certain time interval.

* To calculate `[ΔV_V;Δt_s;tf_s]` (Voltage difference; Time difference; Final time) from 5 known variables:
```rpl
C=15_μF  Vi=2.3_V  Vf=3.2_V  I=10_A  ti=0_μs
@ Expecting [ ΔV=0.9 V Δt=1.35 μs tf=1.35 μs ]
'ROOT(ⒺDC Capacitor Current;[ΔV;Δt;tf];[1_V;1_μs;1_μs])'
```

### Capacitor Charge

* To calculate `[q_C]` (Voltage difference; Time difference; Final time) from 2 known variables:
```rpl
C=20_μF  V=100_V
@ Expecting [ q=2.⁳⁻³ C ]
'ROOT(ⒺCapacitor Charge;[q];[1_C])'
```

### DC Inductor Voltage

These equations approximate the dc voltage induced in an inductor by a change in current in a certain time interval.

![DC Inductor Voltage](img/DCInductorVoltage.bmp)

* To calculate `[ΔIL_A;ILf_A;tf_μs]` (Current difference; Final current; Final time) from 5 known variables:
```rpl
L=100_mH  V=52_V  Δt=32_μs  ILi=23_A  ti=0_μs
@ Expecting [ ΔIL=-1.664⁳⁻² A ILf=22.98336 A tf=32. μs ]
'ROOT(ⒺDC Inductor Voltage;[ΔIL;ILf;tf];[1_A;1_A;1_μs])'
```

### RC Transient

![RC Transient](img/RCTransient.bmp)

* To calculate `[V_V]` (Voltage) from 5 known variables:
```rpl
Vi=0_V  C=50_μF  Vf=10_V  R=100_Ω  t=2_ms
@ Expecting [ V=3.29679 95396 4 V ]
'ROOT(ⒺRC Transient;[V];[1_V])'
```

### RL Transient

![RL Transient](img/RLTransient.bmp)

* To calculate `[I_A]` (Current) from 5 known variables:
```rpl
Vi=0_V  Vf=5_V  R=50_Ω  L=50_mH  t=75_μs
@ Expecting [ I=-2.⁳⁻²³ A ]
@ Failing [ I=7.22565 13671 4⁳⁻³ A ]
@ C#1 NOT OK. MSOLVER gives wrong value BUT algebraic OK: '1/(R_Ω)*((Vf_V)-((Vf_V)-(Vi_V))*EXP(-((t_μs)*(R_Ω)/(L_mH))))'=7.22565 13671 4⁳⁻³ A 
'ROOT([(I_A)=1/(R_Ω)*((Vf_V)-((Vf_V)-(Vi_V))*EXP(-((t_μs)*(R_Ω)/(L_mH))))];[I];[1_A])'
'ROOT(ⒺRL Transient;[I];[1_A])'
```

### Resonant Frequency

* To calculate `[ω0;Qs;Qp;f0]` (Resonant pulsation; Parallel & Series quality factors; Resonant frequency) from 3 known variables:
```rpl
L=500_mH  C=8_μF  R=10_Ω
@ Failing [ ω0=500. r/s Qs=25. Qp=0.04 f0=79.57747 15459 Hz ]
@ C#2 NOT OK. MSOLVER & SOLVER Error: "Inconsistent units" but units are OK and it computes algebraically
'ROOT(ⒺResonant Frequency;[ω0;Qs;Qp;f0];[1_r/s;1;1;1_Hz])'
```

### Plate Capacitor

![Plate Capacitor](img/PlateCapacitor.bmp)

* To calculate `[d_cm;ΔV_V;Ein_(N/C);σ_(μC/m^2)]` (Distance; Voltage; Internal E field; Surface charge density) from 4 known variables:
```rpl
C=25_μF  εr=2.26  A=1_cm^2  Q=75_μC 
@ Failing [ d=8.00418 57871 2⁳⁻⁹ cm σ=750 000 μC/m↑2 Ein=3.74803 89383 6⁳¹⁰ V/m ΔV=3. V ]
@ C#3 NOT OK. MSOLVER: "NO solution?", OK if solve seperately for d, BUT not for σ:"Bad argument type" BUT strangely
@ if I give "1" as a trial value for σ, it works then it will also work seperately for Ein & ΔV BUT indicating SOLVE Error: "Sign Reversal"
'ROOT(ⒺPlate Capacitor;[d;σ;Ein;ΔV];[1_cm;1_(μC/m^2);1_(N/C);1_V])'
```

### Cylindrical Capacitor

![Cylindrical Capacitor](img/CylindricalCapacitor.bmp)

* To calculate `[C_μF;ΔV_V]` (Capacitance; Voltage) from 5 known variables:
```rpl
εr=1  Q=75_μC  ro=1_cm  ri=.999_cm  L=10_cm
@ Expecting [ C=5 560.46819 129 pF ΔV=13 488.07284 21 V ]
'ROOT(ⒺCylindrical Capacitor;[C;ΔV];[1_pF;1_V])'
```

### Solenoid Inductance

![Solenoid Inductance](img/SolenoidInductance.bmp)

* To calculate `[L_mH]` (Inductance) from 4 known variables:
```rpl
μr=2.5  n=40_1/cm  A=0.2_cm^2  h=3_cm
@ Expecting [ L=3.01592 89474 6⁳⁻² mH ]
'ROOT(ⒺSolenoid Inductance;[L];[1_mH])'
```

### Toroid Inductance

![Toroid Inductance](img/ToroidInductance.bmp)

* To calculate `[L_mH]` (Inductance) from 4 known variables:
```rpl
μr=1  N=5000  h=2_cm  ri=2_cm  ro=4_cm
@ Expecting [ L=69.31471 80562 mH ]
@ Error in ri input data of HP50gAUR.pdf
'ROOT(ⒺToroid Inductance;[L];[1_mH])'
```

### Sinusoidal Voltage

* To calculate `[ω_r/s;V_V]` (Pulsation, voltage) from 4 known variables:
```rpl
Vmax=110_V  t=30_μs  f=60_Hz  φ=15_°
@ Expecting [ ω=376.99111 8431 r/s V=29.66992 85671 V ]
'ROOT(ⒺSinusoidal Voltage;[ω;V];[1_r/s;1_V])'
```

### Sinusoidal Current

* To calculate `[I_A;f_Hz]` (Current; Frequency) from 4 known variables:
```rpl
t=3.2_s  Imax=10_A  ω=636_r/s  φ=30_°
@ Expecting [ f=101.22254 3806 Hz I=-0.28436 91656 8 A ]
'ROOT(ⒺSinusoidal Current;[f;I];[1_Hz;1_A])'
```

### Drift Speed & Current Density

* To calculate `[vd_m/s;J_(A/m^2);E_(V/m)]` (Drift speed; Current density; E field) from 5 known variables:
```rpl
@ A='Ⓒπ*(0.1_cm)^2'  ρ='Ⓒqe*n'  σ='Ⓒqe*n*40_(cm^2/(V*s))' This is the working version with numerical input values:
I=1_A  n=8.5e28_(m^-3) A=3.14159 26535 90E-2_cm↑2 ρ=1.36185 01389E10_C/m↑3 σ=54 474 005.556_S/m
@ Expecting [ vd=2.33733 41683 6⁳⁻⁵ m/s J=31.83098 86184 A/cm↑2 E=5.84333 54209⁳⁻³ V/m ]
'ROOT(ⒺDrift Speed & Current Density;[vd;J;E];[1_m/s;1_(A/cm^2);1_(V/m)])'

@ C#4 Second version where the preliminary calculations of input variables A, ρ & σ doesn't work "Invalid algebraic" BUT it should work !
I=1_A  n=8.5e28_(m^-3) r=0.1_cm A='Ⓒπ*r^2' ρ='Ⓒqe*n' σ='Ⓒqe*n*40_(cm^2/(V*s))'
@ Failing [ vd=2.33733 41683 6⁳⁻⁵ m/s J=31.83098 86184 A/cm↑2 E=5.84333 54209⁳⁻³ V/m ]
'ROOT(ⒺDrift Speed & Current Density;[vd;J;E];[1_m/s;1_(A/cm^2);1_(V/m)])'
```

### Electron & Hole Mobilities

In accordance with microscopic Ohm's law, the current density is proportional to the electric field. Holes and electrons therefore move at their average drift speeds during the mean free time between collisions. As long as the electric fields are not very high, the mobilities of holes and electrons are constant.

* To calculate `[μe_(cm^2/(V*s));μe_(cm^2/(V*s));J_(A/m^2);Je_(A/m^2);Jh_(A/m^2);σ_(S/m)]` (Electron & hole mobilities; Current densities; Conductivity) from 5 known variables:
```rpl
τc=4.09365 36801 40e-15_s meeff=1.09312 60456 68e-31 kg mheff=4.55469 18569 5e-31 kg nh=6.0e18_(m^-3) ne=1.04e19_(m^-3) E=6.0e-9_V/m
@ Expecting [ μe=60. cm↑2/(V·s) μh=14.4 cm↑2/(V·s) Je=5.99854 93176 9⁳⁻¹¹ A/m↑2 Jh=8.30568 36706 4⁳⁻¹² A/m↑2 J=6.82911 76847 6⁳⁻¹¹ A/m↑2 σ=1.13818 62807 9⁳⁻² S/m ]
'ROOT(ⒺElectron & Hole Mobilities;[μe;μh;Je;Jh;J;σ];[1_(cm^2/(V*s));1_(cm^2/(V*s));1_(A/m^2);1_(A/m^2);1_(A/m^2);1_(S/m)])'
```

## Fluids

The 30 variables in the Fluids section are:

* `ε`: Roughness (dim.: length)
* `μ`: Dynamic viscosity (dim.: force·time/volume, in SI: N·s/m^3)
* `ρ`: Density (dim.: mass/volume, in SI: kg/m^3)
* `ΔP`: Pressure change (dim.: force/area, in SI: pascal, Pa)
* `Δy`: Height change
* `ΣK`: Total fitting coefficients
* `A`: Cross-sectional area
* `A1`: Initial cross-sectional area
* `A2`: Final cross-sectional area
* `D`: Diameter
* `D1`: Initial diameter
* `D2`: Final diameter
* `f`: Fanning friction factor
* `h`: Depth relative to P0 reference depth
* `hL`: Head loss (dim.: area/time^2, in SI: m^2/s^2)
* `L`: Length
* `M`: Mass flow rate (dim.: mass/time, in SI: kg/s)
* `n`: Kinematic viscosity (dim.: force·time/volume, in SI: N·s/m^3)
* `P`: Pressure at h (dim.: force/area, in SI: pascal, Pa)
* `P0`: Reference pressure (dim.: force/area, in SI: pascal, Pa)
* `P1`: Initial pressure (dim.: force/area, in SI: pascal, Pa)
* `P2`: Final pressure (dim.: force/area, in SI: pascal, Pa)
* `Q`: Volume flow rate (dim.: volume/time, in SI: m^3/s)
* `Reynolds`: Reynolds number
* `v1`: Initial velocitie
* `v2`: Final velocitie
* `vavg`: Average velocity
* `W`: Power input (dim.: energy/time, in SI: watt, W)
* `y1`: Initial height
* `y2`: Final height

### Pressure At Depth

This equation describes hydrostatic pressure for an incompressible fluid. Depth `h` is positive downwards from the reference.

![Pressure At Depth](img/PressureAtDepth.bmp)

* To calculate `[P_kPa]` (Pressure) from 3 known variables:
```rpl
h=100_m  ρ=1025.1817_kg/m^3  P0=1_atm
@ Expecting [ P=1 106.68481 183 kPa ]
'ROOT(ⒺPressure at Depth;[P];[1_kPa])'
```

### Bernoulli Equation

These equations represent the streamlined flow of an incompressible fluid.

![Bernoulli Equation](img/BernoulliEquation.bmp)

* To calculate `[A1_in^2;ΔP_psi;Δy_ft;Q_ft^3/min;M_lb/min;v2_ft/s;A2_in^2;D2_in]` (Various hydrodynamic parameters) from 7 known variables:
```rpl
P2=25_psi  P1=75_psi  y2=35_ft  y1=0_ft  D1=18_in  ρ=64_lb/ft^3  v1=100_ft/s
@ Failing [ A1=254.46900 4941 in↑2 ΔP=-50. psi Δy=35. ft Q=10 602.87520 59 ft↑3/min M=678 584.01317 5 lb/min v2=122.42131 1569 ft/s A2=207.86332 19 in↑2 D2=16.26836 81217 in ]
@ C#5 NOT OK. MSOLVER: "Inconsistent units", OK for A1;ΔP;Δy;Q;M solved one at a time in order NOT OK for v2;A2;D2 "Inconsistent units" while searching for each unknown. It seems that isolating v2 from eqn (1) doesn't work
'ROOT(ⒺBernoulli Equation;[A1;ΔP;Δy;Q;M;v2;A2;D2];[1_in^2;1_psi;1_ft;1_ft^3/min;1_lb/min;1_ft/s;1_in^2;1_in])'
```

Alternate présentation adding one more known value: `v2`
```rpl
P2=25_psi  P1=75_psi  y2=35_ft  y1=0_ft  D1=18_in  ρ=64_lb/ft^3  v1=100_ft/s v2=122.421311569_ft/s
@ Failing [ A1=254.46900 4941 in↑2 ΔP=-50. psi Δy=35. ft Q=10 602.87520 59 ft↑3/min M=678 584.01317 5 lb/min v2=122.42131 1569 ft/s A2=207.86332 19 in↑2 D2=16.26836 81217 in ]
@ C#5 NOT OK. MSOLVER: "Inconsistent units", OK for A1;ΔP;Δy;Q;M solved one at a time NOT OK for A2;D2 wrong values are obtained which is strange because eqn (7) should resolve for A2 BUT the value of Q change inadvertently ... why??
'ROOT(ⒺBernoulli Equation;[A1;ΔP;Δy;Q;M;v2;A2;D2];[1_in^2;1_psi;1_ft;1_ft^3/min;1_lb/min;1_ft/s;1_in^2;1_in])'
```

### Flow with Losses

These equations extend Bernoulli’s equation to include power input (or output) and head loss.

![Flow with Losses](img/FlowwithLosses.bmp)

* To calculate `[A1_in^2;ΔP_psi;Δy_ft;Q_ft^3/min;M_lb/min;v2_ft/s;A2_in^2;D2_in]` (Various hydrodynamic parameters) from 7 known variables:
```rpl
P2=30_psi  P1=65_psi  y2=100_ft  y1=0_ft  ρ=64_lb/ft^3  D1=24_in  hL=2.0_(ft/s)^2  W=25_hp  v1=100_ft/s
@ Expecting [ A1=452.38934 2117 in↑2 ΔP=-35. psi Δy=100. ft Q=18 849.55592 15 ft↑3/min M=1 206 371.57898 lb/min v2=93.12684 14502 ft/s A2=485.77760 7264 in↑2 D2=24.86988 66004 in ]
'ROOT(ⒺFlow with Losses;[A1;ΔP;Δy;Q;M;v2;A2;D2];[1_in^2;1_psi;1_ft;1_ft^3/min;1_lb/min;1_ft/s;1_in^2;1_in])'
```

### Flow In Full Pipes

These equations adapt Bernoulli’s equation for flow in a round, full pipe, including power input (or output) and frictional losses (with the fanning friction factor `f`).

![Flow In Full Pipes](img/FlowInFullPipes.bmp)

* To calculate `[ΔP_Psi;Δy_ft;A_in^2;n_ft^2/s;Q_ft^3/min;M_lb/min;W_hp;Reynolds]` (Various hydrodynamic parameters) from 7 known variables:
```rpl
ρ=62.4_lb/ft^3  D=12_in  vavg=8_ft/s  P2=15_psi  P1=20_psi  y2=40_ft  y1=0_ft  μ=0.00002_lbf*s/ft^2  ΣK=2.25  ε=0.02_in  L=250_ft
@ Failing [ ΔP=-5. psi Δy=40. ft A=113.0973 in↑2 n=1.0312 ft^2/s Q=376.991 ft↑3/min M=23 524.2358 lb/min W=25.8897 hp Reynolds=775 780.5 ]
@ C#6 NOT OK. MSOLVER: "EquationSolver error: Inconsistent units". This uses the new DANNING function.
'ROOT(ⒺFlow In Full Pipes;[ΔP;Δy;A;n;Q;M;W;Reynolds];[1_psi;1_ft;1_in^2;1_ft^2/s;1_ft^3/min;1_lb/min;1_hp;1])'
```

## Forces and Energy

The 37 variables in the Force and Energy section are:

* `α`: Angular acceleration (dim.: angle/time^2, in SI: r/s^2)
* `ω`: Angular velocity (dim.: angle/time, in SI: r/s))
* `ωi`: Initial angular velocitie (dim.: angle/time, in SI: r/s)
* `ωf`: Final angular velocitie (dim.: angle/time, in SI: r/s)
* `ρ`: Fluid density (dim.: mass/volume, in SI: kg/m^3)
* `τ`: Torque (dim.: force·length, in SI: N·m)
* `Θ`: Angular displacement (dim.: angle, in SI: r)
* `A`: Projected area relative to flow
* `ar`: Centripetal acceleration at `r`
* `at`: Tangential acceleration at `r`
* `Cd`: Drag coefficient
* `E`: Energy (dim.: force·length, in SI: joule, J)
* `F`: Force at `r` or `x`, or Spring force ([Hooke’s Law](#Hooke’s Law)), or Attractive
force ([Law of Gravitation](#Law of Gravitation)), or Drag force ([Drag force](#Drag force))
* `I`: Moment of inertia (dim.: mass·length^2, in SI: kg·m^2)
* `k`: Spring constant (dim.: force/length, in SI: N/m)
* `Ki`: Initial kinetic energie (dim.: force·length, in SI: joule, J)
* `Kf`: Final kinetic energie (dim.: force·length, in SI: joule, J)
* `m`: Mass
* `m1`: First mass
* `m2`: Second mass
* `N`: Rotational speed (dim.: turn/time, in SI: rpm)
* `Ni`: Initial rotational speed (dim.: turn/time, in SI: rpm)
* `Nf`: Final rotational speed (dim.: turn/time, in SI: rpm)
* `P`: Instantaneous power (dim.: energy/time, in SI: watt, W)
* `Pavg`: Average power (dim.: energy/time, in SI: watt, W)
* `r`: Radius from rotation axis, or Separation distance ([Law of Gravitation](#Law of Gravitation))
* `t`: Time
* `v`: Velocity
* `vf`: Final velocity
* `v1f`: First final velocity
* `v2f`: Second final velocity
* `vi`: Initial velocity
* `v1i`: Initial velocity of mass `m1`
* `W`: Work (dim.: force·length, in SI: joule, J)
* `x`: Displacement
* `UGi`: Initial gravitational potential energy (dim.: force·length, in SI: joule, J)
* `UGf`: Final gravitational potential energy (dim.: force·length, in SI: joule, J)

### Linear Mechanics

* To calculate `[F_lbf;Ki_ft*lbf;vf_ft/s;Kf_ft*lbf;W_ft*lbf;x_ft;Pavg_hp]` (Various dynamical variables) from 4 known variables:
```rpl
t=10_s  m=50_lb  a=12.5_ft/s^2  vi=0_ft/s
@ Expecting [ F=19.42559 38572 lbf Ki=-1.2⁳⁻²³ ft·lbf vf=125. ft/s Kf=12 140.99616 08 ft·lbf W=12 140.99616 08 ft·lbf x=625. ft Pavg=2.20745 38474 1 hp ]
@ Note Ki is approximately 0
'ROOT(ⒺLinear Mechanics;[F;Ki;vf;Kf;W;x;Pavg];[1_lbf;1_ft*lbf;1_ft/s;1_ft*lbf;1_ft*lbf;1_ft;1_hp])'
```

### Angular Mechanics

* To calculate `[τ_ft*lbf;Ki_ft*lbf;W_ft*lbf;Kf_ft*lbf;at_ft/s^2;Ni_rpm;ωf_r/min;t_min;Nf_rpm;Pavg_hp]` (Various dynamical variables) from 5 known variables:
```rpl
I=1750_lb*in^2  Θ=360_°  r=3.5_in  α=10.5_r/min^2  ωi=0_r/s
@ Expecting [ τ=1.10168 29849 6⁳⁻³ ft·lbf Ki=-2.⁳⁻²³ ft·lbf W=6.92207 83442 6⁳⁻³ ft·lbf Kf=6.92207 83442 6⁳⁻³ ft·lbf at=8.50694 44444 4⁳⁻⁴ ft/s↑2 Ni=-2.⁳⁻²³ rpm ωf=11.48681 38076 r/min t=1.09398 22673 9 min Nf=0.29096 43928 17 rpm Pavg=1.91739 80792 8⁳⁻⁷ hp ]
'ROOT(ⒺAngular Mechanics;[τ;Ki;W;Kf;at;Ni;ωf;t;Nf;Pavg];[1_ft*lbf;1_ft*lbf;1_ft*lbf;1_ft*lbf;1_ft/s^2;1_rpm;1_r/min;1_min;1_rpm;1_hp])'
```

### Centripetal Force

* To calculate `[ω_r/s;F_N;v_m/s;ar_m/s]` (Various dynamical variables) from 3 known variables:
```rpl
m=1_kg  r=5_cm  N=2000_Hz
@ Expecting [ ω=12 566.37061 44 r/s F=7 895 683.52087 N v=628.31853 0718 m/s ar=7 895 683.52087 m/s↑2 ]
'ROOT(ⒺCentripetal Force;[ω;F;v;ar];[1_r/s;1_N;1_m/s;1_m/s^2])'
```

### Hooke’s Law

The force is that exerted by the spring.

![Hooke’s Law](img/Missing name.bmp)

* To calculate `[F_lbf;W_ft*lbf]` (Force; Work) from 2 known variables:
```rpl
k=1725_lbf/in  x=1.25_in
@ Expecting [ F=-2 156.25 lbf W=-112.30468 75 lbf·ft ]
'ROOT(ⒺHooke’s Law;[F;W];[1_lbf;1_ft*lbf])'
```

### 1D Elastic Collisions

![1D Elastic Collisions](img/1DElasticCollisions.bmp)

* To calculate `[v1f_m/s;v2f_m/s]` (Final velocities of mass `m1` & `m2`) from 3 known variables:
```rpl
m1=10_kg  m2=25_kg  v1i=100_m/s
@ Expecting [ v1f=-42.85714 28571 m/s v2f=57.14285 71429 m/s ]
'ROOT(Ⓔ1D Elastic Collisions;[v1f;v2f];[1_m/s;1_m/s])'
```

### Drag Force

* To calculate `[F_N]` (Drag force) from 4 known variables:
```rpl
Cd=0.05  ρ=1000_kg/m^3  A=7.5E6_cm^2  v=35_m/s
@ Expecting [ F=22 968 750. N ]
'ROOT(ⒺDrag Force;[F];[1_N])'
```

### Gravitation Law

* To calculate `[F_N;UGf_J;UGi_J;W_J]` (Gravitational force, work, final & initial potential energy) from 5 known variables:
```rpl
m1=2E15_kg  m2=2E18_kg  r=1000000_km  ri=1000000_km  rf=5000000_km
@ Expecting [ F=266 903.6 N UGf=-5.33807 2⁳¹³ J UGi=-2.66903 6⁳¹⁴ J W=2.13522 88⁳¹⁴ J ]
'ROOT(ⒺGravitation Law;[F;UGf;UGi;W];[1_N;1_J;1_J;1_J])'
```

### Relativity Mass Energy

* To calculate `[E_J]` (Relativistic energy) from 1 known variable:
```rpl
m=9.1E-31_kg
@ Expecting [ E=8.17867 21247⁳⁻¹⁴ J ]
'ROOT(ⒺRelativity Mass Energy;[E];[1_J])'
```

## Gases

The 38 variables in the Gases section are:

* `λ`: Mean free path (dim.: length)
* `ρ`: Flow density (dim.: mass/volume, in SI: kg/m^3)
* `ρr`: Reduced state factor
* `ρri`: Initial reduced state factor
* `ρrf`: Final reduced state factor
* `ρ0`: Stagnation density (dim.: mass/volume, in SI: kg/m^3)
* `A`: Flow area
* `At`: Throat area
* `d`: Molecular diameter
* `k`: Specific heat ratio
* `M`: Mach number
* `m`: Mass
* `MW`: Molecular weight (dim.: mass/mole, in SI: g/gmol)
* `n`: Number of moles, or Polytropic constant ([Polytropic Processes](#Polytropic Processes))
* `P`: Pressure, or Flow pressure ([Isentropic Flow](#Isentropic Flow)) (dim.: force/area, in SI: pascal, Pa)
* `P0`: Stagnation pressure (dim.: force/area, in SI: pascal, Pa)
* `Pc`: Pseudocritical pressure (dim.: force/area, in SI: pascal, Pa)
* `Pi`: Initial pressure (dim.: force/area, in SI: pascal, Pa)
* `Pf`: Final pressure (dim.: force/area, in SI: pascal, Pa)
* `Pr`: Reduced pressure
* `Pri`: Initial reduced ressure
* `Prf`: Final reduced pressure
* `T`: Temperature, or Flow temperature ([Isentropic Flow](#Isentropic Flow))
* `T0`: Stagnation temperature
* `Tc`: Pseudocritical temperature
* `Tr`: Reduced temperature
* `Tri`: Initial reduced temperature
* `trf`: Final reduced temperature
* `Ti`: Initial temperature
* `Tf`: Final temperature
* `V`: Volume
* `Vi`: Initial volume
* `Vf`: Final volume
* `vrms`: Root-mean-square (rms) velocity
* `W`: Work (dim.: force·length, in SI: joule, J)
* `Z`: Gas compressibility correction factor
* `Zi`: Initial gas compressibility correction factor
* `Zf`: Final gas compressibility correction factor

### Ideal Gas

* To calculate `[n_mol;m_kg]` (Number of moles; Mass) from 4 known variables:
```rpl
T=16.85_°C  P=1_atm  V=25_l  MW=36_g/mol
@ Expecting [ n=1.05056 26661 23⁳⁰ mol m=3.78202 55980 42⁳⁻² kg ]
'ROOT(ⒺIdeal Gas;[n;m];[1_mol;1_kg])'
```

### Ideal Gas Law Change

* To calculate `[Vf_l]` (Volume final) from 5 known variables:
```rpl
Pi=1.5_kPa  Pf=1.5_kPa  Vi=2_l  Ti=100_°C  Tf=373.15_K
@ Expecting [ Vf=2. l ]
'ROOT(ⒺIdeal Gas Law Change;[Vf];[1_l])'
```

### Isothermal Expansion

These equations apply to an ideal gas.

* To calculate `[m_kg;W_J]` (Mass; Work) from 5 known variables:
```rpl
Vi=2_l  Vf=125_l  T=300_°C  n=0.25_mol  MW=64_g/mol
@ Expecting [ m=1.6⁳⁻² kg W=4 926.49415 89 J ]
'ROOT(ⒺIsothermal Expansion;[m;W];[1_kg;1_J])'
```

### Polytropic Processes

These equations describe a reversible pressure-volume change of an ideal gas such that `P·Vn` is constant. Special cases include isothermal processes (`n = 1`), isentropic processes (`n = k`, the specific heat ratio), and constant-pressure processes (`n = 0`).

* To calculate `[n_1;Tf_°F]` (Polytropic number; Final temperature) from 5 known variables:
```rpl
Pi=15_psi  Pf=35_psi  Vi=1_ft^3  Vf=0.50_ft^3  Ti=75_°F
@ Expecting [ n=1.22239 24213 4 Tf=346.54537 037 K ]
'ROOT(ⒺPolytropic Processes;[n;Tf];[1;1_K])'
```

### Isentropic Flow

![Isentropic Flow](img/IsentropicFlow.bmp)

The calculation differs at velocities below and above Mach 1. The Mach number is based on the speed of sound in the compressible fluid.

* To calculate `[P_kPa;ρ_kg/m^3;At_cm^2]` (Flow pressure; Flow density; Throat area) from 7 known variables:
```rpl
k=2  M=0.9  T0=26.85_°C  T=373.15_K  ρ0=100_kg/m^3  P0=100_kPa  A=1_cm^2
@ Expecting [ P=154.71213 6111 kPa ρ=124.38333 3333 kg/m↑3 At=0.99280 71853 34 cm↑2 ]
'ROOT(ⒺIsentropic Flow;[P;ρ;At];[1_kPa;1_kg/m^3;1_cm^2])'
```

### Real Gas Law

These equations adapt the ideal gas law to emulate real-gas behavior.

* To calculate `[Z_1;n_mol;m_kg]` (Gas compressibility correction factor; Number of mole; Mass) from 7 known variables:
```rpl
Pc=48_atm  Tc=298_K  P=5_kPa  V=10_l  MW=64_g/mol  T=348.15_K
@ Failing [ Z=0.99977 57972 690 n=1.72768 40576 1⁳⁻² mol m=1.10571 77968 7⁳⁻³ kg ]
@ C#7 NOT OK. MSOLVER: "Inconsistent units". SOLVE for Z alone doesn't work (see below). If I substitute the Z value, SOLVE works for n & m
'ROOT(ⒺReal Gas Law;[Z,n;m];[1;1_mol;1_kg])'
```
```rpl
@ Let me examine the computation of Z alone in its own implicit expression
Pc=48_atm  Tc=298_K  P=5_kPa  T=348.15_K
@ Failing [ Z=0.99977 57972 69 ]
@ Doesn't work. MSOLVER : "Inconsistent units". Very strange because if I substitute P,T,Pc & Tc values, it works, see next trial 
'ROOT([Z=1+(0.31506237-1.04670990/((T_K)/(Tc_K))-0.57832729/((T_K)/(Tc_K))^3)*(0.27*(((P_Pa)/(Pc_Pa))/(Z*((T_K)/(Tc_K)))))+(0.53530771-0.61232032/((T_K)/(Tc_K)))*(0.27*(((P_Pa)/(Pc_Pa))/(Z*((T_K)/(Tc_K)))))^2+0.61232032*0.10488813*(0.27*(((P_Pa)/(Pc_Pa))/(Z*((T_K)/(Tc_K)))))^5/((T_K)/(Tc_K))+0.68157001*(0.27*(((P_Pa)/(Pc_Pa))/(Z*((T_K)/(Tc_K)))))^2/((T_K)/(Tc_K))^3*(1+0.68446549*(0.27*(((P_Pa)/(Pc_Pa))/(Z*((T_K)/(Tc_K)))))^2)*exp(-0.68446549*(0.27*(((P_Pa)/(Pc_Pa))/(Z*((T_K)/(Tc_K)))))^2)];[Z];[1])'
```
```rpl
@ Here it works if I substitute the P,T,Pc & Tc values explicitely. Why this substitution seems to be requires ?
@ Expecting [ Z=0.99977 57972 69 ]
'ROOT([Z=1+(0.31506237-1.04670990/((348.15_K)/(298_K))-0.57832729/((348.15_K)/(298_K))^3)*(0.27*(((5_kPa)/(4863.6_kPa))/(Z*((348.15_K)/(298_K)))))+(0.53530771-0.61232032/((348.15_K)/(298_K)))*(0.27*(((5_kPa)/(4863.6_kPa))/(Z*((348.15_K)/(298_K)))))^2+0.61232032*0.10488813*(0.27*(((5_kPa)/(4863.6_kPa))/(Z*((348.15_K)/(298_K)))))^5/((348.15_K)/(298_K))+0.68157001*(0.27*(((5_kPa)/(4863.6_kPa))/(Z*((348.15_K)/(298_K)))))^2/((348.15_K)/(298_K))^3*(1+0.68446549*(0.27*(((5_kPa)/(4863.6_kPa))/(Z*((348.15_K)/(298_K)))))^2)*exp(-0.68446549*(0.27*(((5_kPa)/(4863.6_kPa))/(Z*((348.15_K)/(298_K)))))^2)];[Z];[1])'
```

### Real Gas State Change

This equation adapts the ideal gas state-change equation to emulate real-gas behavior.

* To calculate `[Zi_1;Zf_1;Vf_l]` (Initial & final gas compressibility correction factor; Final volume) from 7 known variables:
```rpl
Pc=48_atm  Pi=100_kPa  Pf=50_kPa  Ti=348.15_K  Tc=298_K  Vi=10_l  Tf=523.15_K
@ Failing [ Vf=30.17028 92973 l Zi=0.99550 62096 36 Zf=0.99938 68303 14 ]
@ C#7 NOT OK. MSOLVER: "Inconsistent units" even if I provide the value for Zi & Zf (see below)SOLVE for Vf : "Inconsistent units"
'ROOT(ⒺReal Gas State Change;[Zi;Zf;Vf];[1;1;Vf_l])'
```

Same comment as "Real Gas Law" apply here, In the following I substitute all the required values, otherwise this doesn't work
```RPL
@ SOLUTIONS with individual ROOT calls. Since it works the Zi & Zf function (in equation.cc) will be rewritten in closed form
@ Expecting [ Zi=0.99550 62096 36 ]
'ROOT([Zi=1+(0.31506237-1.04670990/((348.15_K)/(298_K))-0.57832729/((348.15_K)/(298_K))^3)*(0.27*(((100_kPa)/(4863.6_kPa))/(Zi*((348.15_K)/(298_K)))))+(0.53530771-0.61232032/((348.15_K)/(298_K)))*(0.27*(((100_kPa)/(4863.6_kPa))/(Zi*((348.15_K)/(298_K)))))^2+0.61232032*0.10488813*(0.27*(((100_kPa)/(4863.6_kPa))/(Zi*((348.15_K)/(298_K)))))^5/((348.15_K)/(298_K))+0.68157001*(0.27*(((100_kPa)/(4863.6_kPa))/(Zi*((348.15_K)/(298_K)))))^2/((348.15_K)/(298_K))^3*(1+0.68446549*(0.27*(((100_kPa)/(4863.6_kPa))/(Zi*((348.15_K)/(298_K)))))^2)*exp(-0.68446549*(0.27*(((100_kPa)/(4863.6_kPa))/(Zi*((348.15_K)/(298_K)))))^2)];[Zi];[1])'

@ Expecting [ Zf=0.99938 68303 14 ]
'ROOT([Zf=1+(0.31506237-1.04670990/((523.15_K)/(298_K))-0.57832729/((523.15_K)/(298_K))^3)*(0.27*(((50_kPa)/(4863.6_kPa))/(Zf*((523.15_K)/(298_K)))))+(0.53530771-0.61232032/((523.15_K)/(298_K)))*(0.27*(((50_kPa)/(4863.6_kPa))/(Zf*((523.15_K)/(298_K)))))^2+0.61232032*0.10488813*(0.27*(((50_kPa)/(4863.6_kPa))/(Zf*((523.15_K)/(298_K)))))^5/((523.15_K)/(298_K))+0.68157001*(0.27*(((50_kPa)/(4863.6_kPa))/(Zf*((523.15_K)/(298_K)))))^2/((523.15_K)/(298_K))^3*(1+0.68446549*(0.27*(((50_kPa)/(4863.6_kPa))/(Zf*((523.15_K)/(298_K)))))^2)*exp(-0.68446549*(0.27*(((50_kPa)/(4863.6_kPa))/(Zf*((523.15_K)/(298_K)))))^2)];[Zf];[1])'
@ => Zi/Zf=0.99611 69983 83  which is the required answer if I deduce it from the HP50g results
```

### Kinetic Theory

These equations describe properties of an ideal gas.

* To calculate `[vrms_m/s;n_mol;m_kg;λ_nm]` (Root-mean-square velocity; Number of mole; Mean free path) from 7 known variables:
```rpl
P=100_kPa  V=2_l  T=300_K  MW=18_g/mol  d=2.5_nm
@ Failing [ vrms=644.76778 7657 m/s n=0.08018 11130 98 mol m=1.44326 00357 69⁳⁻³ kg λ=1.49163 44918 94⁳⁰ nm ]
@ Expecting [ vrms=644.76778 7657 m/s n=0.08018 11130 98 mol m=1.44326 00357 7⁳⁻³ kg λ=1 nm ]
@ C#8 NOT OK. MSOLVER calculates a wrong λ value, SOLVE only calculates separately the 3 first unknowns then the computation of λ is wrong
@ MSOLVER calculates wrong values everuwhere : [ vrms=1 388.08583 078 m/s n=0.0173 mol m=0.00031 14 kg λ=1 nm ]
'ROOT(ⒺKinetic Theory;[vrms;n;m;λ];[1_m/s;1_mol;1_kg;1_nm])'
```


## Heat transfer

The 31 variables in the Heat Transfer section are:

* `α`: Expansion coefficient (dim.: 1/temperature, in SI: K^-1)
* `δ`: Elongation (dim.: length)
* `λ1`: Lower wavelength limits (dim.: length)
* `λ2`: Upper wavelength limits (dim.: length)
* `λmax`: Wavelength of maximum emissive power (dim.: length)
* `ΔT`: Temperature difference
* `A`: Area
* `c`: Specific heat (dim.: energy/(mass·temperature), in SI: J/(kg·K))
* `eb12`: Emissive power in the range λ1 to λ2 (dim.: power/area, in SI: W/m^2)
* `eb`: Total emissive power (dim.: power/area, in SI: W/m^2)
* `f`: Fraction of emissive power in the range λ1 to λ2
* `h`: Convective heat-transfer coefficient (dim.: power/(area·temperature) in SI: W/(m^2·K))
* `h1`: Convective heat-transfer coefficient (dim.: power/(area·temperature) in SI: W/(m^2·K))
* `h3`: Convective heat-transfer coefficient (dim.: power/(area·temperature) in SI: W/(m^2·K))
* `k`: Thermal conductivity (dim.: power/(length·temperature) in SI: W/(m·K))
* `k1`: Thermal conductivity (dim.: power/(length·temperature) in SI: W/(m·K))
* `k2`: Thermal conductivity (dim.: power/(length·temperature) in SI: W/(m·K))
* `k3`: Thermal conductivity (dim.: power/(length·temperature) in SI: W/(m·K))
* `L`: Length
* `L1`: Length
* `L2`: Length
* `L3`: Length
* `m`: Mass
* `Q`: Heat capacity (dim.: energy, in SI: joule, J)
* `qr`: Heat transfer rate (dim.: power=energy/time, in SI: watt, W)
* `T`: Temperature
* `Tc`: Cold surface temperature ([Conduction^](#Conduction)), or Cold fluid temperature ([Convection](1Convection))
* `Th`: Hot surface temperature, or Hot fluid temperature ([Conduction + Convection](#Conduction + Convection))
* `Ti`: Initial temperature
* `Tf`: Final temperature
* `U`: Overall heat transfer coefficient (dim.: power/(area·temperature) in SI: W/(m^2·K))

### Heat Capacity

* To calculate `[Tf_°C;c_kJ/(kg*K)]` (Final temperature; Specific heat) from 4 known variables:
```rpl
ΔT=15_°C  Ti=0_°C  m=10_kg  Q=25_kJ
@ Expecting [ Tf=15. °C c=0.16666 66666 67 kJ/(kg·K) ]
'ROOT(ⒺHeat Capacity;[Tf;c];[1_°C;1_kJ/(kg*K)])'
```

### Thermal Expansion

![Thermal Expansion](img/ThermalExpansion.bmp)

* To calculate `[α_K^-1;Ti_°C]` (Expansion coefficient; Initial temperature) from 4 known variables:
```rpl
ΔT=15_°C  L=10_m  Tf=25_°C  δ=1_cm
@ Expecting [ α=6.66666 66666 7⁳⁻⁵ K⁻¹ Ti=10. °C ]
'ROOT(ⒺThermal Expansion;[α;Ti];[1_K^-1;1_°C])'
```

### Conduction

![Conduction](img/Conduction.bmp)

* To calculate `[qr_W;ΔT_°C]` (Heat transfer rate; Temperature difference) from 5 known variables:
```rpl
Tc=25_°C  Th=75_°C  A=12.5_m^2  L=1.5_cm  k=0.12_W/(m*K)
@ Expecting [ qr=5000 W ΔT=50_°C ]
'ROOT(ⒺConduction;[qr;ΔT];[1_W;1_°C])'
```

### Convection

![Convection](img/Convection.bmp)

* To calculate `[Th_°C;ΔT_°C]` (Hot surface temperature; Temperature difference) from 4 known variables:
```rpl
Tc=26.85_°C  A=200_m^2  h=0.005_W/(m^2*K)  qr=10_W
@ Expecting [ Th=36.85 °C ΔT=10. °C ]
'ROOT(ⒺConvection;[Th;ΔT];[1_°C;1_°C])'
```

### Conduction + Convection

If you have fewer than three layers, give the extra layers a zero thickness and any nonzero conductivity. The two temperatures are fluid temperatures – if instead you know a surface temperature, set the corresponding convective coefficient to 10^999999.

![Conduction + Convection](img/Conduction+Convection.bmp)

* To calculate `[Tc_°C;qr_W;U_W/(m^2*K)]` (Cold surface temperature; Heat transfer rate; Overall heat transfer coefficient) from 11 known variables:
```rpl
ΔT=35_°C  Th=55_°C  A=10_m^2  h1=0.05_W/(m^2*K)  h3=0.05_W/(m^2*K)  L1=3_cm  L2=5_cm  L3=3_cm  k1=0.1_W/(m*K)  k2=.5_W/(m*K)  k3=0.1_W/(m*K)
@ Expecting [ Tc=20. °C qr=8.59950 85995 1 W U=0.02457 00245 7 W/(m↑2·K) ]
'ROOT(ⒺConduction & Convection;[Tc;qr;U];[1_°C;1_W;1_W/(m^2*K)])'
```

### Black Body Radiation

F0λ(λ_m, T_K) is the black body emissive power Function which returns the fraction of total black-body emissive power at temperature `T_K` between wavelengths 0 and `λ_m`. It is the integral of the Planck distribution.

![Black Body Radiation](img/BlackBodyRadiation.bmp)

* To calculate `[λmax_nm;eb_W/m^2;f_1;eb12_W/m^2;q_W]` (Wavelength of maximal emission; Total emissive power; Fraction of emissive power between λ1 & λ2; Emissive power between λ1 & λ2; Heat transfer rate) from 4 known variables:
```rpl
T=1273,15_K  Tmax=1273,15_K  λ1=1000_nm  λ2=600_nm  A=1_cm^2
@ Failing [ λmax=2 276.05231 12 nm eb=148 984.27029 5 W/m↑2 f=3.60929 06376 4e-3 eb12=537.72753 1933 W/m↑2 q=14.89842 70295 W ]
@ C#9 NOT OK: MSOLVER: "Invalid function" due to the integration limits having units (see ISSUE #1307 & #1314), SOLVE computes eb with "Sign reversal", then SOLVE computes eb12 & q correctly
'ROOT(ⒺBlack Body Radiation;[λmax;eb;f;eb12;q];[1_nm;1_W/m^2;1;1_W/m^2;1_W])'
```
The following integral doesn't work:
T=1273,15_K  λ1=1000_nm  λ2=600_nm
'f=15/Ⓒπ^4*∫(Ⓒh*Ⓒc/((λ2_nm)*Ⓒk*(T_K));Ⓒh*Ⓒc/((λ1_nm)*Ⓒk*(T_K));X^3/expm1(X);X)'
But if I calculate the value of the integration limits and apply BASE, I obtain pure number and the integral works:
'Ⓒh*Ⓒc/((λ2_nm)*Ⓒk*(T_K))'=1.88347 62339 3⁳⁻⁸ m/nm and if I apply BASE I obtain 18.83476 23393
'Ⓒh*Ⓒc/((λ1_nm)*Ⓒk*(T_K))'=1.13008 57403 6⁳⁻⁸ m/nm and if I apply BASE I obtain 11.30085 74036
Then the integral computes
'f=15/Ⓒπ^4*∫(11.30085 74036;18.83476 23393;X^3/expm1(X);X)'=3.60929 06376 4e-3


## Magnetism

The 28 variables in the Magnetism section are:

* `α1`: Subtended internal left angle relative to the top ends of the solenoid
* `α2`: Subtended internal right angle relative to the top ends of the solenoid
* `θ1`: Subtended left angle relative to the top ends of the solenoid
* `θ2`: Subtended right angle relative to the ends of the wire
* `θ`: Angle between the line of the magnetic field and the speed of the moving charge
* `μr`: Relative permeability
* `B`: Magnetic field (dim.: mass/(time^2·current), in SI: tesla, T)
* `d`: Separation distance
* `Dpitch`: Pitch of the helicoidal motion (dim.: length)
* `Fba`: Force
* `fc`: Cyclotron frequency (dim.: time^-1, in SI: hertx, Hz)
* `I`: Current (dim.: charge/time, in SI: ampere, A)
* `Ia`: Current (dim.: charge/time, in SI: ampere, A)
* `Ib`: Current (dim.: charge/time, in SI: ampere, A)
* `L`: Length
* `m`: Mass
* `N`: Total number of turns
* `nl`: Number of turns per unit length (dim.: length^-1)
* `q`: Charge of the moving charge (dim.: charge, in SI: coulomb, C)
* `r`: Distance from center of wire
* `Rc`: Radius of the circular part of the motion
* `ri`: Inside radius of toroid
* `ro`: Outside radius of toroid
* `rw`: Radius of wire
* `uB`: Volumic density of magnetic energy (dim.: energy/volume, in SI: J/m^3)
* `v`: Speed of the moving charge
* `T`: Period (dim.: time)
* `VH`: Hall tension (dim.: energy/charge, in SI: volt, V)

#### Straight Wire Infinite

The magnetic field expression differs depending upon whether the point at `r` is inside or outside the wire of radius `rw` and the calculations are done accordingly. The expression for the magnetic field at the distance `r` is approximately valid if the distance is such that `r << L` and therefore also applies for a wire of finite length `L` providing that `r < L/10` (see the following example 2 compared to example 2 of "Straight Wire Finite"). Note that if an electric current passes through a straight wire, one must use the following right-hand rule te determine the direction of the `B` field: when the thumb is pointed in the direction of conventional current (from positive to negative), the curled fingers will then point in the direction of the magnetic field (see fig.).


![Straight Wire Infinite](img/StraightWireInfinite.bmp)

* **Example 1.** Inside the wire, to calculate `[B_T]` (Magnetic field) from 3 known variables:
```rpl
μr=1  rw=0.25_cm  r=0.2_cm  I=25_A
@ Expecting [ B=0.0016 T ]
'ROOT(ⒺStraight Wire Infinite;[B];[1_T])'
```

* **Example 2.** Outside the wire, to calculate `[B_T]` (Magnetic field) from 3 known variables:
```rpl
μr=1  rw=0.25_cm  r=5_cm  I=25_A
@ Expecting [ B=1.⁳⁻⁴ T ]
'ROOT(ⒺStraight Wire Infinite;[B];[1_T])'
```

The code below saves the reference value for comparison with the example 2 in ([B Field Finite Wire](#B Field Finite Wire)):
```rpl
@ Save the reference value for comparison below
B0=B
```

#### Straight Wire Finite

The expression for the magnetic field at the distance `r` depends on the subtended angles `θ1` and `θ2` relative to the ends of the wire of finite length `L`. The magnetic field expression differs depending upon whether the point at `r` is inside or outside the wire of radius `rw` and the calculations are done accordingly.

![B Field From Finite Wire](img/B Field From Finite Wire.bmp)

* **Example 1.** To calculate `[B_T]` (Magnetic field) from 6 known variables:
```rpl
μr=1_1  rw=0.25_cm  r=5_cm  I=25_A  θ1=30_°  θ2=150_°
@ Expecting [ 8.66025 40378 7⁳⁻⁵ T ]
'ROOT(ⒺStraight Wire Finite;[B];[1_T])'
```

* **Example 2.** When 'r << L' which means 'r < L/10', we can verify that the value of B for a infinite wire approximates the exact value calculated for a finite wire of length 'L'.
```rpl
L=3_m  μr=1  rw=0.25_cm  r=5_cm  I=25_A  θ1='atan(r/L/2)'  θ2='180_°-θ1'  
if 'r < L/10' then   
@ Expecting [ B=9.99965 27958 9⁳⁻⁵ T ]
'ROOT(ⒺStraight Wire Finite;[B];[1_T])'
end  
```
Verify relative difference under condition `5_cm << 3_m` with the example 2 of [Straight Wire Infinite](#Straight Wire Infinite)
```rpl
B0 B %Ch
@ Expecting [ -3.47204 10629⁳⁻³ ] 
@ % of relative difference 
```

#### Force Between Wires

The force between wires is positive for an attractive force (for currents having the same sign) and negative otherwise, corresponding to a repulsive force.

![Force Between Wires](img/ForceBetweenWires.bmp)

* To calculate `[Fba_N]` (Magnetic force) from 4 known variables:
```rpl
Ia=10_A  Ib=20_A  μr=1  L=50_cm  d=1_cm
@ Expecting [ Fba=2.00000 00000 1⁳⁻³ N ]
'ROOT(ⒺForce Between Wires;[Fba];[1_N])'
```

#### B Field In Infinite Solenoid

The expression for the magnetic field in the center is approximately valid if the radius of the solenoid < `L` and therefore also applies inside a solenoid of finite length `L`. The right-hand rule applies also here: when the fingers curl around the solenoid in the sense of the current, the thumb points in the direction of the magnetic field (see fig.).

![B Field In Infinite Solenoid](img/BFieldInInfiniteSolenoid.bmp)

* To calculate `[B_T]` (Magnetic field) from 3 known variables:
```rpl
μr=10  nl=5000_m^-1  I=1.25_A
@ Expecting [ B=0.07853 98163 4 T ]
'ROOT(ⒺB Field In Infinite Solenoid;[B];[1_T])'
```
The code below saves the reference value for comparison with the example 2 in ([B Field Finite Solenoid](#B Field Finite Solenoid)):
```rpl
@ Save the reference value for comparison below
B0=B
```

#### B Field In Finite Solenoid

The expression for the magnetic field in the center depends on the subtended internal angles `α1` and `α2` relative to the top ends of the solenoid of finite length `L`.

![B Field In Finite Solenoid](img/B Field In Finite Solenoid.bmp)

* **Example 1.** Inside the wire, to calculate `[B_T]` (Magnetic field) from 5 known variables:
```rpl
μr=10  nl=5_mm^-1  I=1.25_A  α1=150_°  α2=30_°
@ Expecting [ B=6.80174 76159⁳⁻² T ]  
'ROOT(ⒺB Field In Finite Solenoid;[B];[1_T])'
```

* **Example 2.** When `r << L` which means `r < L/10`, we can verify that the value of `B` for a infinite solenoid approximates the exact value calculated for a finite solenoid of length `L`.
   ```rpl
L=3_m  μr=10  r=10_cm  nl=5000_m^-1  I=1.25_A  α2='atan(r/L/2)'  α1='180_°-α2'  
if 'r < L/10' then   
@ Expecting [ B=7.85289 10304 2⁳⁻² T ]
'ROOT(ⒺB Field In Finite Solenoid;[B];[1_T])'
end  
```
Verify relative difference under condition 10_cm << 3_m with the example of ([B Field In Infinite Solenoid](#B Field In Infinite Solenoid))
```rpl
B0 B %Ch
@ Expecting [ -1.38859 9604⁳⁻² ]
@ % of relative difference 
```

#### B Field In Toroid

The magnetic field `B` is calculated in the center of the torroid. The right-hand rule applies also here: when the fingers curl around the outer circle of the torroid following the sense of the current, the thumb points in the direction of the central magnetic field (see fig.).

![B Field In Toroid](img/BFieldInToroid.bmp)

* To calculate `[B_T]` (Magnetic field) from 5 known variables:
```rpl
μr=10  N=50  ri=5_cm  ro=7_cm  I=10_A
@ Expecting [ B=0.01666 66666 67 T ]
'ROOT(ⒺB Field In Toroid;[B];[1_T])'
```

#### Hall Effect

When a conductor carrying a current (to the right) is suddenly plunged into a perpendicular magnetic field (towards the bottom), there is a deviation (towards the top) of the free electrons which unbalances the distribution of electric charges. This transfer of charges from top to bottom gives rise to an electric force which in return balances the magnetic force. And it is precisely the appearance of this new electric field which is revealed by the Hall voltage `VH` measured between the top and the bottom of the conductor (small multimeter on the right).

![Hall Effect](img/Hall Effect VH.bmp)

* To calculate `[VH_V]` (Hall tension) from 5 known variables:
```rpl
n=5e28_(1/m^3)  B=0.641_T  q=1.60217 6634e-19_C  L=2_mm  I=10_A
@ Expecting [ VH=4.00080 73167 3⁳⁻⁷ V ]
@ To be verified since an error in the equation has to be corrected in equation.cc
'ROOT(ⒺHall Effect;[VH];[1_V])'
```

#### Cyclotron Motion

Under the perpendicular magnetic field, the moving charge has a circular trajectory and turns at the cyclotron frequency with the rotation period `T`.

![Cyclotron Motion](img/Cyclotron_BW.bmp)

* To calculate `[Rc_m;fc_Hz;T_s]` (Radius of the circular path; Cyclotron frequency; Period) from 4 known variables:
```rpl
m=1.67262 19259 5e-27_kg  B=0.8_T  q=1.60217 6634e-19_C  v=4.6e7_m/s
@ Expecting [ Rc=0.60028 18834 15 m fc=12 196 149.1501 Hz T=8.19930 93696 5⁳⁻⁸ s ]
'ROOT(ⒺCyclotron Motion;[Rc;fc;T];[1_m;1_Hz;1_s])'
```

#### Helicoidal Motion

Under the magnetic field lines (at angle `θ` with the speed vector), the moving charge has an helicoidal trajectory of pitch `Dpitch`, radius `Rc` and period `T`.

![Helicoidal Motion](img/Helicoidal Motion BW.bmp)

* To calculate `[Rc_m;T_s;Dpitch_m]` (Radius of the circular path; Period, pitch of the helicoidal motion) from 4 known variables:
```rpl
m=1.67262 19259 5e-27_kg  B=0.8_T  q=1.60217 6634e-19_C  v=4.6e7_m/s  θ=30_°
@ Expecting [ Rc=0.30014 09417 08 m T=8.19930 93696 5⁳⁻⁸ s Dpitch=3.26637 26955 m ]
'ROOT(ⒺHelicoidal Motion;[Rc;T;Dpitch];[1_m;1_s;1_m])'
```

#### Volumic Density Magnetic Energy

* To calculate `[uB_(J/m^3)]` (Volumic density of magnetic energy) from 2 known variables:
```rpl
μr=3.0  B=0.8_T 
@ Expecting [ uB=84 882.63631 54 J/m↑3 ]
'ROOT(ⒺVolumic Density Magnetic Energy;[uB];[1_(J/m^3)])'
```


## Motion

The 38 variables in the Motion section are:

* `α`: Angular acceleration (dim.: angle/time^2, in SI: r/s^2)
* `ω`: Angular velocity ([Circular Motion](#Circular Motion)), or Angular velocity at `t` ([Angular Motion](#Angular Motion)) (dim.: angle/time, in SI: r/s)
* `ω0`: Initial angular velocity (dim.: angle/time, in SI: r/s)
* `ρ`: Fluid density (dim.: mass/volume, in SI: kg/m^3)
* `θ`: Angular position at `t` (dim.: angle)
* `θ0`: Initial angular position ([Angular Motion](#Angular Motion)), or Initial vertical angle ([Projectile Motion](#Projectile Motion))
* `φ`: Latitude (dim.: angle)
* `a`: Acceleration
* `Ah`: Projected horizontal area
* `ar`: Centripetal acceleration at `r`
* `Cd`: Drag coefficient
* `fr`: Fraction of the terminal velocity `vt`
* `gloc`: local gravitational acceleration of a planet or star (dim.: length/time^2)
* `gearth`: local gravitational acceleration on Earth (dim.: length/time^2)
* `h`: altitude (dim.: length)
* `hmax`: maximum height of the projectile Motion
* `m`: Mass
* `Mp`: Planet or star mass
* `N`: Rotational speed (dim.: time^-1, in SI: hertz, Hz)
* `R`: Horizontal range ([Projectile Motion](#Projectile Motion)), or Planet or Star radius ([Object in Free Fall](#Object in Free Fall)) & ([Escape and Orbital Velocities](#Escape and Orbital Velocities))
* `r`: Radius
* `rc`: Radius of circular motion
* `t`: Time
* `tf`: Time of flight of a projectile  
* `tfr`: time required to reach the fraction `fr` of the terminal velocity
* `v`: Velocity at `t` ([linear Motion](#linear Motion)), or Tangential velocity at `r` [Circular Motion¸](#Circular Motion)), or Falling velocity at time `t` [Terminal Velocity](#Terminal Velocity)) & ([Buoyancy & Terminal Velocity](#Buoyancy & Terminal Velocity))
* `v0`: Initial velocity
* `ve`: Escape velocity in a gravitational field
* `vcx`: Horizontal (component x) of velocity at `t`
* `vcy`: Vertical (component y) of velocity at `t`
* `vo`: Orbital velocity in a gravitational field
* `Vol`: Volume of the moving object
* `vt`: Terminal velocity reached in a vertical fall
* `x`: Horizontal position at `t`
* `x0`: Initial horizontal position
* `xfr`: Displacement during `tfr`
* `y`: Vertical position at `t`
* `y0`: Initial vertical position

#### Linear Motion

* To calculate `[a_m/s^2;v_m/s]` (Acceleration; Velocity at time `t`) from 4 known variables:
```rpl
x0=0_m  x=100_m  t=10_s  v0=1_m/s 
@ Expecting [ a=1.8 m/s↑2 v=19. m/s ]
'ROOT(ⒺLinear Motion;[a;v];[1_m/s^2;1_m/s])'
```

#### Object In Free Fall

By definition, an object in free fall only experiences local gravitational acceleration `gloc`. This depends on the mass of the star or planet and the distance `r` center to center (where we assume that the position is greater than the radius of the mass). For the Earth, we can calculate an approximate value `gearth` of the acceleration of gravity as a function of latitude `φ` and for an altitude `h` low compared to the Earth's radius (typically: a few thousand meters, valid in commercial aviation).

* To calculate `[t_s;v_ft/s;r_m;gearth_m/s^2]` (Time; Velocity at time `t`; Radius; Acceleration of gravity at latitude `φ` & altitude `h`) from 7 known variables:
```rpl
y0=1000_ft  y=0_ft  v0=0_ft/s  gloc=9.80665_m/s↑2  φ=45_°  h=1000_m  Mp=5.9722e24_kg
@ Expecting [ t=7.88428 18533 5 s v=-253.66926 7182 ft/s r=6 374 488.95425 m gearth=9.80321 00310 8 m/s↑2 ]
'ROOT(ⒺObject In Free Fall;[t;v;r;gearth];[1_s;1_ft/s;1_m;1_m/s^2])'
```

#### Projectile Motion

During the time of flight `tf`, the motion of a projectile follows a symetric parabole of horizontal range `R` and of maximum height `hmax`.

![Projectile Motion](img/Projectile_BW.bmp)

* To calculate `[R_ft;vcx_ft/s;vcy_ft/s;x_ft;y_ft;hmax_ft;tf_s]` (Range, Components `x` & `y` of velocity at time `t`; Positions `x` & `y` at time `t`; Maximum height; Time of flight) from 5 known variables:
```rpl
x0=0_ft  y0=0_ft  θ0=45_°  v0=200_ft/s  t=10_s
@ Expecting [ R=1 243.23800 686 ft vcx=141.42135 6237 ft/s vcy=-180.31912 9327 ft/s x=1 414.21356 237 ft y=-194.48886 5448 ft hmax=310.80950 1716 ft tf=8.79102 02528 1 s ]
'ROOT(ⒺProjectile Motion;[R;vcx;vcy;x;y;hmax;tf];[1_ft;1_ft/s;1_ft/s;1_ft;1_ft;1_ft;1_s])'
```

#### Angular Motion

* To calculate `[ω_r/min;Θ_°]` (Angular velocity at time `t`; Angular position at time `t`) from 4 known variables:
```rpl
θ0=0_°  ω0=0_r/min  α=1.5_r/min^2  t=30_s
@ Expecting [ ω=0.75 r/min θ=10.74295 86587 ° ]
'ROOT(ⒺAngular Motion;[ω;θ];[1_r/min;1_°])'
```

#### Uniform Circular Motion

* To calculate `[ω_r/min;ar_ft/s^2;N_rpm]` (Angular velocity; Centripetal acceleration at `r`; Rotational speed) from 2 known variables:
```rpl
rc=25_in  v=2500_ft/s
@ Expecting [ ω=72 000. r/min ar=3 000 000. ft/s↑2 N=11 459.15590 26 rpm ]
'ROOT(ⒺUniform Circular Motion;[ω;ar;N];[1_r/min;1_ft/s^2;1_rpm])'
```

#### Terminal Velocity

Terminal velocity is the maximum speed attainable by an object as it falls through a fluid like air for instance. It is reached when the sum of the increasing drag force is equal to the downward force of gravity acting on the object (neglecting buoyancy), leading to a zero net force at the resulting terminal velocity.

* **Example 1**. For a falling big mass, to calculate `[vt_ft/s;v_ft/s;tfr_s;xfr_ft]` (Terminal velocity; Velocity at time `t`; Time required to reach the fraction `fr` of `vt`; Displacement during `tfr`) from 6 known variables:
//  '(175.74722 3631_ft/s)*∫(0;10.00590 25332;TANH(t*Ⓒg/(175.74722 3631_ft/s)*(1_s));t)*(1_s)' is the working integral for xfr
  
```rpl
Cd=0.15  ρ=0.025_lb/ft^3  Ah=100000_in^2  m=1250_lb  t=5_s  fr=0.95
@ Expecting [ vt=175.74722 3631 ft/s v=127.18655 2185 ft/s tfr=0.21337 88142 9 s xfr=1 ft ]
@ Failing [ vt=175.74722 3631 ft/s v=127.18655 2185 ft/s tfr=10.00590 25332 s  xfr=1 117.39339 246 ft ]
@ C#10 MSOLVER: works fine for vt, v & tfr  (see second ROOT call below). BUT integration with units seems to work here ISSUE #1314 but gives 
@ wrong value
'ROOT(ⒺTerminal Velocity;[vt;v;tfr;xfr];[1_ft/s;1_ft/s;1_s;1_ft])'
'ROOT(ⒺTerminal Velocity;[vt;v;tfr];[1_ft/s;1_ft/s;1_s])'
```

* **Example 2**. For a human skydiving head first, to calculate `[vt_m/s;v_m/s;tfr_s;xfr_m]` (Terminal velocity; Velocity at time `t`; Time required to reach the fraction `fr` of `vt`; Displacement during `tfr`) from 6 known variables:
//  '(95.13182 74789_m/s)*∫(0;17.76964 17471;TANH(t*Ⓒg/(95.13182 74789_m/s)*(1_s));t)*(1_s)' is the working integral for xfr

```rpl
Cd=0.7  ρ=1.29_kg/m^3  Ah=0.18_m^2  m=75_kg  t=5_s  fr=0.95
@ Failing [ vt=95.13182 74789 m/s  v=45.10777 55851 m/s  tfr=17.76964 17471 s  xfr=1 074.15231 681 m ]
@ C#10 MSOLVER: works fine for vt, v & tfr (see second ROOT call below).BUT integration with units seems to work here ISSUE #1314 but gives
@ wrong value
'ROOT(ⒺTerminal Velocity;[vt;v;tfr;xfr];[1_m/s;1_m/s;1_s;1_m])'
'ROOT(ⒺTerminal Velocity;[vt;v;tfr];[1_m/s;1_m/s;1_s])'
```

#### Buoyancy & Terminal Velocity

Terminal velocity is the maximum speed attainable by an object as it falls through a fluid like air for instance. It is reached when the sum of the increasing drag force plus the buoyancy is equal to the downward force of gravity acting on the object, leading to a zero net force at the resulting terminal velocity. If the terminal velocity is found to be negative, the motion is upward because buoyancy dominates gravity.

* **Example 1**. For a golf ball falling in water, to calculate `[vt_m/s;v_m/s;tfr_s;xfr_m]` (Terminal velocity; Velocity at time `t`; Time required to reach the fraction `fr` of `vt`; Displacement during `tfr`) from 8 known variables:
//  '(0.29459 06011 51 m/s)*∫(0;5.50264 78343 1e-2;TANH(t*Ⓒg/ABS(0.29459 06011 51 m/s)*(1_s));t)*(1_s)' is the working integral for xfr

```rpl
//input data: Cd=0.5  ρ=1.0775_(g/cm^3)  ρf=1000_(kg/m^3)  d=4.282_cm  Ah='Ⓒπ*((d_cm)/2)^2'  Vol='4/3*Ⓒπ*((d_cm)/2)^3'  t=3e-2_s  fr=0.95
Cd=0.5  ρ=1077,5_(kg/m^3)  ρf=1000_(kg/m^3)  d=4.282_cm  Ah=14.40068 68745_cm↑2  Vol=41.10916 07978_cm↑3  t=3e-2_s  fr=0.95
@ Failing [ vt=0.29459 06011 51 m/s  v=0.22419 40616 41 m/s  tfr=5.50264 78343 1e-2 s  xfr=0.01030 03495 63 m ]
@ C#10 MSOLVE works fine for vt, v & tfr (see second ROOT call below). BUT integration with units seems to work here ISSUE #1314 but gives
@ wrong value AND when investigating by hand eqn (1) & input data, it doesn't compute "Invalid Algebraic"
'ROOT(ⒺBuoyancy & Terminal Velocity;[vt;v;tfr;xfr];[1_m/s;1_m/s;1_s;1_m])'
'ROOT(ⒺBuoyancy & Terminal Velocity;[vt;v;tfr];[1_m/s;1_m/s;1_s])'
```

* **Example 2**. For a CO2 bubble in a glass of champagne, to calculate `[vt_m/s;v_m/s;tfr_s;xfr_m]` (Terminal velocity; Velocity at time `t`; Time required to reach the fraction `fr` of `vt`; Displacement during `tfr`) from 8 known variables:
//  '(-1.14234 81034 5_m/s)*∫(0;0.21337 88142 91;TANH(t*Ⓒg/ABS(-1.14234 81034 5_m/s)*(1_s));t)*(1_s)' is the working integral for xfr

```rpl
//input data: Cd=0.01  ρ=1.98_(kg/m^3)  ρf=998_(kg/m^3)  d=0.1_cm  Ah='Ⓒπ*((d_cm)/2)^2'  Vol='4/3*Ⓒπ*((d_cm)/2)^3'  t=0.1_s  fr=0.95
Cd=0.01  ρ=1.98_(kg/m^3)  ρf=998_(kg/m^3)  d=0.1_cm  Ah=7.85398 16339 7e-3_cm↑2  Vol=5.23598 77559 8e-4_cm↑3  t=0.1_s  fr=0.95
@ Failing [ vt=-1.14234 81034 5 m/s  v=-0.79446 37698 69 m/s  tfr=0.21337 88142 91 s  xfr=-0.15488 56277 53 m ]
@ C#10 MSOLVE works fine for vt, v & tfr (see second ROOT call below). BUT integration with units seems to work here ISSUE #1314 but gives
@ wrong value AND when investigating by hand eqn (1) & input data, it doesn't compute "Invalid Algebraic"
'ROOT(ⒺBuoyancy & Terminal Velocity;[vt;v;tfr;xfr];[1_m/s;1_m/s;1_s;1_m])'
'ROOT(ⒺBuoyancy & Terminal Velocity;[vt;v;tfr];[1_m/s;1_m/s;1_s])'
```

#### Escape and Orbital Velocities

The escape velocity is the speed required to completely free oneself from the gravitational field of a star, planet, etc. It is defined as the initial speed allowing you to move away to infinity. The orbital velocity is the speed nneded to maintain a stable circular orbit in a gravitational field.

* To calculate `[ve_m/s;vo_m/s]` (Escape velocity, orbital velocity) from 2 known variables:
    "Escape & Orbital Velocity",  "{ "
```rpl
Mp=1.5E23_lb  R=5000_mi
@ Expecting [ ve=3 485.10855 649 ft/s vo=2 464.34389 346 ft/s ]
'ROOT(ⒺEscape & Orbital Velocity;[ve;vo];[1_ft/s;1_ft/s])'
```


## Optics

The 40 variables in the Optics section are:

* `Δα`: Phase difference between top and bottom rays of light from a single slit (dim.: angle, in SI: r)
* `Δφ`: Phase difference between two rays of light seperated by an angle θ (dim.: angle, in SI: r)
* `Δyint`: Distance between bright fringes on the interference screen
* `Δydiff`: Width of the central maximum of diffraction
* `λ`: Light wavelength
* `θ`: Angle between initial light polarisation direction and polarizer transmission axis ([Malus Law](#Malus Law)), or Angle subtended by two points separated by y on the screen and the middle of one slit ([2 Slits Young Interference](#2 Slits Young Interference)), or two slits ([One Slit Diffraction](#One Slit Diffraction))
* `θ0`: Acceptance angle to enter an optical fiber in the outer medium of refraction index `n0f`
* `θ1`: Angle of incidence in the medium of refraction index n1
* `θ2`: Angle of refraction in the medium of refraction index n2
* `θB`: Brewster angle
* `θc`: Critical angle
* `θr`: Rayleigh's criterion limiting
* `a`: Width of a diffraction slit
* `d`: Distance between two interference slits
* `E0`: Incident light electrtc field
* `f`: Focal length
* `fx0`: Frequency of the incident X-ray (dim.: time^-1, in SI: hertz, HZ)
* `fx`: Frequency of the transmitted X-ray (dim.: time^-1, in SI: hertz, HZ)
* `I`: Transmitted irradiance or polarized light radiance flux (dim.: power/area, in SI: W/m^2)
* `I0`: Incident irradiance or polarized light radiance flux (dim.: power/area, in SI: W/m^2)
* `Ix`: Transmitted irradiance or polarized X rays radiance flux (dim.: power/area, in SI: W/m^2)
* `Ix0`: Incident irradiance or polarized X rayx radiance flux (dim.: power/area, in SI: W/m^2)
* `m`: Magnification
* `n`: Index of refraction
* `n1`: Index of refraction of medium 1
* `n2`: Index of refraction of medium 2
* `nf0`: Index of refraction of the outer medium at the entry of an optic fiber
* `nf1`: Index of refraction of the core medum of an optic fiber
* `nf2`: Index of refraction of the cladding medium of an optic fiber
* `r`: Radius of curvature
* `r1`: Radius of curvature
* `r2`: Radius of curvature
* `u`: Distance to object
* `v`: Distance to image
* `v1`: Speed of light in medium of refraction index n1
* `v2`: Speed of light in medium of refraction index n2
* `vn`: Speed of light in medium of refraction index n
* `L`: distance to the observation screen or to the light detector
* `NA`: Numerical aperture of the optic fiber
* `y`: distance between two image points on the observation screen (perpendicular to the optical axis) or distance between two point source

For reflection and refraction problems, the focal length and radius of curvature are positive in the direction of the outgoing light (reflected or refracted). The object distance is positive in front of the surface. The image distance is positive in the direction of the outgoing light (reflected or refracted). The magnification is positive for an upright image.

#### Refraction Law

![Refraction Law](img/RefractionLaw.bmp)

* To calculate `[θ2_°;v1_m/s;v2_m/s]` (Refraction angle; Speed of light in media of refraction index `n1` & `n2`) from 3 known variables:
```rpl
n1=1  n2=1.333  θ1=45_°
@ Expecting [ θ2=32.03672 30399 ° v1=299 792 457.998 m/s v2=224 900 568.642 m/s ]
'ROOT(ⒺRefraction Law;[θ2;v1;v2];[1_°;1_m/s;1_m/s])'
```

#### Critical Angle

![Critical Angle](img/CriticalAngle.bmp)

* To calculate `[θc_°;v1_m/s;v2_m/s]` (Critical angle; Speed of light in media of refraction index `n1` & `n2`) from 2 known variables:
```rpl
n1=1  n2=1.5
@ Expecting [ θc=41.81031 48958 ° v1=299 792 457.998 m/s v2=199 861 638.66 m/s ]
'ROOT(ⒺCritical Angle;[θc;v1;v2];[1_°;1_m/s;1_m/s])'
```

#### Fiber Optic

![Fiber Optic](img/Fiber Optic BW.bmp)

* To calculate `[vf0_m/s;vf1_m/s;vf2_m/s;NA;θ0_°]` (Speed of light in media of refraction index `n1f` & `n2f`; Numerical aperture; Acceptance angle) from 3 known variables:
```rpl
nf0=1.2  nf1=1.5  nf2=1.45
@ Expecting [ vf0=249 827 048.333 m/s vf1=199 861 638.66 m/s vf2=206 753 419.294 m/s NA=0.32004 77394 95 θ0=18.66581 19909 ° ]
'ROOT(ⒺFiber Optic;[vf0;vf1;vf2;NA;θ0];[1_m/s;1_m/s;1_m/s;1;1_°])'

#### Brewster’s Law

The Brewster angle is the angle of incidence at which the reflected wave is completely polarized.

![Brewster’s Law](img/Brewster’sLaw.bmp)

* To calculate `[θB_°;θ2_°;v1_m/s;v2_m/s]` (Brewster input angle; Refraction angle; Speed of light in media of refraction index `n1` & `n2`) from 2 known variables:
```rpl
n1=1  n2=1.5
@ Faiing [ v1=299 792 458 m/s v2=199 861 638.667 m/s θB=56.30993 2474 ° θB=33.69006 7526 ° ]
@ C#11 NOT OK MSOLVER: "Inconsistent units". SOLVE computation of θ2 alone fails.
'ROOT(ⒺBrewster’s Law;[v1;v2;θB;θ2];[1_m/s;1_m/s;1_°;1_°])'
```

#### Spherical Reflection

![Spherical Reflection](img/SphericalReflection.bmp)

* To calculate `[m;f_cm]` (Magnification; Focal distance) from 3 known variables:
```rpl
u=10_cm  v=300_cm  r=19.35_cm
@ Expecting [ m=-30. f=9.67741 93548 4 cm ]
'ROOT(ⒺSpherical Reflection;[m;f];[1;1_cm])'
```

#### Spherical Refraction

![Spherical Refraction](img/SphericalRefraction.bmp)

* To calculate `[n2;v1_m/s;v2_m/s]` (Refraction index in medium 2, speed of light in media of refraction index `n1`, `n2`) from 4 known variables:
```rpl
u=8_cm  v=12_cm  r=2_cm  n1=1
@ Expecting [ n2=1.5 v1=299 792 457.998 m/s v2=199 861 638.66 m/s ]
'ROOT(ⒺSpherical Refraction;[n2;v1;v2];[1;1_m/s;1_m/s])'
```

#### Thin Lens

![Thin Lens](img/ThinLens.bmp)

* To calculate `[f_cm;v_cm;m;vn_m/s]` (Focal distance, distance to image, magnification, speed of light in media of refraction index `n`) from 4 known variables:
```rpl
r1=5_cm  r2=20_cm  n=1.5  u=50_cm
@ Expecting [ f=13.33333 33333 cm v=18.18181 81818 cm m=-0.36363 63636 36 vn=199 861 638.66 m/s ]
'ROOT(ⒺThin Lens;[f;v;m;vn];[1_cm;1_cm;1;1_m/s])'
```

#### Rayleigh’s Criterion  

![Rayleigh’s Criterion](img/Rayleigh's Criterion_BW.bmp)

* To calculate `[θr_°;y_m]` (Rayleigh’s criterion angle; Distance between two point source) from 3 known variables:
```rpl
λ=550_nm  d=9_mm  L=18.7_km
@ Failing [ θr=4.27056 28265e-3 °  y=1.39381 16503 9 m ]
@ c#12 NOT OK MSOLVER: "Inconsistent units". SOLVE computation of θr & y alone fails: "Inconsistent units"
'ROOT(ⒺRayleigh’s Criterion;[θr;y];[1_°;1_m])'
```

#### Malus Law

If lineraly polarized light is incident on a perfect linear polarizer the transmitted light is the component at angle `θ` between the light polarisation direction and the polarizer transmission axis. The Malus law is given in terms of light irradiances. A relavistic version of the laws applies for X rays and more energetic electromagnetic radiations (with loss up to 10% in irradiance). The decrease in frequency (`fx < fx0`) and therefore in energy (`h·fx`) of a transmitted photon is due to the movement of the interacting electron of the polarizer (Compton scattering).

![Malus Law](img/Malus Law BW.bmp)

* To calculate `[I_(W/m^2);Ix_(W/m^2),E0_V/m]` (Polarized light radiance flux; Polarized radiance flux of emitted Xrays; Electric field) from 5 known variables:
```rpl
θ=30_°  I0=10_(W/m^2)  fx0=3e17_Hz  fx=2.7e17_Hz  I0x=0.1_(W/m^2)
@ Failing [ I=7.5 W/m↑2  Ix=0.06751 63889 34 W/m↑2   E0=86.80210 98203 V/m ]
@ c#13 NOT OK MSOLVER: "Inconsistent units". Expression of Ex needs to be corrected. But if I try with the correct expression 'I0x_(W/m^2)*((fx_Hz)/(fx0_Hz))*(1+Ⓒλc/Ⓒc*((fx0_Hz)-(fx_Hz)))*(COS(θ_°))^2' I get "Internal error, please report".
'ROOT(ⒺMalus Law;[I;Ix;E0];[1_(W/m^2);1_(W/m^2);1_V/m])'
```

#### 2 Slits Young Interference

* To calculate `[Δφ_r;I_(W/m^2);y_m;Δyint_m]` (Phase difference between two rays of light seperated by an angle θ; Irradiance of emitted light; Distance between two image points on the observation screen, Distance between bright fringes) from 5 known variables:
```rpl
L=2_m  d=800._μm  λ=600_nm  θ='ASIN(0.6*(λ_nm)/(d_μm))'  Imax=10_(W/m^2)
@ Expecting [ Δφ=3.76991 11843 1 r I=6.54508 49718 7 W/m↑2 y=9.00000 09112 5⁳⁻⁴ m Δyint=1.5⁳⁻³ m ]
'ROOT(Ⓔ2 Slits Young Interference;[Δφ;I;y;Δyint];[1_r;1_(W/m^2);1_m;1_m])'
```

#### One Slit Diffraction

* To calculate `[Δα_r;I_(W/m^2);y_m;Δydiff_m]` (Phase difference between top and bottom rays of light; Irradiance of emitted light; Distance between two image points on the observation screen; Width of the central maximum) from 5 known variables:
```rpl
L=3_m  a=1000._μm  λ=600_nm  θ='ASIN(0.3*(λ_nm)/(a_μm))'  Imax=10_(W/m^2)
@ Expecting [ Δα=1.88495 55921 5 r I=7.36839 72932 2 W/m↑2 y=5.40000 00874 8⁳⁻⁴ m Δydiff=3.6⁳⁻³ m ]
'ROOT(ⒺOne Slit Diffraction;[Δα;I;y;Δydiff];[1_r;1_(W/m^2);1_m;1_m])'
```


## Oscillations

The 25 variables in the Oscillations section are:

* `ω`: Angular frequency (dim.: angle/time, in SI: r!s)
* `ω0`: Natural angular frequency (dim.: angle#time, in SI: r!s)
* `ωu`: Underdamped angular frequency (dim.: angle/time, in SI: r!s)
* `γ`: Reduced damping coefficient (dim.: angle/time, in SI: r/s)
* `φ`: Phase angle
* `θ`: Cone angle
* `a`: Acceleration at `t`
* `E`: Total energy (kinetic plus potential)
* `f`: Frequency (dim.: time^-1; in SI: hertz, Hz)
* `G`: Shear modulus of elasticity (dim.: pressure, in SI: pascal, Pa)
* `h`: Cone height
* `I`: Moment of inertia (dim.: mass·area, in SI: kg·m^2)
* `J`: Polar moment of inertia (dim.: length^4; in SI: m^4)
* `k`: Spring constant (dim.: force/length, in SI: N/m)
* `L`: Length of pendulum
* `m`: Mass
* `Q`: Quality factor
* `t`: Time
* `T`: Period, Period dor small amplitude (Simple Pendulum)
* `Treal`: Real period for large amplitude (Simple Pendulum)
* `v`: Velocity at `t`
* `x`: Displacement at `t`
* `xm`: Displacement amplitude
* `xh`: Displacement amplitude of harmonic motion
* `xp`: Resulting displacement amplitude of driven & damped oscillations

#### Mass-Spring System

![Mass-Spring System](img/Mass-SpringSystem.bmp)

* To calculate `[ω_(r/s);T_s;f_Hz]` (Angular frequency; Oscillation period; Frequency) from 2 known variables:
```rpl
k=20_N/m  m=5_kg
@ Expecting [ ω=2. r/s T=3.14159 26535 9 s f=0.31830 98861 84 Hz ]
'ROOT(ⒺMass‐Spring System;[ω;T;f];[1_(r/s);1_s;1_Hz])'
```

#### Simple Pendulum

![Simple Pendulum](img/SimplePendulum.bmp)

* To calculate `[ω_(r/s);Treal_s;T_s;f_Hz]` (Angular frequency; Real period for large amplitude; Period for small amplitude; Frequency) from 2 known variables:
```rpl
L=15_cm  θmax=80_°
@ Failing [ ω=8.08564 57173 6 r/s  Treal=0.88361 42622 96 s   T=0.77707 89775 87 s  f=1.28687 04840 1 Hz ]
@ c#14 NOT OK MSOLVER & SOLVE: INFINITE LOOP & "Invalid function". However Treal can be calculated alone with θmax or (θmax_°) by the following 
@ sum: Treal='2*Ⓒπ*√((L_cm)/Ⓒg)*(Σ(x;0;5;((2·x)!÷((2↑x)·x!)²)²·sin((θmax_°)÷2)↑(2·x)))'
'ROOT(ⒺSimple Pendulum;[ω;Treal;T;f];[1_(r/s);1_s;1_s;1_Hz])'
```
```rpl
@ Without Treal, MSOLVER works fine for the remaining 3 unknowns:
L=15_cm  θmax=80_°
@ Expecting [ ω=8.08564 57173 6 r/s T=0.77707 89775 87 s f=1.28687 04840 1 Hz ]
'ROOT(ⒺSimple Pendulum;[ω;T;f];[1_(r/s);1_s;1_Hz])'
```

#### Conical Pendulum

![Conical Pendulum](img/ConicalPendulum.bmp)

* To calculate `[θ_°;ω_r/s;T_s;f_Hz]` (Cone angle; Angular frequency; Oscillation period; Frequency) from 2 known variables:
```rpl
L=25_cm  h=20_cm
@ Expecting [ θ=36.86989 76458 ° ω=7.00237 45972 3 r/s T=0.89729 35137 83 s f=1.11446 25305 3 Hz ]
'ROOT(ⒺConical Pendulum;[θ;ω;T;f];[1_°;1_r/s;1_s;1_Hz])'
```

#### Torsional Pendulum

![Torsional Pendulum](img/TorsionalPendulum.bmp)

* To calculate `[ω_r/s;T_s;f_Hz]` (Angular frequency; Oscillation period; Frequency) from 4 known variables:
```rpl
G=1000_kPa  J=17_mm^4  L=26_cm  I=50_kg*m^2
@ Expecting [ ω=1.14354 37497 9⁳⁻³ r/s T=5 494.48616 051 s f=1.82000 64042 2⁳⁻⁴ Hz ]
'ROOT(ⒺTorsional Pendulum;[ω;T;f];[1_r/s;1_s;1_Hz])'
```

#### Simple Harmonic

* To calculate `[x_cm;v_cm/s;a_m/s^2;m_kg;E_J]` (Displacement; Velocity & Acceleration at `t`; Mass; Total energy) from 5 known variables:
    "'(x_cm)=(xm_cm)*COS((ω0_(r/s))*(t_s)+(φ_°))' "
    "'(v_(m/s))=-(ω0_(r/s))/(1_r)*(xm_m)*SIN((ω0_(r/s))*(t_s)+(φ_°))' "
    "'(a_(m/s^2))=-(ω0_(r/s))^2/(1_r)^2*(xm_m)*COS((ω0_(r/s))*(t_s)+(φ_°))' "
    "'(ω0_(r/s))=(1_r)*√((k_(N/m))/(m_kg))' "
    "'(E_J)=(1/2)*(m_kg)*((ω0_(r/s))/(1_r)*(xm_m))^2' "
```rpl
xm=10_cm  ω0=15_r/s  φ=25_°  t=25_μs  k=10_N/m  
@ Failing [ x=9.06149 24146 7 cm  v=-63.44371 46156 cm/s  a=-2 038.83579 33 cm/s↑2  m=4.44444 44444 4⁳⁻² kg  E=1.125 J ]
@ C#15 NOT OK MSOLVER: "No solution?" Solve shows on sreen an error of sign in the acceleration eqn see ISSUE #1353
'ROOT(ⒺSimple Harmonic;[x;v;a;m;E];[1_cm;1_cm/s;1_m/s^2;1_kg;1_J])'
```

#### Underdamped Oscillations

We are considering here a damped mass-spring oscillator having the natural angular frequency `ω0`. The corresponding differential equation : `−k*x − b*dx/dt = m*d^2x/dt^2` describes the underdamped oscillations.

* To calculate `[m_kg;γ_(r/s);ωu_(r/s);x_cm;v_cm/s;a_m/s^2;E_J;Q]` (Mass; Reduced damping coefficient; Underdamped angular frequency; Displacement; Velocity & Acceleration at `t`; Mass; Total energy at `t`; Quality factor) from 6 known variables:
```rpl
xm=10_cm  ω0=15_r/s  φ=25_°  t=25_μs  k=10_N/m  b=0.2_(kg/s)  
@ Failing [ m=4.44444 44444 4⁳⁻² kg  γ=0.2 r/s  ωu=14.99966 6663 r/s  x=9.06146 97962 2 cm  v=-64.34829 19812 cm/s  a=-2 026.14705 038_cm/s↑2  E=0.24809 02514 79 J  Q=75.  ]
@ C#16 NOT OK MSOLVER: "No solution?". SOLVE for m alone doesn't work, howevwr algebraic expressions work and can be computed.
'ROOT(ⒺUnderdamped Oscillations;[m;γ;ωu;x;v;a;E;Q];[1_kg;1_(r/s);1_(r/s);1_cm;1_cm/s;1_m/s^2;1_J;1])'
```

The code below saves the reference value for comparison with the example in [Driven Damped Oscillations](#Driven Damped Oscillations):
```rpl
@ Save the reference value for comparison below
E0=E
```

#### Driven Damped Oscillations

We are considering here a damped mass-spring oscillator where the external driving force is of the form `Fdriving = Fd*cos(ω*t)` acting at the angular frequency `ω`. The corresponding differential equation : `−k*x − b*dx/dt + Fd*cos(ω*t) = m*d^2x/dt^2` describes the driven damped oscillations. When the driving frequency `ω` comes close to the natural frequency `ω0` this is the onset of resonance with amplitude increase and the total energy accumulates up to a possible catastrophy when the structure is overcome (see fig)

![Driven Damped Oscillations](img/Driven Damped Oscillations2_BW.bmp)

* To calculate `[m_kg;γ_(r/s);ωu_(r/s);φ_°;xp_m;x_cm;v_cm/s;a_m/s^2;E_J;Q]` (Mass; Reduced damping coefficient; Underdamped angular frequency; Phase angle; Resulting amplitude; Displacement; Velocity & Acceleration at `t`; Total energy at `t`; Quality factor) from 9 known variables which correspond to the values of the previous section:
```rpl
ω=14.9_r/s  ω0=15_r/s  θ=25_°  t=25_μs  k=10_N/m  b=0.2_(kg/s)  xh=10_cm  Fd=0.9_N
@ Failing [ m=4.44444 44444 4⁳⁻² kg  γ=0.2_r/s  ωu=14.99966 6663 r/s  φ=-44.90402 72598 °  xp=4.79694 88212 1 m  x=348.95077 8802 cm  v=4 979.30583 593 cm/s  a=-774.84786 538 m/s↑2  E=115.97995 9919 J  Q=75.  ]
@ C17 NOT OK MSOLVER: "No solution?". It is however perfectly determined as I was able to calculate. HINT#1: solving separately for γ:OK for m (in eqn 2) NOT OK: "No solution ?" which means unable to isolate: m='(k_(N/m))/((ω0_(r/s))/(1_r))^2'. HINT#2 solving separately for xp, Q, v, a & E leads to "Inconsistent units" even if I assign all intermediate results. But I insist on the fact that algebraic expression of each variable are OK and computes to expected correct values
'ROOT(ⒺDriven Damped Oscillations;[m;γ;ωu;φ;xp;x;v;a;E;Q];[1_kg;1_(r/s);1_(r/s);1_°;1_m;1_cm;1_cm/s;1_m/s^2;1_J;1])'
```

Verify relative difference with the total energy of the case ([Underdamped Oscillations](#Underdamped Oscillations))
```rpl
E0 E %Ch
@ Expecting [ 46 649.10006 64 ]
@ % of relative augmentation which illustrates the huge energy gain due to the driving force acting very near the resonance frequency.
```


## Plane Geometry

The 21 variables in the Plane Geometry section are:

* `β`: Central angle of polygon
* `θ`: Vertex angle of polygon
* `A`: Area
* `b`: Base length ([Rectangle](#Rectangle)) & ([Triangle](#Triangle)), or Length of semiaxis in x direction ([Ellipse](#Ellipse))
* `C`: Circumference
* `d`: Distance to rotation axis in y direction
* `h`: Height (Rectangle, Triangle), or Length of semiaxis in y direction (Ellipse)
* `I`: Moment of inertia about x axis (dim.: mass·area, in SI: kg·m^2)
* `Ix`: Moment of inertia about x axis (dim.: mass·area, in SI: kg·m^2)
* `Id`: Moment of inertia in x direction at distance d (dim.: mass·area, in SI: kg·m^2)
* `Iy`: Moment of inertia about y axis (dim.: mass·area, in SI: kg·m^2)
* `J`: Polar moment of inertia at centroid (dim.: mass·area, in SI: kg·m^2)
* `L`: Side length of polygon
* `n`: Number of sides
* `P`: Perimeter
* `r`: Radius
* `ri`: Inside radius
* `ro`: Outside radius
* `rs`: Distance to side of polygon
* `rv`: Distance to vertex of polygon
* `v`: Horizontal distance to vertex

#### Circle

![Circle](img/Circle.bmp)

* To calculate `[C_cm;A_cm^2;I_mm^4;J_mm^4;Id_mm^4]` (Circonference; Area; Moment of inertia about `x` axis; Polar moment of inertia at centroid; Moment of inertia in `x` direction at distance `d`) from 2 known variables:
```rpl
r=5_cm  d=1.5_cm
@ Expecting [ C=31.41592 65359 cm A=78.53981 63397 cm↑2 I=4 908 738.52123 mm↑4 J=9 817 477.04247 mm↑4 Id=6 675 884.38888 mm↑4 ]
'ROOT(ⒺCircle;[C;A;I;J;Id];[1_cm;1_cm^2;1_mm^4;1_mm^4;1_mm^4])'
```

#### Ellipse

![Ellipse](img/Ellipse.bmp)

* To calculate `[C_cm;A_cm^2;I_mm^4;J_mm^4;Id_mm^4]` (Circonference; Area; Moment of inertia about `x` axis; Polar moment of inertia at centroid; Moment of inertia in `x` direction at distance `d`) from 3 known variables:
```rpl
b=17.85_μm  h=78.9725_μin  d=.00000012_ft
@ Expecting [ C=7.98046 33593 6⁳⁻³ cm A=1.12485 79868⁳⁻⁶ cm↑2 I=1.13150 61302 6⁳⁻¹⁰ mm↑4 J=9.07327 72104 7⁳⁻⁹ mm↑4 Id=1.13301 09695 2⁳⁻¹⁰ mm↑4 ]
'ROOT(ⒺEllipse;[C;A;I;J;Id];[1_cm;1_cm^2;1_mm^4;1_mm^4;1_mm^4])'
```

#### Rectangle

![Rectangle](img/Rectangle.bmp)

* To calculate `[C_cm;A_cm^2;I_km^4;J_km^4;Id_km^4]` (Perimeter; Area; Moment of inertia about `x` axis; Polar moment of inertia at centroid; Moment of inertia in `x` direction at distance `d`) from 3 known variables:
```rpl
b=4_chain  h=7_rd  d=39.26_in
@ Expecting [ P=23 134.32 cm A=28 327 994.9568 cm↑2 I=2.92569 11916 2⁳⁻⁷ km↑4 J=1.82109 34968 2⁳⁻⁶ km↑4 Id=2.95386 09978 8⁳⁻⁷ km↑4 ]
'ROOT(ⒺRectangle;[P;A;I;J;Id];[1_cm;1_cm^2;1_km^4;1_km^4;1_km^4])'
```

#### Regular Polygon

![Regular Polygon](img/RegularPolygon.bmp)

* To calculate `[P_cm;A_cm^2;rs_cm;rv_cm;θ_°;β_°]` (Perimeter; Area; Distance to side; Distance to vertex of polygon; Vertex and central angles of polygon) from 2 known variables:
```rpl
n=8  L=0.5_yd
@ Expecting [ P=365.76 cm A=10 092.95006 19 cm↑2 rs=55.18892 20358 cm rv=59.73605 87541 cm θ=135. ° β=45. ° ]
'ROOT(ⒺRegular Polygon;[P;A;rs;rv;θ;β];[1_cm;1_cm^2;1_cm;1_cm;1_°;1_°])'
```

#### Circular Ring

![Circular Ring](img/CircularRing.bmp)

* To calculate `[A_cm^2;I_mm^4;J_mm^4;Id_mm^4]` (Area; Moment of inertia about `x` axis; Polar moment of inertia at centroid; Moment of inertia in `x` direction at distance `d`) from 3 known variables:
```rpl
ro=4_μm  ri=25_Å  d=0.1_mil
@ Expecting [ A=5.02654 62822 5e-7 cm↑2 I=2.01061 92983e-10 mm↑4 J=4.02123 85965 9e-10 mm↑4 Id=5.25354 58977 5e-10 mm↑4 ]
@ Note these answers are OK, the manual HP50g_AUR (p.5-46) is definitively in error (as calculated  directly)
'ROOT(ⒺCircular Ring;[A;I;J;Id];[1_cm^2;1_mm^4;1_mm^4;1_mm^4])'
```

#### Triangle

![Triangle](img/Triangle.bmp)

* To calculate `[b_in;A_in^2;Ix_in^4;Iy_in^4;J_in^4;Id_in^4]` (Base length; Area; Moment of inertia about `x` and `y` axis; Polar moment of inertia at centroid; Moment of inertia in `x` direction at distance `d`) from 4 known variables:
```rpl
h=4.33012781892_in  v=2.5_in  P=15_in  d=2_in
@ Expecting [ b=4.99999 90762 4 in A=10.82531 75473 in↑2 Ix=11.27637 66118 in↑4 Iy=11.27636 82785 in↑4 J=22.55274 48902 in↑4 Id=54.57764 6801 in↑4 ]
MSOLVE doesn't calculate anymore "Inconsistent units"
'ROOT(ⒺTriangle;[b;A;Ix;Iy;J;Id];[1_in;1_in^2;1_in^4;1_in^4;1_in^4;1_in^4])'
```


## Solid geometry
The 12 variables in the Solid Geometry section are:

* `A`: Total surface area
* `b`: Base length
* `d`: Distance to rotation axis in z direction
* `h`: Height in z direction ([Cone](#Cone)) & ([Cylinder](#Cylinder)), or Height in y direction ([Parallelepiped](#Parallelepiped))
* `I`: Moment of inertia about x axis (dim.: mass·area, in SI: kg·m^2)
* `Ixx`: Moment of inertia about x axis (dim.: mass·area, in SI: kg·m^2)
* `Id`: Moment of inertia in x direction at distance d (dim.: mass·area, in SI: kg·m^2)
* `Izz`: Moment of inertia about z axis (dim.: mass·area, in SI: kg·m^2)
* `m`: Mass
* `r`: Radius
* `t`: Thickness in z direction (dim.: length)
* `V`: Volume

#### Cone

![Cone](img/Cone.bmp)

* To calculate `[V_cm^3;A_cm^2;Ixx_kg*m^2;Izz_kg*m^2;Id_kg*m^2]` (Volume; Area; Moment of inertia about `x` axis & `z` axis; Moment of inertia in `x` direction at distance `d`) from 4 known variables:
```rpl
r=7_cm  h=12.5_cm  m=12.25_kg  d=3.5_cm
@ Expecting [ V=641.40850 0108 cm↑3 A=468.99530 2857 cm↑2 Ixx=0.01618 14843 75 kg·m↑2 Izz=1.80075⁳⁻² kg·m↑2 Id=3.11877 34375⁳⁻² kg·m↑2 ]
'ROOT(ⒺCone;[V;A;Ixx;Izz;Id];[1_cm^3;1_cm^2;1_kg*m^2;1_kg*m^2;1_kg*m^2])'
```

#### Cylinder

![Cylinder](img/Cylinder.bmp)

* To calculate `[V_in^3;A_in^2;Ixx_lb*in^2;Izz_lb*in^2;Id_lb*in^2]` (Volume; Area; Moment of inertia about `x` axis & `z` axis; Moment of inertia in `x` direction at distance `d`) from 4 known variables:
```rpl
r=8.5_in  h=65_in  m=12000_lb  d=2.5_in
@ Expecting [ V=14 753.70449 94 in↑3 A=3 925.42002 066 in↑2 Ixx=4 441 750. lb·in↑2 Izz=433 500. lb·in↑2 Id=4 516 750. lb·in↑2 ]
'ROOT(ⒺCylinder;[V;A;Ixx;Izz;Id];[1_in^3;1_in^2;1_lb*in^2;1_lb*in^2;1_lb*in^2])'
```

#### Parallelepiped

![Parallelepiped](img/Parallelepiped.bmp)

* To calculate `[V_in^3;A_in^2;I_lb*in^2;Id_lb*in^2]` (Volume; Area; Moment of inertia about `x` axis; Moment of inertia in `x` direction at distance `d`) from 4 known variables:
```rpl
b=36_in  h=12_in  t=72_in  m=83_lb  d=7_in
@ Expecting [ V=31 104. in↑3 A=7 776. in↑2 I=36 852. lb·in↑2 Id=40 919. lb·in↑2 ]
'ROOT(ⒺParallelepiped;[V;A;I;Id];[1_in^3;1_in^2;1_lb*in^2;1_lb*in^2])'
```

#### Sphere

![Sphere](img/Sphere.bmp)

* To calculate `[I_kg*m^2;r_cm;V_cm^3;A_cm^2]` (Moment of inertia about `x` axis; Radius; Volume; Area) from 3 known variables:
```rpl
d=14_cm  m=3.75_kg  Id=486.5_lb*in^2
@ Expecting [ I=6.88691 91393 3⁳⁻² kg·m↑2 r=8.43593 34656 2 in V=41 208.72679 31 cm↑3 A=5 769.57188 64 cm↑2 ]
'ROOT(ⒺSphere;[I;r;V;A];[1_kg*m^2;r_cm;1_cm^3;1_cm^2])'
```


## Solid State Devices
The 54 variables in the Solid State Devices section are:

* `αF`: Forward common-base current gain
* `αR`: Reverse common-base current gain
* `γ`: Body factor (dim.: (energy/charge)^0.5, in SI: V^0.5)
* `λ`: Modulation parameter (dim.: charge/energy, in SI: V^-1)
* `μn`: Electron mobility (dim.: speed/(electric field), in SI: m^2/(V·s))
* `φp`: Fermi potential (dim.: energy/charge, in SI: volt, V)
* `ΔL`: Length adjustment ([PN Step Junctions](#PN Step Junctions)), or Channel encroachment ([NMOS Transistors](#NMOS Transistors)) (dim.: length)
* `ΔW`: Width adjustment ([PN Step Junctions](#PN Step Junctions)), or Width contraction ([NMOS Transistors](#NMOS Transistors))
* `a`: Channel thickness
* `Aj`: Effective junction area
* `BV`: Breakdown voltage (dim.: energy/charge, in SI: volt, V)
* `Cj`: Junction capacitance per unit area (dim.: charge^2/(energy·area); in SI: F/m^2)
* `Cox`: Silicon dioxide capacitance per unit area (dim.: charge^2/(energy·area); in SI: F/m^2)
* `E1`: Breakdown-voltage field factor (dim.: force/charge, in SI: N/C=V/m)
* `Emax`: Maximum electric field (dim.: force/charge, in SI: N/C=V/m)
* `G0`: Channel conductance (dim.: resistance^-1, in SI: siemens, S)
* `gds`: Output conductance (dim.: resistance^-1, in SI: siemens, S)
* `gm`: Transconductance (dim.: resistance^-1, in SI: siemens, S)
* `I`: Diode current (dim.: charge/time, in SI: ampere, A)
* `IB`: Total base current (dim.: charge/time, in SI: ampere, A)
* `IC`: Total collector current (dim.: charge/time, in SI: ampere, A)
* `ICEO`: Collector current (collector-to-base open) (dim.: charge/time, in SI: ampere, A)
* `ICO`: Collector current (emitter-to-base open) (dim.: charge/time, in SI: ampere, A)
* `ICS`: Collector-to-base saturation current (dim.: charge/time, in SI: ampere, A)
* `ID, IDS`: Drain current (dim.: charge/time, in SI: ampere, A)
* `IE`: Total emitter current (dim.: charge/time, in SI: ampere, A)
* `IES`: Emitter-to-base saturation current (dim.: charge/time, in SI: ampere, A)
* `IS`: Transistor saturation current (dim.: charge/time, in SI: ampere, A)
* `J`: Current density (dim.: current/area, in SI: A/m^2)
* `Js`: Saturation current density (dim.: current/area, in SI: A/m^2)
* `L`: Drawn mask length ([PN Step Junctions](#PN Step Junctions)), or Drawn gate length ([NMOS Transistors](#NMOS Transistors)), or Channel length ([JFETs](#JFETs))
* `Le`: Effectives gate length
* `NA`: P-side doping ([PN Step Junctions](#PN Step Junctions)), or Substrate doping ([NMOS Transistors](#NMOS Transistors)) (dim.: 1/volume, in SI: m^-3)
* `ND`: N-side doping ([PN Step Junctions](#PN Step Junctions)), or N-channel doping ([JFETs](#JFETs)) (dim.: 1/volume, in SI: m^-3)
* `ni`: Silicon density (dim.: 1/volume, in SI: m^-3)
* `T`: Temperature
* `tox`: Gate silicon dioxide thickness
* `Va`: Applied voltage (dim.: energy/charge, in SI: volt, V)
* `VBC`: Base-to-collector voltage (dim.: energy/charge, in SI: volt, V)
* `VBE`: Base-to-emitter voltage (dim.: energy/charge, in SI: volt, V)
* `Vbi`: Built-in voltage (dim.: energy/charge, in SI: volt, V)
* `VBS`: Substrate voltage (dim.: energy/charge, in SI: volt, V)
* `VCEsat`: Collector-to-emitter saturation voltage (dim.: energy/charge, in SI: volt, V)
* `VDS`: Applied drain voltage (dim.: energy/charge, in SI: volt, V)
* `VDsat`: Saturation voltage (dim.: energy/charge, in SI: volt, V)
* `VGS`: Applied gate voltage (dim.: energy/charge, in SI: volt, V)
* `Vt`: Threshold voltage (dim.: energy/charge, in SI: volt, V)
* `Vt0`: Threshold voltage (at zero substrate voltage) (dim.: energy/charge, in SI: volt, V)
* `W`: Drawn mask width ([PN Step Junctions](#PN Step Junctions)), or Drawn width ([NMOS Transistors](#NMOS Transistors)), or Channel width ([JFETs](#JFETs))
* `We`: Effective width
* `xd`: Depletion-region width
* `xdmax`: Depletion-layer width
* `xj`: Junction depth

#### PN Step Junctions

These equations for a silicon PN-junction diode use a “two-sided step-junction” model–the doping density changes abruptly at the junction. The equation assume the current density is determined by minority carries injected across the depletion region and the PN junction is rectangular in its layout, The temperature should be between 77 and 500 K.

![PN Step Junctions](img/PNStepJunctions.bmp)

* To calculate `[ni_m^-3;Vbi_V;xd_μ;Cj_pF/cm^2;Emax_V/cm;BV_V;J_A/cm^2;Aj_cm^2;I_mA]` (Silicon density; Built-in voltage; Depletion-region width; Junction capacitance per unit area; Maximum electric field; Breakdown voltage; Current density; Effective junction area; Diode current) from 11 known variables:
```rpl
ND=1E22_cm^-3  NA=1E15_cm^-3  T=300_K  Js=1e-6_μA/cm^2  Va=-20_V  E1=3.3E5_V/cm  W=10_μ  ΔW=1_μ  L=10_μ  ΔL=1_μ  xj=2_μ
@ Expecting [ ni=9.64987 39813 5⁳¹⁵ (m↑3)⁻¹ Vbi=1.01379 86500 1 V xd=5.25726 60039 9 μ Cj=2 004.17545 83 pF/cm↑2 Emax=79 941.92659 86 V/cm BV=358.08260 5833 V J=-1.00000 00000 2⁳⁻¹² A/cm↑2 Aj=2.57097 33552 9⁳⁻⁶ cm↑2 I=-2.57097 338⁳⁻¹⁵ mA ]
'ROOT(ⒺPN Step Junctions;[ni;Vbi;xd;Cj;Emax;BV;J;Aj;I];[1_m^-3;1_V;1_μ;1_pF/cm^2;1_V/cm;1_V;1_A/cm^2;1_cm^2;1_mA])'
```

#### NMOS Transistor

These equations for a silicon NMOS transistor use a two-port network model. They include linear and nonlinear regions in the device characteristics and are based on a gradual-channel approximation (the electric fields in the direction of current flow are small compared to those perpendicular to the flow). The drain current and transconductance calculations differ depending on whether the transistor is in the linear, saturated, or cutoff region. The equations assume the physical geometry of the device is a rectangle, second-order length-parameter effects are negligible, shot-channel, hot-carrier, and velocity-saturation effects are negligible, and subthreshold currents are negligible.

![NMOS Transistor](img/NMOSTransistor.bmp)

* To calculate `[ni_(cm^-3);We_μ;Le_μ;Cox_pF/cm^2;γ_V^.5;φp_V;Vt_V;VDsat_V;IDS_mA;gds_S;gm_mA/V]` (Silicon density; Effective width; Effectives gate length; Silicon dioxide capacitance per unit area; Body factor; Fermi potential; Threshold voltage; Saturation voltage; Drain current; Output conductance; Transconductance) from 13 known variables:
```rpl
"'(γ_(V^(1/2)))='√((2*Ⓒεsi*Ⓒε0)*Ⓒqe*(NA_(cm^-3)))/(Cox_(pF/cm^2))' "
tox=700_Å  NA=1e15_1/cm^3  μn=600_(cm^2)/(V*s)  T=26.85_°C  Vt0=0.75_V  VGS=5_V  VBS=0_V  VDS=5_V  W=25_μ  ΔW=1_μ  L=4_μ  ΔL=0.75_μ  λ=0.05_1/V
@ Failing [ ni=1.03628 29628⁳²⁶ (m↑3)⁻¹ We=23. μ Le=2.5 μ Cox=49 330.47498 38 pF/cm↑2 γ=0.37247 98153 35 V↑(¹/₂) φp=-0.29855 54642 25 V
@ Vt=0.75 V VDsat=4.25 V IDS=IDS=2.97832 74271 5 mA gds=0.00014 89163 71 S gm=1.42391 28564 6 mA/V ]
@ C#17 NOT OK MSOLVE: "Error: constant ?". SOLVE seperately fails for γ & for Vt. BUT algebraic expressions are OK to compute.
'ROOT(ⒺNMOS Transistor;[ni;We;Le;Cox;γ;φp;Vt;VDsat;IDS;gds;gm];[1_m^-3;1_μ;1_μ;1_pF/cm^2;1_V^(1/2);1_V;1_V;1_V;1_mA;1_S;1_mA/V])'

```
#### Bipolar Transistors

These equations for an NPN silicon bipolar transistor are based on large-signal models developed by J.J. Ebers and J.L. Moll. The offset-voltage calculation differs depending on whether the transistor is saturated or not. The equations also include the special conditions when the emitter-base or collector-base junction is open, which are convenient for measuring transistor parameters.

![Bipolar Transistors](img/BipolarTransistors.bmp)

* To calculate `[VBE_V;IS_nA;ICO_nA;ICEO_nA;IE_mA;IC_mA;IB_mA;VCEsat_V]` (Base-to-emitter voltage; Transistor saturation current; Collector current (emitter-to-base open); Collector current (collector-to-base open); Total emitter current; Total collector current; Total base current; Collector-to-emitter saturation voltage) from 7 known variables:
```rpl
'((Ⓒk*(T_K)))/Ⓒqe*LN(((1+(IC_mA)/(IB_mA)*(1-αR)))/(αR*(1-(IC_mA)/(IB_mA)*((1-αF)/αF))))'
IES=1e-5_nA  ICS=2e-5_nA  T=26.85_°C  αF=0.98  αR=0.49  IC=1_mA  VBC=-10_V
Failing [ VBE=0.65531 72748 12 V IS=0.00000 98 nA ICO=0.00001 0396 nA ICEO=0.00051 98 nA IE=-1.02040 81632 5 mA IB=0.02040 81632 55 mA VCEsat=0.02962 86227 01 V VCE=...]
@ C#18 NOT OK MSOLVER: "Inconsistent units". Algebraic conforms to HP50g. SOLVE alone works but for IC=1. mA "Sign reversal".
@ Then SOLVE fails surprisingly for the IB (simple eqn). And fails again for VCEsat "Argument outside domain" which doesn't compute
'ROOT(ⒺBipolar Transistors;[VBE;IS;ICO;ICEO;IE;IC;IB;VCEsat];[1_V;1_nA;1_nA;1_nA;1_mA;1_mA;1_mA;1_V])'
```

#### JFETs

These equations for a silicon N-channel junction field-effect transistor (JFET) are based on the single-sided stepjunction
approximation, which assumes the gates are heavily doped compared to the channel doping,. The draincurrent
calculation differs depending on whether the gate-junction depletion-layer thickness is less than or greater
than the channel thickness. The equations assume the channel is uniformly doped and end effects (such as contact,
drain, and source resistances) are negligible.

![JFETs](img/JFETs.bmp)

* To calculate `[ni_(cm^-3);Vbi_V;xdmax_μ;G0_S;ID_mA;VDsat_V;Vt_V;gm_mA/V]` (Silicon density; Built-in voltage; Depletion-layer width; Channel conductance; Drain current; Saturation voltage; Threshold voltage; Transconductance) from 8 known variables:
```rpl
'(ID_mA)='(G0_S)*((VDS_V)-((2/3)*√((2*Ⓒεsi*Ⓒε0)/(Ⓒqe*(ND_(cm^-3))*(a_μ)^2)))*(((Vbi_V)-(VGS_V)+(VDS_V))^(3/2)-((Vbi_V)-(VGS_V))^(3/2)))'
ND=1e16_1/cm^3  W=6_μ  a=1_μ  L=2_μ  μn=1248_cm^2/(V*s)  VGS=-4_V  VDS=4_V  T=26.85_°C
@ Expecting [ ni=9.64987 39813 5⁳¹⁵ (m↑3)⁻¹ Vbi=0.35808 22815 91 V xdmax=1.04848 1973 μ G0=0.00059 98549 32 S ... ]
@ C#19 NOT OK MSOLVER: "No solution ?" idem if SOLVE for ni. SOLVE separately works for Vbi;xdmax;G0 but failed for ID. Algebraic calculation
@ also fails "Inconsistent units", cannot go further.
'ROOT(ⒺJFETs;[ni;Vbi;xdmax;G0;ID;VDsat;Vt;gm];[1_(cm^-3);1_V;1_μ;1_S;1_mA;1_V;1_V;1_mA/V])'
```


## Stress Analysis
The 28 variables in the Stress Analysis section are:

* `δ`: Elongation
* `ϵ`: Normal strain
* `γ`: Shear strain (dim.: angle)
* `φ`: Angle of twist
* `σ`: Normal stress (dim.: pressure, in SI: pascal, Pa)
* `σ1`: Maximum principal normal stress (dim.: pressure, in SI: pascal, Pa)
* `σ2`: Minimum principal normal stress (dim.: pressure, in SI: pascal, Pa)
* `σavg`: Normal stress on place of maximum shear stress (dim.: pressure, in SI: pascal, Pa)
* `σx`: Normal stress in x direction (dim.: pressure, in SI: pascal, Pa)
* `σx1`: Normal stress in rotated-x direction (dim.: pressure, in SI: pascal, Pa)
* `σy`: Normal stress in y direction (dim.: pressure, in SI: pascal, Pa)
* `σy1`: Normal stress in rotated-y direction (dim.: pressure, in SI: pascal, Pa)
* `τ`: Shear stress (dim.: pressure, in SI: pascal, Pa)
* `τmax`: Maximum shear stress (dim.: pressure, in SI: pascal, Pa)
* `τx1y1`: Rotated shear stress (dim.: pressure, in SI: pascal, Pa)
* `τxy`: Shear stress (dim.: pressure, in SI: pascal, Pa)
* `θ`: Rotation angle
* `θp1`: Angle to plane of maximum principal normal stress
* `θp2`: Angle to plane of minimum principal normal stress
* `θs`: Angle to plane of maximum shear stress
* `A`: Area
* `E`: Modulus of elasticity (dim.: pressure, in SI: pascal, Pa)
* `G`: Shear modulus of elasticity (dim.: pressure, in SI: pascal, Pa)
* `J`: Polar moment of inertia (dim.: length^4, in SI: m^4)
* `L`: Length
* `P`: Load (dim.: force, in SI: newton, N)
* `r`: Radius
* `T`: Torque (dim.: force·length, in SI: N·m)

#### Normal Stress

![Normal Stress](img/NormalStress.bmp)

* To calculate `[δ_in;ε;σ_psi]` (Elongation; Normal strain; Normal stress) from 4 known variables:
```rpl
P=40000_lbf  L=1_ft  A=3.14159265359_in^2  E=10e6_psi
@ Failing [ σ=12 732.39544 73 psi δ=0.00127 32395 45 ft 45 ε=0.00127 32395 45 ] 
@ C#20 NOT OK MSOLVER: "Inconsistent units". SOLVE individually works for σ but indicate for δ "Inconsistent units" for ε "Sign reversal" NOTE: HP50g_AUR has an error in the calculated value of δ
'ROOT(ⒺNormal Stress;[σ;δ;ε];[1_psi;1_in;1])'
```

#### Shear Stress

![Shear Stress](img/ShearStress.bmp)

* To calculate `[T_ft*lbf;γ_°;φ_°]` (Torque; Shear strain; Angle of twist) from 5 known variables:
```rpl
    "Shear Stress",  "{ "
L=6_ft  r=2_in  J=10.4003897419_in^4  G=12000000_psi  τ=12000_psi
@ Failing [ T=5 200.19487 095 ft·lbf γ=0.05729 57795 13 ° φ=2.06264 80624 7 ° ] 
@ C#21 NOT OK MSOLVER: "Inconsistent units". SOLVE individually works for T, but angle have to be in radians: corrections of eqns to be checked
'ROOT(ⒺShear Stress;[T;γ;φ];[1_ft*lbf;1_°;1_°])'
```

#### Stress On An Element

Stresses and strains are positive in the directions shown in the picture.

![Stress On An Element](img/StressOnAnElement.bmp)

* To calculate `[σx1_kPa;σy1_kPa;τx1y1_kPa]` (Normal stress in rotated-`x` direction; Normal stress in rotated-`y` direction; Rotated shear stress) from 4 known variables:
```rpl
    "Shear Stress",  "{ "
σx=15000_kPa  σy=4755_kPa  τxy=7500_kPa  θ=30_°
@ Failing [ σx1=18 933.94052 84 kPa σy=4 755 kPa τx1y1=-686.21513 0886 kPa ] 
@ C#22 NOT OK MSOLVER: "Inconsistent units". SOLVE individually works for σx1;σy1 but correction to eqn (3) to be checked
'ROOT(ⒺStress On An Element;[σx1;σy1;τx1y1];[1_kPa;1_kPa;1_kPa])'
```

#### Mohr’s Circle

![Mohr’s Circle](img/Mohr’sCircle.bmp)

* To calculate `[σ1_psi;σ2_psi;θp1_°;θp2_°;τmax_psi;θs_°;σavg_psi]` (Maximum principal normal stress; Minimum principal normal stress; Angle to plane of maximum principal normal stress; Angle to plane of minimum principal normal stress; Maximum shear stress; Angle to plane of maximum shear stress; Normal stress on place of maximum shear stress) from 3 known variables:
```rpl
σx=-5600_psi  σy=-18400_psi  τxy=4800_psi
@ Expecting [ σ1=-1 755.00122 011 psi σ2=-22 244.99877 99 psi θp1=18.43494 88229 ° θp2=108.43494 8823 ° τmax=10 244.99877 99 psi θs=-26.56505 11771 ° σavg=-12 000. psi ][ σ1=-1 755.00122 011 psi σ2=-22 244.99877 99 psi θp1=18.43494 88229 ° θp2=108.43494 8823 ° τmax=10 244.99877 99 psi θs=-26.56505 11771 ° σavg=-12 000. psi ] 
'ROOT(ⒺMohr’s Circle;[σ1;σ2;θp1;θp2;τmax;θs;σavg];[1_psi;1_psi;1_°;1_°;1_psi;1_°;1_psi])'
```


## Waves
The 36 variables in the Waves section are:

* `β`: Sound intensity level (dim.: intensity in logarithmic scale, in SI: decibel, dB)
* `Δp`: Sound pressure variafion around atmospheric pressure (dim.: force/area, in SI: pascal, Pa)
* `Δpm`: Amplitude of sound pressure variafion around atmospheric pressure (dim.: force/area, in SI: pascal, Pa)
* `λ`: Wavelength
* `ω`: Angular frequency (dim.: angle/time. in SI: r/s^2)
* `ρ`: Volumic density of medium (dim.: mass/volume, in SI: kg/m^3)
* `θcone`: Angle at the tip of the cone formed by the supersonic shockwave
* `μ`: Linear density of a string (dim.: mass/length, in SI: kg/m)
* `as`: Acceleration at `x` and `t` of vibrating particles ([Longitudinal Waves](#Longitudinal Waves)), or Acceleration at x and t of air particles ([Sound Waves](#Sound Waves))
* `ay`: Acceleration at `x` and `t` of vibrating particles for transversal waves
* `B`: Bulk modulus of elasticity ([Sound Waves](#Sound Waves)) (dim.: pressure, in SI: Pa), or Magnetic field ([Electromagnetic Waves](#Electromagnetic Waves)) (dim.: mass/(time^2·current), in SI: T)
* `E`: Electric field at `x` and `t` of an electromagnetic wave (dim.: force/charge, in SI: N/C=V/m)
* `Em`: Electric Field amplitude (dim.: force/charge, in SI: N/C=V/m)
* `f`: Frequency (dim.: 1/time; in SI: hertz, Hz)
* `favg`: Frequency average (dim.: 1/time; in SI: hertz, Hz)
* `fbeat`: Beat frequency (dim.: 1/time; in SI: hertz, Hz)
* `f0`: Emission frequency (dim.: 1/time; in SI: hertz, Hz)
* `ffixed-fixed`: frequency of harmonics on a string fixed at both ends (dim.: 1/time; in SI: hertz, Hz)
* `ffixed-free`: frequency of harmonics on a string fixed at one end and free at the other end (dim.: 1/time; in SI: hertz, Hz)
* `fopen-open`: frequency of harmonics in a tube open at both ends (dim.: 1/time; in SI: hertz, Hz)
* `fopen-close`: frequency of harmonics in a tube open at one end and close at the other end (dim.: 1/time; in SI: hertz, Hz)
* `I`: Sound intensity (dim.: energy/(area·time), in SI: W/m^2)
* `k`: Angular wave number (dim.: agle/length, in SI: r/m)
* `ninteger`: Harmonics number being an integer number
* `nodd`: Harmonics number being an odd number
* `Ps`: Power of the source (dim.: energy/time, in SI: watt, W)
* `s`: Longitudinal displacement at `x` and `t` of vibrating particles ([Longitudinal Waves](#Longitudinal Waves)), or Longitudinal displacement at `x` and `t` of air particles ([Sound Waves](#Sound Waves))
* `sm`: Longitudinal displacement amplitude of vibrating particles ([Longitudinal Waves](#Longitudinal Waves)), or Longitudinal displacement amplitude of air particles ([Sound Waves](#Sound Waves))
* `t`: Time
* `u`: Mass or flow velocity ([Mach Number](#Mach Number))
* `v`: Velocity of the propagating sound in medium ([Sound Waves](#Sound Waves)), or Wave speed ([Transverse Waves](#Transverse Waves)) & ([Longitudinal Waves](#Longitudinal Waves))
* `vs`: Velocity at `x` and `t` of vibrating particles ([Longitudinal Waves](#Longitudinal Waves)), or Velocity at `x` and `t` of air particles ([Sound Waves](#Sound Waves))
* `vsair`: Velocity of the propagating sound in the ait as a function of temperature
* `vy`: Velocity at `x` and `t` of vibrating particles for transversal waves
* `x`: Position
* `y`: Transverse displacement at `x` and `t` 

#### Transverse Waves

* To calculate `[f_Hz;λ_cm;v_cm/s;y_cm;vy_cm/s;ay_cm/s^2]` (Frequency; Wavelength; Transverse displacement at x and t; Wave speed; Velocity & Acceleration at `x` and `t` of vibrating particles) from 5 known variables:
```rpl
ym=6.37_cm  k=32.11_r/cm  x=0.03_cm  ω=7000_r/s  t=1_s
@ Expecting [ f=1 114.08460 164 Hz λ=0.19567 69015 cm v=218.00062 2859 cm/s y=2.87199 87256 6 cm vy=-39 800.72414 76 cm/s ay=312 130 000. cm/s↑2 ] 
'ROOT(ⒺTransverse Waves;[f;λ;v;y;vy;ay];[1_Hz;1_cm;1_cm/s;1_cm;1_cm/s;1_cm/s^2])'
```

#### Longitudinal Waves

* To calculate `[s_cm;λ_cm;f_Hz;v_m/s;vs_cm/s;as_cm/s^2]` (Frequency; Wavelength; Transverse displacement at x and t; Wave speed; Velocity & Acceleration at `x` and `t` of vibrating particles) from 5 known variables:
```rpl
sm=6.37_cm  k=32.11_r/cm  x=0.03_cm  ω=7000_r/s  t=1_s
@ Expecting [ s=5.68581 77353 7 cm λ=0.19567 69015 cm f=1 114.08460 164 Hz v=2.18000 62285 9 m/s vs=20 103.99107 96 cm/s as=278 605 069.033 cm/s↑2 ] 
'ROOT(ⒺLongitudinal Waves;[s;λ;f;v;vs;as];[1_cm;1_cm;1_Hz;1_m/s;1_cm/s;1_cm/s^2])'
```

#### Sound Waves

* To calculate `[v_m/s;f_Hz;I_W/m^2;β_dB;s_cm;vs_cm/s;as_cm/s^2;Δpm_Pa;Δp_Pa;Ps_W]` (Wave speed; Frequency; Sound intensity; Sound intensity level in dB; Longitudinal displacement, Velocity & Acceleration at `x` and `t` of vibrating particles; Amplitude of sound pressure variafion around atmospheric pressure; Sound pressure variafion; Power of the source) from 9 known variables:
```rpl
sm=10_cm  ω=6000_r/s  B=12500_kPa  ρ=65_kg/m^3   x=2_cm   t=0_s  r=10_m   k=13.6821_r/m  φ=2_r
@ Expecting [ v=438.52900 9654 m/s f=954.92965 8551 Hz I=5.13078 94129 5⁳⁹ W/m↑2 s=7.63005 85995 9 cm vs=-38 783.49400 98 cm/s as=274 682 109.585 cm/s↑2 Δpm=17 102 631.3765 Pa Δp=11 054 996.6924 Pa Ps=6.44754 01307 3⁳¹² W ]
@ Failing [ v=438.52900 9654 m/s f=954.92965 8551 Hz I=5.13078 94129 5⁳⁹ W/m↑2 s=7.63005 85995 9 cm vs=-38 783.49400 98 cm/s as=-274 682 109.585 cm/s↑2 Δpm=17 102 631.3765 Pa Δp=11 054 996.6924 Pa Ps=6.44754 01307 3⁳¹² W ]
@ C#23 NOT OK : the sign of the acceleration is in error, here as < 0 see ISSUE # 1353
'ROOT(ⒺSound Waves;[v;f;I;s;vs;as;Δpm;Δp;Ps];[1_m/s;1_Hz;1_(W/(m^2));1_cm;1_cm/s;1_cm/s^2;1_Pa;1_Pa;1_W])'
```

#### Doppler Effect

In the classical Doppler effect it is assumed that the speed of the observer and the source are lower than the speed of sound in the air. The speed of the receiver relative to the air `vr` is added to the speed of sound in the air if the receiver is moving towards the source, subtracted if the receiver is moving away from the source whatever the movement of the source. The speed of the source relative to the air `va` is subtracted from the speed of sound in the air if the source is moving towards the receiver, added if the source is moving away from the receiver whatever the movement of the receiver.

#### Mach Number

For an object moving at a supersonic speed, the shockwave describes a cone having the angle `θcone` at its tip where the opposite side is the distance travelled by the sound and the hypothenuse is the distance travelled by the object.

#### String Standing Waves

A string being fixed or free at its ends admits only discrete harmonics as standing waves on the string. A string being fixed (or free) at both ends admits all integer harmonics. A string being being fixed at one end and free at the other end admits only all odd integer harmonics.

#### Sound Wave Harmonics

A tube being open or closed at its ends admits only discrete harmonics as standing waves of the sound in the air within the tube. A tube being open (or closed) at both ends admits all integer harmonics. A tube being being open at one end and closed at the other end admits only all odd integer harmonics.

#### Beat Acoustics

In acoustics, a beat is an interference pattern between two sounds of slightly different frequencies, perceived as a periodic variation in amplitude whose rate is the difference of the two frequencies. The sum of two unit-amplitude sine waves can be expressed as a carrier wave of frequency `favg` whose amplitude is modulated by an envelope wave of frequency `fbeat`.

#### Electromagnetic Waves


## Relativity
The 107 variables in the Relativity section are:

* `α`: Light Doppler effect, light arrival angle in the rest frame
* `β`: Velocity relativistic speed ratio
* `βe`: Escape velocity relativistic speed ratio
* `βp`: Plane velocity relativistic speed ratio
* `βg`: Ground velocity relativistic speed ratio
* `Δt`: Proper time interval (Time Dilation), or Duration of the circumnavigation trip at latitude `φ` for a non-rotating planet (Circumnavigating Airplanes)
* `Δtp`: Dilated time interval
* `ΔtpG`: Gravitationally dilated time interval
* `Δx`: Proper space interval
* `Δxp`: Contracted space interval
* `λ1, λ2`: Wavelength of the photon as measured by the observer at positions `R1` and `R2`
* `λ∞`: Wavelength of the photon as measured by the observer at infinity
* `γ`: Lorentz factor
* `γ21`: Factor of combined special and general relativity effects
* `γv1, γv2`: Lorentz factor for velocities `v1` and `v2`
* `γG`: Lorentz factor associated to gravitational dilation
* `γG1, γG2`: Lorentz factor associated to gravitational dilation at heights `h1` and `h2`
* `ω`: Angular velocity of a rotating planet during a day (dim.: angle/time, in SI: r/s)
* `φ`: Latitude (dim.: angle)
* `θ` : Aberration of light, emission angle in the frame at rest
* `θp` : Aberration of light, emission angle in the moving frame
* `ΔτWE`: Time difference between westward and eastward flights
* `ΔτE`: Flight time in the eastward direction
* `ΔτW`: Flight time in the westward direction
* `Δτg`: Elapsed time variation due to the ground tangential velocity
* `ΔτpE`: Elapsed time variation due to the plane altitude and velocity in the eastward direction
* `ΔτpW`: Elapsed time variation due to the plane altitude and velocity in the westward direction
* `As`: Schwarzschild  black hole surface area
* `ax`: X component of the acceleration
* `ay`: Y component of the acceleration
* `az`: Z component of the acceleration
* `apx`: Transformed x component of the acceleration observed in the moving frame
* `apy`: Transformed y component of the acceleration observed in the moving frame
* `apz`: Transformed z component of the acceleration observed in the moving frame
* `Bx`: X component of the magnetic field (dim.: mass/(time^2·current), in SI: tesla, T)
* `By`: Y component of the magnetic field (dim.: mass/(time^2·current), in SI: tesla, T)
* `Bz`: Z component of the magnetic field (dim.: mass/(time^2·current), in SI: tesla, T)
* `Bpx`: Transformed x component of the magnetic field (dim.: mass/(time^2·current), in SI: tesla, T)
* `Bpy`: Transformed y component of the magnetic field (dim.: mass/(time^2·current), in SI: tesla, T)
* `Bpz`: Transformed z component of the magnetic field (dim.: mass/(time^2·current), in SI: tesla, T)
* `E`: Total energy
* `Ep`: Transformed total energy
* `E0`: Rest energy associated to the rest mass
* `Ex`: X component of the electric field (dim.: force/charge, in SI: N/C=V/m)
* `Ey`: Y component of the electric field (dim.: force/charge, in SI: N/C=V/m)
* `EZ`: Z component of the electric field (dim.: force/charge, in SI: N/C=V/m)
* `Epx`: Transformed x component of the electric field (dim.: force/charge, in SI: N/C=V/m)
* `Epy`: Transformed y component of the electric field (dim.: force/charge, in SI: N/C=V/m)
* `Epz`: Transformed z component of the electric field (dim.: force/charge, in SI: N/C=V/m)
* `f` = Light Doppler effect, frequency received in the frame at rest (dim.: 1/time, in SI: hertz, Hz)
* `fp` = Light Doppler effect, frequency emitted in the moving frame (dim.: 1/time, in SI: hertz, Hz)
* `fs`: Wave frequency of the source (dim.: 1/time, in SI: hertz, Hz)
* `frl`: Relativistic longitudinal Doppler frequency (dim.: 1/time, in SI: hertz, Hz)
* `frt`: Relativistic transverse Doppler frequency (dim.: 1/time, in SI: hertz, Hz)
* `h1, h2`: Height of points 1 and 2
* `hp`: Altitude of plane (dim.: length)
* `K`: Kinetic energy
* `M`: Mass of planet
* `m0`: Rest mass
* `MGu`: Reduced gravitational mass given in geometrized units (dim.: length)
* `Mxearth`: Mass factor as a multiple of Earth mass
* `Mxsun`: Mass factor as a multiple of Sun mass
* `MxSagA`: Mass factor as a multiple of SagitariusA* mass
* `Pθ` = Light relativistic beaming, angular distribution of photon in the moving frame from a source isotropic and stationary
* `p`: Momentum (dim.: mass·speed, in SI: kg·m/s)
* `px`: X component of the momentum (dim.: mass·speed, in SI: kg·m/s)
* `py`: Y component of the momentum (dim.: mass·speed, in SI: kg·m/s)
* `pz`: Z component of the momentum (dim.: mass·speed, in SI: kg·m/s)
* `ppx`: Transformed x component of the momentum (dim.: mass·speed, in SI: kg·m/s)
* `ppy`: Transformed y component of the momentum (dim.: mass·speed, in SI: kg·m/s)
* `ppz`: Transformed z component of the momentum (dim.: mass·speed, in SI: kg·m/s)
* `PBH`: Black hole evaporation power (dim.: energy/time, in SI: watt, W)
* `R`: Planet radius
* `R1, R2`: Radius to points 1 and 2
* `Rem`: Radius at which the photon is emitted
* `rs`: Schwarzschild  radius
* `rxearth`: Radius factor as a multiple of Earth radius
* `SBH`: Black hole entropy (dim.: energy/temperature, in SI: J/K)
* `Tday`: Day duration
* `tev`: Evaporation time of a black hole
* `TH`: Black hole temperature
* `txyr`: Duration factor as a multiple of a year
* `ux`: X component of the velocity
* `uy`: Y component of the velocity
* `uz`: Z component of the velocity
* `upx`: Transformed x component of the velocity observed in the moving frame
* `upy`: Transformed y component of the velocity observed in the moving frame
* `upz`: Transformed z component of the velocity observed in the moving frame
* `v1, v2`: Velocity at heights `h1` and `h2`
* `vg`: Tangential speed of the ground of a rotatong planet
* `Vs`: Schwarzschild  volume of a black hole
* `Vxsun`: Volume factor as a multiple of Sun volume
* `x`: X coordinate in the rest frame
* `y`: Y coordinate in the rest frame
* `z`: Z coordinate in the rest frame
* `xp`: Transformed X coordinate in the moving frame
* `yp`: Transformed Y coordinate in the moving frame
* `zp`: Transformed Z coordinate in the moving frame
* `v`: Velocity along the x axis
* `ve`: Escape velocity
* `z`: Gravitational redshift parameter

The relativistic transformations are parametrized by the real constant `v` representing a velocity confined to the x-direction. The respective inverse transformation is then parameterized by the negative of this velocity.

#### Lorentz Transformation

The primed reference frame `[xp yp zp]` is travelling with velocity `v` in the positive x direction. Therefore, the y and z coordinates of the rest frame remain unchanged.

#### Time Dilation

The dilation comes from the fact that the Lorentz factor `γ` is greater or equal to one and the proper time interval is multiplied by this factor.

#### Space Contraction

The contraction comes from the fact that the Lorentz factor `γ` is greater or equal to one and the proper space interval is divided by this factor.

#### Velocity Superposition

These expressions replace the usual Galilean addition of velocities. It can be checked that superposing with `v = c` leads to `upx = c`, hence the impossibility to superpose velocities to go beyond the velocity limit `c`. Since the velocity `v` is confined to the x-direction, the y and z components of velocity remain unchanged.

#### Acceleration Superposition

Even if the velocity `v` is confined to the x-direction, all components of the observed acceleration are transformed in the moving frame.

#### E & B Fields Transformation

#### Longitudinal Doppler Effect

#### Transverse Doppler Effect

#### Light Propagation

#### Energy & Momentum

The total relativistic energy `E` and the norm of the momentum `p` form the invariant `mo·c^2` which remains the same in all frames. The kinetic energy `K` is the difference between the total relativistic energy `E` and the rest energy `E0 = mo·c^2`.

#### Ultrarelativistic Cases

#### Gravitational Time Dilation

#### Gravitational Redshift

#### Circumnavigating Airplanes

It is assumed that the planes are circumnavigating at the same altitude `h`, same latitude `φ` and the during same flight duration `Δt` measured in the airplanes. The ground is rotating with the planet at the angular frequency `ω`. The Schwarzschild metric is taken into account. The calculation should formally invokes an integral for the elapsed proper time along a path and is approximated here to the first order in speed since the velocities of the plane and the planet surface are << `c` (slow airplanes over a slow rotating planet).

#### Clocks at different heights

It is assumed that the two clocks are at rest with respect to the ground at a latitude `φ` and are rotating with the planet at the angular frequency `ω`. The clocks are at their respective heights `h1` and `h2` for inxtance at the top and bottom of a mountain. For simplicity, the planet is assumed to have a spherical distribution. The Schwarzschild metric is taken into account. The calculation should formally invokes an integral for the elapsed proper time along a path and is approximated here to the first order in speed since the tangential velociies at height `h1` and `h2` are << `c` (slow rotating planet).

#### B H Schwarzschild Geometry

#### B H Thermodynamics

## Modern Physics
The 37 variables in the Modern Physics section are:

* `β`: Velocity relativistic speed ratio
* `φ`: Work function of the substance (dim.: energy, in SI: eV)
* `θ`: Diffused photon angle of deflection (Photoelectric Effect, Compton Effect), or Angle between incident photon and cristallographic plane (DeBroglie Wave)
* `γ`: Lorentz relativistic factor
* `%rFr12`: Relative % of change between distribution fractions integrated from f1 to f2
* `%rFr34`: Relative % of change between distribution fractions integrated from f3 to f4
* `λ`: Wavelength of the incident photon (Compton Effect)
* `λp`: Wavelength of the diffused photon
* `A`: Area, or Total activity (Radioactivity) (dim.: number of decay/time, in SI: becquerel, Bq)
* `d`: Interatomic distance between cristallographic plane
* `En`: Electron energy in level n (dim.: energy, in SI: eV)
* `Enp`: Electron energy in level np (dim.: energy, in SI: eV)
* `ebfafb`: Emissive power in the range fa to fb (dim.: energy/time·area, in SI: W/m^2)
* `eb`: Total emissive power for the entire Planck spectrum (dim.: energy/time·area, in SI: W/m^2)
* `Eph`: Photon energy (dim.: energy; in SI: eV)
* `f`: Frequency of the photon (dim.: 1/time, in SI: hertz, Hz)
* `f0`: Threshold frequency (dim.: 1/time, in SI: hertz, Hz)
* `f1, f2`: Lower and upper frequency limits of integration (dim.: 1/time, in SI: hertz, Hz)
* `f3, f4`: Lower and upper frequency limits of integration (dim.: 1/time, in SI: hertz, Hz)
* `fpeak`: Frequency of maximum emissive power for the Planck distribution (dim.: 1/time, in SI: hertz, Hz)
* `Frfafb`: Fraction of total emissive power in the range fa to fb
* `FrPl12`: Fraction of Planck emissive power in the range f1 to f2
* `FrWn12`: Fraction of Wien emissive power in the range f1 to f2
* `FrRJ12`: Fraction of Rayleign-Jeans emissive power in the range f1 to f2
* `FrPl34`: Fraction of Planck emissive power in the range f3 to f4
* `FrWn34`: Fraction of Wien emissive power in the range f3 to f4
* `FrRJ34`: Fraction of Rayleign-Jeans emissive power in the range f3 to f4
* `K`: Kinetic energy of ejected electron (dim.: energy; in SI: eV)
* `Kmax`: Maximum kinetic energy of photoelectron (dim.: energy; in SI: eV)
* `p`: Momentum (dim.: mas·speed, in SI: kg·m/s)
* `m`: Mass
* `q`: Heat transfer rate (dim.: energy/time, in SI: watt, W)
* `r`: Radius of electron orbit in level n
* `T`: Temperature
* `V0`: Stopping potential (dim.: energy/charge, In SI: volt,V)

#### Planck & Wien Comparison

In this section, two comparisons are done between the Planck and Wien spectral distributiona. Based on a incomplete thermodynamic argument, the latter is an approximation of the true Planck law describing the spectral distribution for the light emitted by a black-body. The choice of temperature `T` determines the frequency ranges for integration between `f1` and `f2`, or between `f3` and `f4`. One shall determine in which frequency interval both distribution differs notably or agree. The asymptotic agreement for large frequency is clearly illustrated in the picture. The user is free to choose one or the other comparison fractions (replacing it in `Frfafb`) to compute the corresponding enissive power and the heat transfer rate from the black-body.

![Planck & Wien Comparison](img/Planck&Wien_Distributions.bmp)

#### Planck & Rayleigh-Jeans Comparison

In this section, two comparisons are done between the Planck and Rayleigh-Jeans spectral distributiona. Based on the equipartition theorem argument, the latter is an approximation of the true Planck law describing the spectral distribution for the light emitted by a black-body. The choice of temperature `T` determines the frequency ranges for integration between `f1` and `f2`, or between `f3` and `f4`. One shall determine in which frequency interval both distribution agree or differs considerably, leading to a divergence called UV catastrophy corresponding to unphysical fractions greather than one. The asymptotic agreement for small frequency is clearly illustrated in the picture. The user is free to choose one or the other comparison fractions (replacing it in `Frfafb`) to compute the corresponding enissive power and the heat transfer rate from the black-body.

![Planck & Rayleigh-Jeans Comparison](img/Planck&Rayleigh-Jeans_Distributions.bmp)

#### Photoelectric Effect

Einstein explained the photoelectric effect with the energy quantification of the electromagnetic wave. The photoelectron is then emitted only if the energy `E` of the incident photon is greather or equal to the work function `φ` of the material. A simple energy budget equation determines the maximum kinetic energy `Kmax` of the photoelectron.

#### Compton Effect

In the Compton effect, both energy and momentum are conserved during the collision of the incident photon and the electron, which underlines the fact that the photon must henceforth be considered as a particle. When a high frequency `f` (or energy `E = hf`) photon scatters due to an interaction with a charged particle, there is a decrease in the energy of the photon emitted at an angle `θ` and thus, an increase in its wavelength `λp`. The energy of the ejected electron `Kmax` is relativist.

#### De Broglie Wave

At all scales where measurements have been possible, matter exhibits wave-like behavior (e.g. Young interference for protons, even for molecules). More precisely, a beam of neutron can be diffracted just like a beam of light or a water wave as it is the case in the Bragg diffraction. Here, the interference is constructive when the phase difference between the matter wave reflected off different atomic planes at an angle `θ` is a multiple of 2π, giving the Following condition: `2·d·SINθ = n·λ` with `n` integer and where `θ` is mesured between the refected waves trajectory and the atomic plane.

#### Bohr Atomic Model

## Nuclear Physics
The 39 variables in the Nuclear Physics section are:

* `λ`: Decay constant (Radioactivity) (dim.: 1/time, in SI: s^-1)
* `A`: Total activity (Radioactivity) (dim.: number of decay/time, in SI: becquerel, Bq), or Mass number (Nuclear Physics)
* `AXα`: Mass number of the radionuclide X undergoing α decay
* `AYα`: Mass number of the daughter nuclide Y from α decay
* `AXβ⊖`: Mass number of the radionuclide X undergoing β- decay
* `AYβ⊖`: Mass number of the daughter nuclide Y from β- decay
* `AXβ⊕`: Mass number of the radionuclide X undergoing β+ decay
* `AYβ⊕`: Mass number of the daughter nuclide Y from β+ decay
* `Aa`: Mass number of the incident nuclide or particle a
* `Ab`: Mass number of the product nuclide or particle b
* `AX`: Mass number of the reactant nuclide X
* `AY`: Mass number of the product nuclide Y
* `A0`: Initial total activity (dim.: number of decay/time, in SI: becquerel, Bq)
* `EB`: Nuclear binding energy (dim.: energy, in SI: MeV)
* `EBse`: Semiempirical nuclear binding energy (dim.: energy, in SI: MeV)
* `ma`: Mass of incident nuclide or particle a (In SI: u)
* `mb`: Mass of outgoing nuclide or particle b (In SI: u)
* `mX`: Mass of reactant neutral atom X (In SI: u)
* `mY`: Mass of product neutral atom X (In SI: u)
* `MW`: Mass weight (dim.: mass/mol, in SI: g/mol)
* `N`: Number of nucleid (Radioactivity), or Number of neutron (Nuclear Physics)
* `N0`: Initial number of nucleid
* `Q`: Net energy balance of a nuclear reaction (dim.: energy, in SI: MeV)
* `Qα`: Net energy balance of an α decay (dim.: energy, in SI: MeV)
* `Qβ⊖`: Net energy balance of a β- decay (dim.: energy, in SI: MeV)
* `Qβ⊕`: Net energy balance of a β+ decay (dim.: energy, in SI: MeV)
* `R`: Radius of the nucleus having `A` nucleons
* `T1/2`: Half-life of radionuclide (dim.: time)
* `Z`: Number of proton
* `ZXα`: Proton number of the radionuclide X undergoing α decay
* `ZYα`: Proton number of the daughter nuclide Y from α decay
* `ZXβ⊖`: Proton number of the radionuclide X undergoing β- decay
* `ZYβ⊖`: Proton number of the daughter nuclide Y from β- decay
* `ZXβ⊕`: Proton number of the radionuclide X undergoing β+ decay
* `ZYβ⊕`: Proton number of the daughter nuclide Y from β+ decay
* `Za`: Proton number of the incident nuclide or charge of the incident particle a
* `Zb`: Proton number of the incident nuclide or or charge of the product particle b
* `ZX`: Proton number of the reactant nuclide X
* `ZY`: Proton number of the product nuclide Y

For all nuclear reactions, including nuclear decays, we have charge conservation `Zp = Z` and mass number conservation `Mp = M` (the same number of nucleons). This therefore requires assigning numbers A and Z to incident particles `a` and `b` whether they are for example gamma photons or neutrinos (`A = 0`, Z = 0`), positrons (`A = 0`, Z = +1`) or others. The reaction energy `Q` is always calculated with the mass-energy equivalence by the mass difference between the reactants and the products. For instance, spontaneous decays and fissions are always exothermic `Q > 0` while some nuclear reactions like can also be endothermic `Q < 0` (like nucleosynthesis of heavy elements). The mass difference should be computed with at least 8 significant digits because reactants and products masses can be of the same order of magnitude. 

#### Radioactivity

#### Radius & Binding Energy

#### α Decay

#### β⊖ Decay

#### β⊕ Decay

#### General Nuclear Reaction
