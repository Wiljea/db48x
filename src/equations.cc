// ****************************************************************************
//  equations.cc                                                  DB48X project
// ****************************************************************************
//
//   File Description:
//
//      Implementation of the equations library
//
//
//
//
//
//
//
//
// ****************************************************************************
//   (C) 2024 Christophe de Dinechin <christophe@dinechin.org>
//   This software is licensed under the terms outlined in LICENSE.txt
// ****************************************************************************
//   This file is part of DB48X.
//
//   DB48X is free software: you can redistribute it and/or modify
//   it under the terms outlined in the LICENSE.txt file
//
//   DB48X is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// ****************************************************************************

#include "equations.h"

#include "expression.h"
#include "grob.h"
#include "renderer.h"
#include "solve.h"
#include "variables.h"

RECORDER(equations,         16, "Equation objects");
RECORDER(equations_error,   16, "Error on equation objects");


// ============================================================================
//
//   Equation definitions
//
// ============================================================================

static const cstring basic_equations[] =
// ----------------------------------------------------------------------------
//   List of basic equations
// ----------------------------------------------------------------------------
//   clang-format off
{
    // ------------------------------------------------------------------------
    //   Columns and beams
    // ------------------------------------------------------------------------

    "Columns and Beams", nullptr,

    "Elastic Buckling",  "{ "
    "  '(Pcr_kN)=(Ⓒπ²*(E_kPa)*(A_cm^2))/sq((K*(L_m))/(r_cm))' "
    "  '(Pcr_kN)=(Ⓒπ²*(E_kPa)*(I_mm^4))/sq(K*(L_m))' "
    "  '(σcr_kPa)=(Pcr_kN)/(A_cm^2)' "
    "  '(r_cm)= sqrt((I_mm^4)/(A_cm^2))' "
    "}",

    "Eccentric Columns", "{"
    "  '(σmax_kPa)=((P_kN)/(A_cm^2))*(1+((ε_cm)*(c_cm))/sq(r_cm)*inv(cos(K/2*((L_m)/(r_cm))*sqrt((P_kN)/((E_kPa)*(A_cm^2)))*1_r)))'"
    "  '(r_cm)=sqrt((I_mm^4)/(A_cm^2))'"
    "}",

    "Simple Deflection", "{"
    "  '(y_cm)=((P_kN)*((L_m)-(a_m)))*(x_m)/(6*(L_m)*(E_kPa)*(I_mm^4))*(x²+(L-a)²-L²)-((M_N*m)*x)/(E*I)*((c_m)-x²/(6*L)-L/3-c²/(2*L))-((w_N/m)*x)/(24*E*I)*(L³+x²*(x-2*L))'"
    "}",

    "Simple Slope", "{"
    "  '(θ_°)=(1_r)*(((P_N)*((L_m)-(a_m)))/(6*L*(E_kPa)*(I_mm^4))*(3*(x_m)²+(L-a)²-L²)-(M_N*m)/(E*I)*(c-(x²-c²)/(2*L)-L/3)-(w_N/m)/(24*E*I)*(L³+x²*(4*x-6*L)))'"
    "}",

    "Simple Moment", "{"
    "  '(Mx_N*m)=(P_N)*IFTE((x_m)≤(a_m);((L_m)-(a_m))*(x_m);(L-x)*a)/L+(M_N*m)*IFTE(x≤(c_m);x;x-L)/L+((w_N/m)*x*(L-x))/2'"
    "}",

    "Simple Shear", "{"
    "  '(V_N)=((P_N)*((L_m)-(a_m)))/L+(M_N*m)/L+((w_N/m)*(L-2*(x_m)))/2'"
    "}",

    "Cantilever Deflection", "{"
    "  '(y_m)=((P_N)*IFTE((x_m)≤(a_m);x;a)²)/(6*(E_kPa)*(I_mm^4))*IFTE(x≤a;x-3*a;a-3*x)+(M_N*m)*IFTE(x≤(c_m);x²;c*(2*x-c))/(2*E*I)-((w_N/m)*x²)/(24*E*I)*(6*L²-4*L*x+x²)'"
    "}",

    "Cantilever Slope", "{"
    "  '(θ_°)=(1_r)*((P_N)*IFTE((x_m)≤(a_m);x*(x-2*a);-a²)/(2*(E_kPa)*(I_mm^4))+(M_N*m)*IFTE(x≤(c_m);x;c)/(E*I)-((w_N/m)*x)/(6*E*I)*(3*L²-3*L*x+x²))'"
    "}",

    "Cantilever Moment", "{"
    "  '(Mx_N*m)=IFTE((x_m)≤(a_m);1;0)*(P_N)*((x_m)-(a_m))+IFTE(x≤(c_m);1;0)*(M_N*m)-(w_N/m)/2*(L²-2*L*x+x²)'"
    "}",


    "Cantilever Shear", "{"
    "  '(V_n)=IFTE((x_m)≤(a_m);1;0)*(P_N)+(w_N/m)*((L_m)-(x_m))'"
    "}",

    "Electricity", nullptr,
// 67 eqns
    "Coulomb's Law & E Field",  "{ "
    "  '(F_N)=1/(4*Ⓒπ*Ⓒε0*εr)*((q1_C)*(q2_C)/(r_m)^2)' "
    "  '(E_(N/C))=(F_N)/(qtest_C)' "
    "}",

    "E Field Infinite Line",  "{ "
    "  '(E_(N/C))=1/(2*Ⓒπ*Ⓒε0*εr)*((λ_(C/m))/(r_m)' "
	"  '(λ_(C/m))=(Q_C)/(L_m)' "
    "}",

    "E Field Finite Line",  "{ "
    "  '(E_(N/C))=1/(4*Ⓒπ*Ⓒε0*εr)*((λ_(C/m))/(r_m)*((SIN(θ1_r)-SIN(θ2_r))' "
	"  '(λ_(C/m))=(Q_C)/(L_m)' "
    "}",

    "E Field Infinite Plate",  "{ "
    "  '(E_(N/C))=(σ_(μC/cm^2))/(2*Ⓒε0*εr)' "
	"  '(σ_(μC/m^2))=(Q_μC)/(A_(cm^2))' "
    "}",

    "Ohm's Law & Power",  "{ "
    "  '(V_V)=(I_A)*(R_Ω)' "
	"  '(P_W)=(V_V)*(I_A)' "
	"  '(P_W)=(I_A)^2*(R_Ω)' "
	"  '(P_W)=(V_V)^2/(R_Ω)' "
    "}",

    "Volt Divider",  "{ "
    "  '(V1_V)=(V_V)*((R1_Ω)/((R1_Ω)+(R2_Ω)))' "
    "}",

    "Current Divider",  "{ "
    "  '(I1_A)=(I_A)*((R2_Ω)/((R1_Ω)+(R2_Ω)))' "
    "}",

    "Wire Resistance",  "{ "
    "  '(R_Ω)=(ρ_(Ω*m))*(L_m)/(A_(m^2))' "
    "}",

    "Resistivity & Conductivity",  "{ "
    "  '(ρ_(Ω*m))=(ρ0_(Ω*m))*(1+(αT_K^-1)*((T_K)-(T0_K)))' "
	"  '(σ_(S/m))=1/(ρ_(Ω*m))' "
    "}",

    "Series & Parallel R",  "{ "
    "  '(Rs_Ω)=(R1_Ω)+(R2_Ω)'"
	"  '1/(Rp_Ω)=1/(R1_Ω)+1/(R2_Ω)' "
    "}",

    "Series & Parallel C",  "{ "
    "  '1/(Cs_μF)=1/(C1_μF)+1/(C2_μF)' "
	"  '(Cp_μF)=(C1_μF)+(C2_μF)' "
    "}",

    "Series & Parallel L",  "{ "
    "  '(Ls_mH)=(L1_mH)+(L2_mH)' "
	"  '1/(Lp_mH)=1/(L1_mH)+1/(L2_mH)' "
    "}",

    "Capacitive Energy",  "{ "
    "  '(E_J)=(1/2)*(C_μF)*(V_V)^2' "
	"  '(E_J)=(1/2)*(q_μC)*(V_V)' "
	"  '(E_J)=(q_μC)^2/(2*(C_μF))' "
    "}",

    "Volumic Density Electric Energy",  "{ "
    "  '(uE_(J/m^3))=(1/2)*Ⓒε0*εr*(E_(V/m))^2' "
	"}",

    "Inductive Energy",  "{ "
    "  '(E_J)=(1/2)*(L_mH)*(I_A)^2' "
    "}",

    "RLC Current Delay",  "{ "
    "  'TAN(φs_°)=((XL_Ω)-(XC_Ω))/(R_Ω)' "
	"  'TAN(φp_°)=(1/(XC_Ω)-1/(XL_Ω))*(R_Ω)' "
	"  '(XC_Ω)=1/((ω_(r/s))*(C_μF))' "
	"  '(XL_Ω)=(ω_(r/s))*(L_mH)' "
	"  '(ω_(r/s))=2*(Ⓒπ_r)*(f_Hz)' "
    "}",

    "DC Capacitor Current",  "{ "
    "  '(I_A)=(C_μF)*((ΔV_V)/(Δt_s))' "
	"  '(ΔV_V)=(Vf_V)-(Vi_V)' "
	"  '(Δt_μs)=(tf_μs)-(ti_μs)' "
    "}",

    "Capacitor Charge",  "{ "
    "  '(q_C)=(C_μF)*(V_V)' "
    "}",

    "DC Inductor Voltage",  "{ "
    "  '(V_V)=-(L_mH)*((ΔI_A)/(Δt_μs))' "
	"  '(ΔI_A)=(If_A)-(Ii_A)' "
	"  '(Δt_μs)=(tf_μs)-(ti_μs)' "
    "}",

    "RC Transient",  "{ "
    "  '(V_V)=(Vf_V)-((Vf_V)-(Vi_V))*EXP((-(t_ms))/((R_Ω)*(μC_F)))' "
    "}",

    "RL Transient",  "{ "
    "  '(I_A)=1/(R_Ω)*((Vf_V)-((Vf_V)-(Vi_V))*EXP((-(t_μs))/((R_Ω)*(L_mH)))' "
    "}",

    "Resonant Frequency",  "{ "
    "  'Qs=1/(R_Ω)*√((L_mH)/(C_μF))' "
	"  'Qp=(R_Ω)*√((C_μF)/(L_mH))' "
	"  '(ω0_(r/s))=2*(Ⓒπ_r)*(f0_Hz) "
	"  '(ω0_(r/s))=1_r/√((L_mH)*(C_μF))' "
    "}",

    "Plate Capacitor",  "{ "
    "  '(C_μF)=Ⓒε0*εr*(A_(cm^2))/(d_cm)' "
	"  '(ΔV_V)=(Ein_(V/m))*(d_cm)' "
    "  '(Ein_(N/C))=(σ_(μC/cm^2))/(Ⓒε0*εr)' "
	"  '(σ_(μC/m^2))=(Q_μC)/(A_(cm^2))' "
    "}",

    "Cylindrical Capacitor",  "{ "
    "  '(C_μF)=2*Ⓒπ*Ⓒε0*εr*(L_cm)/(LN((ro_cm)/(ri_cm)))' "
	"  '(ΔV_V)=(Q_μC)*(LN((ro_cm)/(ri_cm)))/(2*Ⓒπ*Ⓒε0*εr*(L_cm))' "
    "}",

    "Solenoid Inductance",  "{ "
    "  '(L_mH)=Ⓒμ0*μr*(n_cm^-1)^2*(A_(cm^2))*(h_cm)' "
    "}",

    "Toroid Inductance",  "{ "
    "  '(L_mH)=Ⓒμ0*μr*N^2*(h_cm)/(2*Ⓒπ)*LN((ro_cm)/(ri_cm))' "
    "}",

    "Sinusoidal Voltage",  "{ "
    "  '(V_V)=(Vmax_V)*SIN((ω_(r/s))*(t_μs)+(φ_°))' "
	"  '(ω_(r/s))=2*(Ⓒπ_r)*(f_Hz)' "
    "}",

    "Sinusoidal Current",  "{ "
    "  '(I_A)=(Imax_A)*SIN((ω_(r/s))*(t_s)+(φ_°))' "
	"  '(ω_(r/s))=2*(Ⓒπ_r)*(f_Hz)' "
    "}",

    // Example of the following in https://en.wikipedia.org/wiki/Drift_velocity#Numerical_example 
    "Drift Speed % Current Density",  "{ "
    "  '(vd_(m/s))=(I_A)/((n_(m^-3))*Ⓒqe*(A_(cm^2)))' "
    "  '(J_(A/m^2))=(vd_(m/s))*(ρ_(C/m^3)) "
    "  '(J_(A/m^2))=(σ_(S/m))*(E_(V_m))' "
    "}",
	
    // Example of the following in https://en.wikipedia.org/wiki/Electron_mobility#Examples 
    "Electron % Hole Mobilities",  "{ "
    "  '(J_(A/m^2))=(Je_(A/m^2))+(Jh_(A/m^2))' "
    "  '(Je_(A/m^2))=Ⓒqe*(ne_(m^-3))*(μe_(cm^2/(V*s)))*(E_(V/m)) "
    "  '(Jh_(A/m^2))=Ⓒqe*(nh_(m^-3))*(μh_(cm^2/(V*s)))*(E_(V/m)) "
    "  '(μe_(cm^2/(V*s)))=Ⓒqe*(τc_s)/(meeff_kg)' "
    "  '(μh_(cm^2/(V*s)))=Ⓒqe*(τc_s)/(mheff_kg)' "
    "  '(σ_(S/m)=Ⓒqe*((μe_(cm^2/(V*s)))*(ne_(m^-3))+(μh_(cm^2/(V*s)))*(nh_(m^-3)))' "
    "}",


    // ------------------------------------------------------------------------
    //   Physics
    // ------------------------------------------------------------------------

    "Phys",     nullptr,

    "RelativityMassEnergy",             "'(E_J)=(m_kg)*Ⓒc^2'",
    "IdealGas",                         "'(P_Pa)*(V_m^3)=(n_mol)*ⒸR*(T_K)'"
};
//   clang-format on


static runtime &invalid_equation_error()
// ----------------------------------------------------------------------------
//    Return the error message for invalid equations
// ----------------------------------------------------------------------------
{
    return rt.invalid_equation_error();
}


static symbol_p equation_label(symbol_r sym)
// ----------------------------------------------------------------------------
//   Simplify equations to show then in menu label
// ----------------------------------------------------------------------------
{
    if (sym)
    {
        size_t   len    = 0;
        utf8     source = sym->value(&len);
        if (object_p obj = object::parse(source, len))
            if (expression_p expr = obj->as<expression>())
                if (symbol_p ssym = expr->as_symbol(false))
                    return ssym;
    }
    return sym;
}


const equation::config equation::equations =
// ----------------------------------------------------------------------------
//  Define the configuration for the equations
// ----------------------------------------------------------------------------
{
    .menu_help      = "",
    .help           = "",
    .prefix         = L'Ⓔ',
    .type           = ID_equation,
    .first_menu     = ID_EquationsMenu00,
    .last_menu      = ID_EquationsMenu99,
    .name           = ID_EquationName,
    .value          = ID_EquationValue,
    .command        = ID_EquationSolver,
    .file           = "config/equations.csv",
    .builtins       = basic_equations,
    .nbuiltins      = sizeof(basic_equations) / sizeof(*basic_equations),
    .error          = invalid_equation_error,
    .label          = equation_label
};



// ============================================================================
//
//   Menu implementation
//
// ============================================================================

PARSE_BODY(equation)
// ----------------------------------------------------------------------------
//    Skip, the actual parsing is done in the symbol parser
// ----------------------------------------------------------------------------
{
    return do_parsing(equations, p);
}


EVAL_BODY(equation)
// ----------------------------------------------------------------------------
//   Equations always evaluate to their value
// ----------------------------------------------------------------------------
{
    object_g value = o->value();
    return rt.push(+value) ? OK : ERROR;
}


RENDER_BODY(equation)
// ----------------------------------------------------------------------------
//   Render the equation into the given buffer
// ----------------------------------------------------------------------------
{
    equation_g eq = o;
    do_rendering(equations, o, r);
    if (!r.editing() && Settings.ShowEquationBody())
    {
        if (object_g obj = eq->value())
        {
            r.put(':');
            obj->render(r);
        }
    }
    return r.size();
}


GRAPH_BODY(equation)
// ----------------------------------------------------------------------------
//   Render "normally"
// ----------------------------------------------------------------------------
{
    equation_g eq = o;
    if (Settings.ShowEquationBody())
    {
        if (object_g val = eq->value())
        {
            size_t namelen = 0;
            utf8 name = eq->name(&namelen);
            if (symbol_g namesym = symbol::make(name, namelen))
            {
                if (grob_g valg = val->graph(g))
                {
                    coord vv = g.voffset;
                    g.voffset = 0;
                    if (grob_g nameg = object::do_graph(+namesym, g))
                    {
                        coord nv = g.voffset;
                        g.voffset = 0;
                        grob_g r = expression::infix(g,
                                                     nv, nameg,
                                                     0, ":",
                                                     vv, valg);
                        return r;
                    }
                }
            }
        }
    }
    return object::do_graph(o, g);
}


HELP_BODY(equation)
// ----------------------------------------------------------------------------
//   Help topic for equations
// ----------------------------------------------------------------------------
{
    return o->do_instance_help(equation::equations);
}


MENU_BODY(equation_menu)
// ----------------------------------------------------------------------------
//   Build a equations menu
// ----------------------------------------------------------------------------
{
    return o->do_submenu(equation::equations, mi);
}


HELP_BODY(equation_menu)
// ----------------------------------------------------------------------------
//   Show the help for the given equation menu
// ----------------------------------------------------------------------------
{
    return o->do_menu_help(equation::equations, o);
}


MENU_BODY(EquationsMenu)
// ----------------------------------------------------------------------------
//   The equations menu is dynamically populated
// ----------------------------------------------------------------------------
{
    return equation::do_collection_menu(equation::equations, mi);
}


utf8 equation_menu::name(id type, size_t &len)
// ----------------------------------------------------------------------------
//   Return the name for a menu entry
// ----------------------------------------------------------------------------
{
    return do_name(equation::equations, type, len);
}


COMMAND_BODY(EquationName)
// ----------------------------------------------------------------------------
//   Put the name of a equation on the stack
// ----------------------------------------------------------------------------
{
    int key = ui.evaluating;
    if (constant_p cst = equation::do_key(equation::equations, key))
        if (equation_p eq = cst->as<equation>())
            if (rt.push(eq))
                return OK;
    if (!rt.error())
        rt.type_error();
    return ERROR;
}


INSERT_BODY(EquationName)
// ----------------------------------------------------------------------------
//   Put the name of a equation in the editor
// ----------------------------------------------------------------------------
{
    int key = ui.evaluating;
    return ui.insert_softkey(key, " Ⓔ", " ", false);
}


HELP_BODY(EquationName)
// ----------------------------------------------------------------------------
//   Put the help for a given equation name
// ----------------------------------------------------------------------------
{
    int key = ui.evaluating;
    if (constant_p cst = equation::do_key(equation::equations, key))
        if (equation_p eq = cst->as<equation>())
            return eq->help();
    return utf8("Equations");
}


COMMAND_BODY(EquationValue)
// ----------------------------------------------------------------------------
//   Put the value of a equation on the stack
// ----------------------------------------------------------------------------
{
    int key = ui.evaluating;
    if (constant_p cst = equation::do_key(equation::equations, key))
        if (equation_p eq = cst->as<equation>())
            if (object_p value = eq->value())
                if (rt.push(value))
                    return OK;
    if (!rt.error())
        rt.type_error();
    return ERROR;
}


INSERT_BODY(EquationValue)
// ----------------------------------------------------------------------------
//   Insert the value of a equation
// ----------------------------------------------------------------------------
{
    int key = ui.evaluating;
    if (object_p cstobj = equation::do_key(equation::equations, key))
        if (equation_p eq = cstobj->as<equation>())
            if (object_p value = eq->value())
                return ui.insert_object(value, " ", " ");
    return ERROR;
}


HELP_BODY(EquationValue)
// ----------------------------------------------------------------------------
//   Put the help for a given equation value
// ----------------------------------------------------------------------------
{
    return EquationName::do_help(nullptr);
}



COMMAND_BODY(EquationSolver)
// ----------------------------------------------------------------------------
//   Solve for a given equation
// ----------------------------------------------------------------------------
{
    int key = ui.evaluating;
    if (constant_p cst = equation::do_key(equation::equations, key))
        if (equation_p eq = cst->as<equation>())
            if (directory::store_here(static_object(ID_Equation), eq))
                if (const SolvingMenu *smenu =
                    object::static_object(ID_SolvingMenu)->as<SolvingMenu>())
                    return smenu->object::evaluate();
    if (!rt.error())
        rt.type_error();
    return ERROR;
}


INSERT_BODY(EquationSolver)
// ----------------------------------------------------------------------------
//   Insert the code in a program to solve a library equation
// ----------------------------------------------------------------------------
{
    int key = ui.evaluating;
    if (constant_p cst = equation::do_key(equation::equations, key))
        if (equation_p eq = cst->as<equation>())
            if (rt.push(eq))
                return ui.insert_object(eq, " ", " StEQ SolvingMenu ");
    return ERROR;
}


HELP_BODY(EquationSolver)
// ----------------------------------------------------------------------------
//   Put the help for a given equation value
// ----------------------------------------------------------------------------
{
    return EquationName::do_help(nullptr);
}


COMMAND_BODY(LibEq)
// ----------------------------------------------------------------------------
//   Evaluate a library equation
// ----------------------------------------------------------------------------
{
    return equation::lookup_command(equation::equations, false);
}
