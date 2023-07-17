// ****************************************************************************
//  complex.cc                                                    DB48X project
// ****************************************************************************
//
//   File Description:
//
//      Complex numbers
//
//      There are two representations for complex numbers:
//      - rectangular representation is one of X;Y, X+ⅈY, X-ⅈY, X+Yⅈ or X-Yⅈ
//      - polar representation is X∡Y
//
//      Some settings control how complex numbers are rendered
//
//
//
// ****************************************************************************
//   (C) 2023 Christophe de Dinechin <christophe@dinechin.org>
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

#include "complex.h"

#include "arithmetic.h"
#include "compare.h"
#include "functions.h"
#include "parser.h"
#include "renderer.h"
#include "runtime.h"


// ============================================================================
//
//   Generic complex operations
//
// ============================================================================
//
//  The generic operations optimize for the most efficient operation
//  if there is a difference between rectangular and polar
//

SIZE_BODY(complex)
// ----------------------------------------------------------------------------
//   Size of a complex number
// ----------------------------------------------------------------------------
{
    object_p p = object_p(payload(o));
    object_p e = p->skip()->skip();
    return byte_p(e) - byte_p(o);
}


algebraic_g complex::re() const
// ----------------------------------------------------------------------------
//   Return real part in a format-independent way
// ----------------------------------------------------------------------------
{
    if (type() == ID_rectangular)
        return rectangular_p(this)->re();
    return polar_p(this)->re();
}


algebraic_g complex::im() const
// ----------------------------------------------------------------------------
//   Return imaginary part in a format-independent way
// ----------------------------------------------------------------------------
{
    if (type() == ID_rectangular)
        return rectangular_p(this)->im();
    return polar_p(this)->im();
}


algebraic_g complex::mod() const
// ----------------------------------------------------------------------------
//   Return modulus in a format-independant way
// ----------------------------------------------------------------------------
{
    if (type() == ID_polar)
        return polar_p(this)->mod();
    return rectangular_p(this)->mod();
}


algebraic_g complex::arg() const
// ----------------------------------------------------------------------------
//   Return argument in a format-independant way
// ----------------------------------------------------------------------------
{
    if (type() == ID_polar)
        return polar_p(this)->arg();
    return rectangular_p(this)->arg();
}


complex_g complex::conjugate() const
// ----------------------------------------------------------------------------
//   Return complex conjugate in a format-independent way
// ----------------------------------------------------------------------------
{
    return make(type(), x(), -y());
}


complex_p complex::make(id type, algebraic_r x, algebraic_r y)
// ----------------------------------------------------------------------------
//   Build a complex of the right type
// ----------------------------------------------------------------------------
{
    if (!x.Safe() || !y.Safe())
        return nullptr;
    if (type == ID_polar)
        return polar::make(x, y);
    return rectangular::make(x, y);
}


rectangular_p complex::make(int re, int im)
// ----------------------------------------------------------------------------
//   Build a simple complex constant
// ----------------------------------------------------------------------------
{
    return rectangular_p(make(ID_rectangular,
                              integer::make(re), integer::make(im)));
}


complex_g operator-(complex_r x)
// ----------------------------------------------------------------------------
//  Unary minus
// ----------------------------------------------------------------------------
{
    if (x->type() == object::ID_polar)
    {
        polar_p p = polar_p(complex_p(x));
        return complex_p(polar::make(-p->mod(), p->arg()));
    }
    rectangular_p r = rectangular_p(complex_p(x));
    return rectangular::make(-r->re(), -r->im());
}


complex_g operator+(complex_r x, complex_r y)
// ----------------------------------------------------------------------------
//   Complex addition - Don't even bother doing it in polar form
// ----------------------------------------------------------------------------
{
    return rectangular::make(x->re() + y->re(), x->im() + y->im());
}


complex_g operator-(complex_r x, complex_r y)
// ----------------------------------------------------------------------------
//   Complex subtraction - Always in rectangular form
// ----------------------------------------------------------------------------
{
    return rectangular::make(x->re() - y->re(), x->im() - y->im());
}


complex_g operator*(complex_r x, complex_r y)
// ----------------------------------------------------------------------------
//   If both are in rectangular form, rectangular, otherwise polar
// ----------------------------------------------------------------------------
{
    object::id xt = x->type();
    object::id yt = y->type();
    if (xt != object::ID_rectangular || yt != object::ID_rectangular)
        return polar::make(x->mod() * y->mod(), x->arg() + y->arg());

    rectangular_p xx = rectangular_p(complex_p(x));
    rectangular_p yy = rectangular_p(complex_p(y));
    algebraic_g xr = xx->re();
    algebraic_g xi = xx->im();
    algebraic_g yr = yy->re();
    algebraic_g yi = yy->im();
    return rectangular::make(xr*yr-xi*yi, xr*yi+xi*yr);
}


complex_g operator/(complex_r x, complex_r y)
// ----------------------------------------------------------------------------
//   Like for multiplication, it's slighly cheaper in polar form
// ----------------------------------------------------------------------------
{
    object::id xt = x->type();
    object::id yt = y->type();
    if (xt != object::ID_rectangular || yt != object::ID_rectangular)
        return polar::make(x->mod() / y->mod(), x->arg() - y->arg());

    rectangular_p xx = rectangular_p(complex_p(x));
    rectangular_p yy = rectangular_p(complex_p(y));
    algebraic_g a = xx->re();
    algebraic_g b = xx->im();
    algebraic_g c = yy->re();
    algebraic_g d = yy->im();
    algebraic_g r = sq::run(c) + sq::run(d);
    return rectangular::make((a*c+b*d)/r, (b*c-a*d)/r);
}


polar_g complex::as_polar() const
// ----------------------------------------------------------------------------
//   Switch to polar form if preferred for computation
// ----------------------------------------------------------------------------
{
    if (type() == ID_rectangular)
    {
        rectangular_g r = rectangular_p(this);
        return polar::make(r->mod(), r->arg());
    }
    return polar_p(this);
}


rectangular_g complex::as_rectangular() const
// ----------------------------------------------------------------------------
//   Switch to rectangular form if preferred for computation
// ----------------------------------------------------------------------------
{
    if (type() == ID_polar)
    {
        polar_g r = polar_p(this);
        return rectangular::make(r->re(), r->im());
    }
    return rectangular_p(this);
}


PARSE_BODY(complex)
// ----------------------------------------------------------------------------
//   Parse the various forms of complex number
// ----------------------------------------------------------------------------
//   We accept the following formats:
//   a. (1;3)           Classic RPL
//   b. (1 3)           Classic RPL
//   c. 1ⅈ3             ⅈ as a separator
//   d. 1+ⅈ3            ⅈ as a prefix
//   e. 1-ⅈ3
//   f. 1+3ⅈ            ⅈ as a postfix
//   g. 1-3ⅈ
//   h. 1∡30            ∡ as a separator
//
//   Cases a-g generate a rectangular form, case i generates a polar form
//   Cases c-h can be surrounded by parentheses as well
//
//   In case (a), we do not accept (1,3) which classic RPL would accept,
//   because in DB48X 1,000.000 is a valid real number with thousands separator
{
    gcutf8      src    = p.source;
    size_t      max    = p.length;
    id          type   = ID_object;

    // Find the end of the possible complex number and check parentheses
    utf8        first  = src;
    utf8        last   = first;
    utf8        ybeg   = nullptr;
    size_t      xlen   = 0;
    size_t      ylen   = 0;
    bool        paren  = false;
    bool        signok = false;
    bool        ineq   = false;
    char        sign   = 0;
    while (size_t(last - first) < max)
    {
        unicode cp = utf8_codepoint(last);

        // Check if we have an opening parenthese
        if (last == first && cp == '(')
        {
            paren = true;
            first++;
        }

        // Check if found a '+' or '-' (cases d-g)
        else if (signok && (cp == '+' || cp == '-'))
        {
            if (sign)
            {
                // Cannot have two signs
                rt.syntax_error().source(last);
                return WARN;
            }
            sign = cp;
            ybeg = last + 1;
            if (type != ID_polar)
                xlen = last - first;
        }

        // Check if we found the ⅈ sign
        else if (cp == I_MARK)
        {
            // Can't have two complex signs
            if (type != ID_object)
            {
                rt.syntax_error().source(last);
                return WARN;
            }
            type = ID_rectangular;

            // Case of ⅈ as a separator (c)
            if (!sign)
            {
                ybeg = last + utf8_size(cp);
                xlen = last - first;
            }
            // Case of prefix ⅈ (d or e)
            else if (last == ybeg)
            {
                ybeg = last + utf8_size(cp);
            }
            // Case of postfix ⅈ (f or g)
            else
            {
                ylen = last - ybeg;
            }
        }

        // Check if we found the ∡ sign
        else if (cp == ANGLE_MARK)
        {
            // Can't have two complex signs, or have that with a sign
            if (type != ID_object || sign)
            {
                rt.syntax_error().source(last);
                return WARN;
            }
            type = ID_polar;

            // Case of ∡ as a separator (h)
            ybeg = last + utf8_size(cp);
            xlen = last - first;
        }

        // Check if we found a space or ';' inside parentheses
        else if (paren && (cp == ' ' || cp == ';'))
        {
            // Can't have two complex signs
            if (type != ID_object)
            {
                rt.syntax_error().source(last);
                return WARN;
            }
            type = ID_rectangular;
            ybeg = last + 1;
            xlen = last - first;
        }

        // Check if we found characters that we don't expect in a complex
        else if (cp == '"' || cp == '{' || cp == '[' || cp == L'«')
        {
            return SKIP;
        }

        // Check if we have equations in our complex
        else if (cp == '\'')
        {
            ineq = !ineq;
        }

        // Check if we have two parentheses
        else if (paren && !ineq && cp == '(')
        {
            rt.syntax_error().source(last);
            return WARN;
        }

        // Check if we found the end of the complex number
        else if (cp == ' ' || cp == '\n' || cp == '\t' ||
                 cp == ')' || cp == '}')
        {
            break;
        }

        // We can have a sign except after exponent markers
        signok = cp != 'e' && cp != 'E' && cp != L'⁳';

        // Loop on next characters
        last += utf8_size(cp);
    }

    // If we did not find the necessary structure, just skip
    if (type == ID_object || !xlen || !ybeg)
        return SKIP;

    // Check if we need to compute the length of y
    if (!ylen)
        ylen = last - ybeg;

    // Compute size that we parsed
    size_t parsed = last - first + paren;

    // Parse the first object
    gcutf8 ysrc = ybeg;
    size_t xsz = xlen;
    algebraic_g x = algebraic_p(object::parse(first, xlen));
    if (!x)
        return ERROR;
    if (xlen != xsz)
    {
        rt.syntax_error().source(utf8(src) + xlen);
        return ERROR;
    }

    // Parse the second object
    size_t ysz = ylen;
    algebraic_g y = algebraic_p(object::parse(ysrc, ylen));
    if (!y)
        return ERROR;
    if (ylen != ysz)
    {
        rt.syntax_error().source(utf8(ysrc) + ylen);
        return ERROR;
    }
    if (sign == '-')
    {
        y = neg::run(y);
        if (!y)
            return ERROR;
    }

    // Build the resulting complex
    complex_g result = complex::make(type, x, y);
    p.out = complex_p(result);
    p.end = parsed;

    return OK;
}



// ============================================================================
//
//   Specific code for rectangular form
//
// ============================================================================

algebraic_g rectangular::mod() const
// ----------------------------------------------------------------------------
//   Compute the modulus in rectangular form
// ----------------------------------------------------------------------------
{
    algebraic_g r = re();
    algebraic_g i = im();
    return hypot::evaluate(r, i);
}

algebraic_g rectangular::arg() const
// ----------------------------------------------------------------------------
//   Compute the argument in rectangular form
// ----------------------------------------------------------------------------
{
    algebraic_g r = re();
    algebraic_g i = im();
    return atan2::evaluate(i, r);
}


bool rectangular::is_zero() const
// ----------------------------------------------------------------------------
//   A complex in rectangular form is zero iff both re and im are zero
// ----------------------------------------------------------------------------
{
    return re()->is_zero(false) && im()->is_zero(false);
}


RENDER_BODY(rectangular)
// ----------------------------------------------------------------------------
//   Render a complex number in rectangular form
// ----------------------------------------------------------------------------
{
    algebraic_g re = o->re();
    algebraic_g im = o->im();
    re->render(r);
    r.need_sign();
    im->render(r);
    r.put(unicode(I_MARK));
    return r.size();
}



// ============================================================================
//
//   Polar-specific code
//
// ============================================================================

algebraic_g polar::re() const
// ----------------------------------------------------------------------------
//   Compute the real part in polar form
// ----------------------------------------------------------------------------
{
    algebraic_g m = mod();
    algebraic_g a = arg();
    return m * cos::run(a);
}

algebraic_g polar::im() const
// ----------------------------------------------------------------------------
//   Compute the imaginary part in polar form
// ----------------------------------------------------------------------------
{
    algebraic_g m = mod();
    algebraic_g a = arg();
    return m * sin::run(a);
}


bool polar::is_zero() const
// ----------------------------------------------------------------------------
//   A complex in polar form is zero iff modulus is zero
// ----------------------------------------------------------------------------
{
    return mod()->is_zero(false);
}


polar_p polar::make(algebraic_r mr, algebraic_r ar)
// ----------------------------------------------------------------------------
//   Build a normalized polar from given modulus and argument
// ----------------------------------------------------------------------------
{
    if (!mr.Safe() || !ar.Safe())
        return nullptr;
    algebraic_g pi = integer::make(4) * atan::run(integer::make(1));
    algebraic_g a = ar;
    algebraic_g m = mr;
    if (m->is_negative(false))
    {
        a = a + pi;
        m = neg::run(m);
    }
    if (a.Safe() && !a->is_strictly_symbolic())
    {
        if ((a < -pi)->as_truth(false) || (a > pi)->as_truth(false))
            a = (a+pi) % (pi+pi) - pi;
    }
    if (!a.Safe() || !m.Safe())
        return nullptr;
    return rt.make<polar>(m, a);
}


PARSE_BODY(polar)
// ----------------------------------------------------------------------------
//   Parse a complex number in polar form
// ----------------------------------------------------------------------------
{
    // This is really handled in the parser for 'rectangular'
    return SKIP;
}


RENDER_BODY(polar)
// ----------------------------------------------------------------------------
//   Render a complex number in rectangular form
// ----------------------------------------------------------------------------
{
    algebraic_g m = o->mod();
    algebraic_g a = o->arg();
    m->render(r);
    r.put(unicode(ANGLE_MARK));
    a->render(r);
    return r.size();
}


COMMAND_BODY(ImaginaryUnit)
// ----------------------------------------------------------------------------
//   Push a unit complex number on the stack
// ----------------------------------------------------------------------------
{
    object_g i = complex::make(0, 1);
    if (!i)
        return ERROR;
    if (!rt.push((i)))
        return ERROR;
    return OK;
}


COMMAND_BODY(RealToComplex)
// ----------------------------------------------------------------------------
//   Take two values in x and y and turn them into a complex
// ----------------------------------------------------------------------------
{
    object_g re = rt.stack(1);
    object_g im = rt.stack(0);
    if (!re || !im)
        return ERROR;
    if (!(re->is_real() || re->is_strictly_symbolic()) ||
        !(im->is_real() || im->is_strictly_symbolic()))
    {
        rt.type_error();
        return ERROR;
    }
    complex_g z = rectangular::make(algebraic_p(re.Safe()),
                                    algebraic_p(im.Safe()));
    if (!z.Safe() || !rt.drop())
        return ERROR;
    if (!rt.top(z))
        return ERROR;
    return OK;
}


COMMAND_BODY(ComplexToReal)
// ----------------------------------------------------------------------------
//   Take a complex value and convert it into two real values
// ----------------------------------------------------------------------------
{
    object_g z = rt.top();
    if (!z)
        return ERROR;
    if (!z->is_complex())
    {
        rt.type_error();
        return ERROR;
    }
    if (!rt.top(complex_p(z.Safe())->re()))
        return ERROR;
    if (!rt.push(object_p(complex_p(z.Safe())->im())))
        return ERROR;
    return OK;
}


COMMAND_BODY(ToRectangular)
// ----------------------------------------------------------------------------
//  Convert the top-level complex to rectangular form
// ----------------------------------------------------------------------------
{
    object_g x = rt.top();
    if (!x)
        return ERROR;
    if (!x->is_complex())
    {
        rt.type_error();
        return ERROR;
    }
    complex_g z = complex_p(x.Safe());
    if (z->type() == ID_polar)
    {
        z = rectangular::make(z->re(), z->im());
        if (!rt.push(object_p(complex_p(z.Safe()))))
            return ERROR;
    }
    return OK;
}


COMMAND_BODY(ToPolar)
// ----------------------------------------------------------------------------
//  Convert the top-level complex to polar form
// ----------------------------------------------------------------------------
{
    object_g x = rt.top();
    if (!x)
        return ERROR;
    if (!x->is_complex())
    {
        rt.type_error();
        return ERROR;
    }
    complex_g z = complex_p(x.Safe());
    if (z->type() == ID_rectangular)
    {
        z = polar::make(z->mod(), z->arg());
        if (!rt.push(object_p(complex_p(z.Safe()))))
            return ERROR;
    }
    return OK;
}



// ============================================================================
//
//   Implementation of complex functions
//
// ============================================================================

COMPLEX_BODY(sqrt)
// ----------------------------------------------------------------------------
//   Complex implementation of sqrt
// ----------------------------------------------------------------------------
{
    id zt = z->type();
    if (zt == ID_polar)
    {
        // Computation is a bit easier in polar form
        polar_r p = (polar_r) z;
        algebraic_g mod = p->mod();
        algebraic_g arg = p->arg();
        algebraic_g two = integer::make(2);
        return polar::make(sqrt::run(mod), arg / two);
    }

    rectangular_r r = (rectangular_r) z;
    algebraic_g a = r->re();
    algebraic_g b = r->im();
    algebraic_g znorm = abs::run(algebraic_p(z));
    algebraic_g two = algebraic_p(integer::make(2));
    algebraic_g re = sqrt::run((znorm + a) / two);
    algebraic_g im = sqrt::run((znorm - a) / two);
    if (b->is_negative(false))
        im = neg::run(im);
    else if (b->is_strictly_symbolic())
        im = sign::run(im) * im;
    return rectangular::make(re, im);
}


COMPLEX_BODY(cbrt)
// ----------------------------------------------------------------------------
//   Complex implementation of cbrt
// ----------------------------------------------------------------------------
{
    polar_g p = z->as_polar();
    if (!p.Safe())
        return nullptr;
    algebraic_g mod = p->mod();
    algebraic_g arg = p->arg();
    algebraic_g three = integer::make(3);
    return polar::make(cbrt::run(mod), arg / three);
}


COMPLEX_BODY(sin)
// ----------------------------------------------------------------------------
//   Complex implementation of sin
// ----------------------------------------------------------------------------
{
    // sin(z) = (exp(iz) - exp(-iz)) / 2i
    complex_g i = complex::make(0,1);
    complex_g iz = i * z;
    complex_g niz = -iz;
    iz = complex::exp(iz);
    niz = complex::exp(niz);
    return (iz - niz) / complex::make(0,2);
}

COMPLEX_BODY(cos)
// ----------------------------------------------------------------------------
//   Complex implementation of cos
// ----------------------------------------------------------------------------
{
    // cos(z) = (exp(iz) + exp(-iz)) / 2
    complex_g i = complex::make(0,1);
    complex_g iz = i * z;
    complex_g niz = -iz;
    iz = complex::exp(iz);
    niz = complex::exp(niz);
    return (iz + niz) / complex::make(2,0);
}


COMPLEX_BODY(tan)
// ----------------------------------------------------------------------------
//   Complex implementation of tan
// ----------------------------------------------------------------------------
{
    // tan(z) = -i * (exp(iz) - exp(-iz)) / (exp(iz) + exp(-iz))
    complex_g i = complex::make(0,1);
    complex_g iz = i * z;
    complex_g niz = -iz;
    iz = complex::exp(iz);
    niz = complex::exp(niz);
    return complex::make(0,-1) * (iz - niz) / (i + niz);
}


COMPLEX_BODY(asin)
// ----------------------------------------------------------------------------
//   Complex implementation of asin
// ----------------------------------------------------------------------------
{
    // asin(z) = i log(sqrt(1 - z^2) - iz)
    complex_g sq = z * z;
    complex_g one = complex::make(1,0);
    sq = complex::sqrt(one - sq);
    complex_g i = complex::make(0, 1);
    complex_g iz = i * z;
    return i * complex::log(sq - iz);
}


COMPLEX_BODY(acos)
// ----------------------------------------------------------------------------
//   Complex implementation of acos
// ----------------------------------------------------------------------------
{
    // acos(z) = -i log(z + i sqrt(1 - z^2))
    complex_g sq = z * z;
    complex_g one = complex::make(1,0);
    sq = complex::sqrt(one - sq);
    complex_g ni = complex::make(0,-1);
    return ni * complex::log(z - ni* sq);
}


COMPLEX_BODY(atan)
// ----------------------------------------------------------------------------
//   Complex implementation of atan
// ----------------------------------------------------------------------------
{
    // atan(z) = -i/2 ln((i-z) / (i + z))
    complex_g i = complex::make(0,1);
    return complex::log((i - z) / (i + z)) / complex_g(complex::make(0,2));
}


COMPLEX_BODY(sinh)
// ----------------------------------------------------------------------------
//   Complex implementation of sinh
// ----------------------------------------------------------------------------
{
    // sinh(z) = (exp(z) - exp(-z)) / 2
    return (complex::exp(z) - complex::exp(-z)) / complex_g(complex::make(2,0));
}

COMPLEX_BODY(cosh)
// ----------------------------------------------------------------------------
//   Complex implementation of cosh
// ----------------------------------------------------------------------------
{
    // cosh(z) = (exp(z) + exp(-z)) / 2
    complex_g two = complex::make(2,0);
    return (complex::exp(z) - complex::exp(-z)) / two;
}


COMPLEX_BODY(tanh)
// ----------------------------------------------------------------------------
//   Complex implementation of tanh
// ----------------------------------------------------------------------------
{
    // tanh(z) = (exp(2z) - 1) / (exp(2z) + 1)
    complex_g e2z = complex::exp(z + z);
    complex_g one = complex::make(1,0);
    return (e2z - one) /  (e2z + one);
}


COMPLEX_BODY(asinh)
// ----------------------------------------------------------------------------
//   Complex implementation of asinh
// ----------------------------------------------------------------------------
{
    // asinh(z) = ln(z + sqrt(z^2 + 1))
    complex_g one = complex::make(1, 0);
    return complex::log(z + complex::sqrt(z*z + one));
}


COMPLEX_BODY(acosh)
// ----------------------------------------------------------------------------
//   Complex implementation of acosh
// ----------------------------------------------------------------------------
{
    // asinh(z) = ln(z + sqrt(z^2 - 1))
    complex_g one = complex::make(1, 0);
    return complex::log(z + complex::sqrt(z*z - one));
}


COMPLEX_BODY(atanh)
// ----------------------------------------------------------------------------
//   Complex implementation of atanh
// ----------------------------------------------------------------------------
{
    // atanh(z) = 1/2 ln((1+z) / (1-z))
    complex_g one = complex::make(1, 0);
    complex_g two = complex::make(2, 0);
    return complex::log((one + z) / (one - z)) / two;
}


COMPLEX_BODY(log1p)
// ----------------------------------------------------------------------------
//   Complex implementation of log1p
// ----------------------------------------------------------------------------
{
    rt.type_error();
    return z;
}

COMPLEX_BODY(expm1)
// ----------------------------------------------------------------------------
//   Complex implementation of expm1
// ----------------------------------------------------------------------------
{
    rt.type_error();
    return z;
}


COMPLEX_BODY(log)
// ----------------------------------------------------------------------------
//   Complex implementation of log
// ----------------------------------------------------------------------------
{
    // log(a.exp(ib)) = log(a) + i b
    algebraic_g mod = z->mod();
    algebraic_g arg = z->arg();
    return rectangular::make(log::run(mod), arg);
}

COMPLEX_BODY(log10)
// ----------------------------------------------------------------------------
//   Complex implementation of log10
// ----------------------------------------------------------------------------
{
    algebraic_g ten = integer::make(10);
    algebraic_g zero = integer::make(0);
    complex_g logten = rectangular::make(log::run(ten), zero);
    return complex::log(z) / logten;
}


COMPLEX_BODY(log2)
// ----------------------------------------------------------------------------
//   Complex implementation of log2
// ----------------------------------------------------------------------------
{
    algebraic_g two = integer::make(2);
    algebraic_g zero = integer::make(0);
    complex_g logtwo = rectangular::make(log::run(two), zero);
    return complex::log(z) / logtwo;
}


COMPLEX_BODY(exp)
// ----------------------------------------------------------------------------
//   Complex implementation of exp
// ----------------------------------------------------------------------------
{
    // exp(a+ib) = exp(a)*exp(ib)
    algebraic_g re = z->re();
    algebraic_g im = z->im();
    return polar::make(exp::run(re), im);
}


COMPLEX_BODY(exp10)
// ----------------------------------------------------------------------------
//   Complex implementation of exp10
// ----------------------------------------------------------------------------
{
    algebraic_g ten = integer::make(10);
    algebraic_g zero = integer::make(0);
    complex_g logten = rectangular::make(log::run(ten), zero);
    return complex::exp(logten * z);
}


COMPLEX_BODY(exp2)
// ----------------------------------------------------------------------------
//   Complex implementation of exp2
// ----------------------------------------------------------------------------
{
    algebraic_g two = integer::make(2);
    algebraic_g zero = integer::make(0);
    complex_g logtwo = rectangular::make(log::run(two), zero);
    return complex::exp(logtwo * z);
}


COMPLEX_BODY(erf)
// ----------------------------------------------------------------------------
//   Complex implementation of erf
// ----------------------------------------------------------------------------
{
    rt.unimplemented_error();
    return z;
}

COMPLEX_BODY(erfc)
// ----------------------------------------------------------------------------
//   Complex implementation of erfc
// ----------------------------------------------------------------------------
{
    rt.unimplemented_error();
    return z;
}

COMPLEX_BODY(tgamma)
// ----------------------------------------------------------------------------
//   Complex implementation of tgamma
// ----------------------------------------------------------------------------
{
    rt.unimplemented_error();
    return z;
}

COMPLEX_BODY(lgamma)
// ----------------------------------------------------------------------------
//   Complex implementation of lgamma
// ----------------------------------------------------------------------------
{
    rt.unimplemented_error();
    return z;
}
