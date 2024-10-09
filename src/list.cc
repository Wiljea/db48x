// ****************************************************************************
//  list.cc                                                       DB48X project
// ****************************************************************************
//
//   File Description:
//
//     Implementation of RPL lists
//
//
//
//
//
//
//
//
// ****************************************************************************
//   (C) 2022 Christophe de Dinechin <christophe@dinechin.org>
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

#include "list.h"

#include "algebraic.h"
#include "array.h"
#include "compare.h"
#include "expression.h"
#include "grob.h"
#include "locals.h"
#include "parser.h"
#include "polynomial.h"
#include "precedence.h"
#include "program.h"
#include "renderer.h"
#include "runtime.h"
#include "symbol.h"
#include "unit.h"
#include "utf8.h"
#include "variables.h"

#include <stdio.h>
#include <stdlib.h>


RECORDER(list, 16, "Lists");
RECORDER(list_parse, 16, "List parsing");
RECORDER(list_error, 16, "Errors processing lists");

object::result list::list_parse(id      type,
                                parser &p,
                                unicode open,
                                unicode close)
// ----------------------------------------------------------------------------
//   Generic parser for sequences (list, program, etc)
// ----------------------------------------------------------------------------
//   We have to be careful here, because parsing sub-objects may allocate new
//   temporaries, which itself may cause garbage collection. So we need to keep
//   track of our current usage of the scratchpad to allow for recursive parsing
//   of complex objects, like { { A B { C D } } }
{
    // We have to be careful that we may have to GC to make room for list
    gcutf8   s           = p.source;
    size_t   max         = p.length;
    object_g infix       = nullptr;
    object_g prefix      = nullptr;
    object_g postfix     = nullptr;
    object_g obj         = nullptr;
    object_g xroot_arg   = nullptr;
    int      precedence  = p.precedence;
    int      lowest      = precedence;
    uint     arity       = 0;
    uint     arg         = 0;
    size_t   objcount    = 0;
    size_t   non_alg     = 0;
    size_t   non_alg_len = 0;
    bool     xroot       = false;

    // The IFTE command is special in that we don't evaluate its arguments
    bool        ifte     = false;

    // The `|` operator (where) is special when parsing parentheses
    bool        iswhere  = false;

    record(list, "Parse %lc%lc precedence %d length %u [%s]",
           open, close, precedence, max, utf8(s));

    // Check if we have the opening marker
    unicode cp  = 0;
    if (open)
    {
        cp = utf8_codepoint(s);
        if (cp != open)
            return SKIP;
        s = utf8_next(s);
        cp = 0;                 // Do not accept "'" as an empty equation
    }

    scribble scr;
    while (utf8_more(p.source, s, max))
    {
        cp = utf8_codepoint(s);
        if (cp == close)
        {
            s = utf8_next(s);
            break;
        }
        if (precedence && (cp == '\'' || cp == ')' || cp == ';'))
            break;
        if (utf8_whitespace(cp))
        {
            s = utf8_next(s);
            continue;
        }

        // Parse an object
        size_t done        = (byte *) s - (byte *) p.source;
        size_t length      = max > done ? max - done : 0;
        id     postfix_cmd = id(0);

        // For algebraic objects, check if we have or need parentheses
        if (precedence && length)
        {
            if (precedence > 0)
            {
                // Check if we see parentheses, or if we have `sin sin X`
                bool parenthese = (cp == '(' || arity > 1) && !infix;
                if (parenthese  || infix || prefix)
                {
                    int     childp  = infix      ? int(infix->precedence() + 1)
                                    : parenthese ? int(LOWEST)
                                                 : int(SYMBOL);

                    if (infix && cp == '(' && infix->type() == ID_Where)
                    {
                        iswhere = true;
                        if (utf8_more(p.source, s, max))
                            s = utf8_next(s);
                    }
                    if (iswhere)
                    {
                        parenthese = false;
                        childp = int(LOWEST);
                    }

                    parser  child(p, s, childp);
                    unicode iopen  = parenthese ? '(' : 0;
                    unicode iclose = parenthese ? ')' : 0;
                    id      ctype  = type == ID_unit ? ID_expression : type;

                    if (!infix && arity > 1)
                    {
                        if (arg)
                            iopen = 0;
                        arg++;
                        if (arg < arity)
                            iclose = ';';
                    }

                    record(list_parse, "%+s starting at offset %u '%s'",
                           parenthese ? "Parenthese" : "Child",
                           utf8(s) - utf8(p.source),
                           utf8(s));

                    auto rc = list_parse(ctype, child, iopen, iclose);
                    if (rc != OK)
                        return rc;
                    obj = child.out;
                    if (!obj)
                        return ERROR;

                    // A sum must have a name as first argument
                    if (arg == 1 && prefix)
                    {
                        id pty = prefix->type();
                        if (pty == ID_Sum || pty == ID_Product)
                        {
                            if (obj->as_quoted<symbol>() == nullptr)
                            {
                                rt.missing_variable_error()
                                    .source(+s+1, child.length-2);
                                return ERROR;
                            }
                        }
                    }

                    s = s + child.length;
                    record(list_parse,
                           "Child parsed as %t length %u",
                           object_p(obj), child.length);
                    precedence = -precedence; // Stay in postfix mode
                    cp = utf8_more(p.source, s, max) ? utf8_codepoint(s) : 0;
                    length = 0;

                    if (iswhere)
                    {
                        utf8 start = +s - child.length;
                        while (utf8_whitespace(cp) &&
                               utf8_more(p.source, s, max))
                        {
                            s = utf8_next(s);
                            cp = utf8_codepoint(s);
                        }
                        if (cp == ')')
                        {
                            iswhere = false;
                        }
                        else if (cp != ';')
                        {
                            rt.unterminated_error().source(start, +s-start);
                            return ERROR;
                        }
                        s = utf8_next(s);
                        cp = utf8_more(p.source,s,max) ? utf8_codepoint(s) : 0;
                    }
                }

                // Check to see if we have a sign
                else if (cp == '-' || cp == '+')
                {
                    if (cp == '+')
                    {
                        s = utf8_next(s);
                        continue;
                    }
                    obj = static_object(ID_neg);
                    length = 1;
                }

            }
            if (precedence < 0)
            {
                // Check special postfix notations
                switch(cp)
                {
                case L'²':
                    postfix_cmd = ID_sq;
                    break;
                case L'³':
                    postfix_cmd = ID_cubed;
                    break;
                case '!':
                    postfix_cmd = ID_fact;
                    break;
                case L'⁻':
                    if (utf8_codepoint(utf8_next(s)) == L'¹')
                        postfix_cmd = ID_inv;
                    break;
                default:
                    break;
                }
                if (postfix_cmd)
                {
                    utf8 cur = utf8(s);
                    (obj ? postfix : obj) = command::static_object(postfix_cmd);
                    length = postfix_cmd == ID_inv
                        ? utf8_next(utf8_next(cur)) - cur
                        : utf8_next(cur) - cur;
                    precedence = -precedence; // Stay in postfix mode
                }
            }
        }

        if (!obj && length)
        {
            obj = object::parse(s, length, precedence);
            record(list_parse,
                   "Item parsed as %t length %u arity %u",
                   object_p(obj), length, arity);
        }
        if (!obj)
            return ERROR;

        if (precedence && !postfix_cmd)
        {
            // We are parsing an equation
            if (precedence > 0)
            {
                // We just parsed an algebraic, e.g. 'sin', etc
                // stash it and require parentheses for arguments
                id type = obj->type();
                ifte = type == ID_IFTE;
                if (!is_algebraic(type) && !ifte)
                {
                    if (objcount)
                    {
                        rt.prefix_expected_error().source(s, length);
                        return ERROR;
                    }
                    non_alg = +s - +p.source;
                    non_alg_len = length;
                }

                // TODO: A symbol could be a function, need to deal with that
                if (is_command(type))
                {
                    prefix = obj;
                    if (prefix)
                    {
                        arity = prefix->arity();
                        arg = 0;

                        // xroot is the one command that parses args backwards
                        if (type == ID_xroot)
                            xroot = true;
                    }
                    if (arity)
                    {
                        obj = nullptr;
                        precedence = -SYMBOL;
                    }
                    else
                    {
                        prefix = nullptr;
                    }
                }
            }
            else if (int objprec = obj->precedence())
            {
                // We just parsed an infix, e.g. +, -, etc
                // stash it, or exit loop if it has lower precedence
                if (objprec < lowest)
                    break;
                if (objprec < FUNCTIONAL)
                {
                    infix = obj;
                    precedence = -objprec;
                    obj = nullptr;
                }
                else if (!infix)
                {
                    // Implicit multiplication
                    infix = static_object(ID_mul);
                }
            }
            else
            {
                precedence = -precedence;
            }
        }

        if (xroot && obj && !xroot_arg)
        {
            xroot_arg = obj;
            xroot = false;
            obj = nullptr;
            precedence = -precedence;
        }

        if (obj)
        {
            // Copy the parsed object to the scratch pad (may GC)
            do
            {
                record(list_parse, "Copying %t to scratchpad", object_p(obj));
                objcount++;

                size_t objsize = obj->size();

                // For equations, copy only the payload
                if (precedence && (!ifte || arg == 0))
                    if (expression_p eq = obj->as<expression>())
                        obj = eq->objects(&objsize);

                if (!rt.append(obj, objsize))
                    return ERROR;

                if (xroot_arg)
                {
                    obj = xroot_arg;
                    xroot_arg = nullptr;
                    arg++;
                }
                else if (prefix)
                {
                    if (arity > 1 && arg < arity)
                    {
                        precedence = -precedence;
                        break;
                    }
                    obj = prefix;
                    prefix = nullptr;
                }
                else if (postfix)
                {
                    obj = postfix;
                    postfix = nullptr;
                }
                else
                {
                    obj = infix;
                    infix = nullptr;
                }
            } while (obj);

            if (iswhere)
            {
                infix = object::static_object(ID_Where);
                precedence = -WHERE;
            }
        }

        // Jump past what we parsed
        s = s + length;

        // For equations switch between infix and prefix
        precedence = -precedence;
    }

    record(list, "Exiting parser at %s infix=%t prefix=%p",
           utf8(s), object_p(infix), object_p(prefix));


    // If we still have a pending opcode here, syntax error (e.g. '1+`)
    if (infix || prefix)
    {
        object_g obj = prefix ? prefix : infix;
        rt.command(obj);
        if (objcount == 0 && scr.growth() == 0)
        {
            if (!rt.append(obj))
                return ERROR;
            objcount = 1;
        }
        else
        {
            rt.argument_expected_error();
            return ERROR;
        }
    }

    if (non_alg && objcount != 1)
    {
        rt.syntax_error().source(p.source + non_alg, non_alg_len);
        return ERROR;
    }

    // Check that we have a matching closing character
    if (close && cp != close)
    {
        if (cp != ';')
            record(list_error,
                   "Missing terminator, got %u (%c) not %u (%c) at %s",
                   cp, cp, close, close, utf8(s));
        rt.unterminated_error().source(p.source, +s - +p.source);
        return ERROR;
    }

    // Create the object
    gcbytes scratch = scr.scratch();
    size_t  alloc   = scr.growth();
    size_t  parsed  = utf8(s) - utf8(p.source);

    // Check for the case of an empty equation
    if (alloc == 0 && type == ID_expression)
    {
        record(list_error, "Empty equation");
        rt.syntax_error().source(p.source);
        return ERROR;
    }

    // Create the object
    p.length = parsed;
    p.out    = rt.make<list>(type, scratch, alloc);

    record(list_parse, "Parsed as %t length %u", object_p(p.out), parsed);

    // Return success
    return p.out ? OK : ERROR;
}


intptr_t list::list_render(renderer &r, unicode open, unicode close) const
// ----------------------------------------------------------------------------
//   Render the list into the given buffer
// ----------------------------------------------------------------------------
{
    // Check if we need an indent in the body
    bool   unnest      = false;
    bool   first       = true;
    list_g list        = this;
    id     lty         = type();
    bool   need_indent = lty == ID_program;
    settings::SaveShowEquationBody sseb(false);

    for (object_p obj : *list)
    {
        id oty = obj->type();
        switch (oty)
        {
        case ID_array:
            if (first && lty == oty)
                unnest = true;
            // fallthrough
        case ID_list:
        case ID_program:
        case ID_locals:
        case ID_comment:
        case ID_IfThen:
        case ID_IfThenElse:
        case ID_DoUntil:
        case ID_WhileRepeat:
        case ID_StartStep:
        case ID_ForNext:
        case ID_ForStep:
        case ID_IfErrThen:
        case ID_IfErrThenElse:
            need_indent = true;
            break;
        default:
            break;
        }
        if (need_indent)
            break;
        first = false;
    }

    // Write the header, e.g. "{ "
    bool crpgm = need_indent && !unnest;
    if (open)
    {
        r.put(open);
        if (!unnest)
        {
            r.indent();
            r.wantSpace();
        }
    }
    if (crpgm)
        r.wantCR();

    // Loop on all objects inside the list
    for (object_p obj : *list)
    {
        id oty = obj->type();
        if (oty == ID_program || oty == ID_list || oty == ID_array)
        {
            if (open == 1 || !unnest)
                r.wantCR();
            if (unnest && open == 1)
                r.put("   ");
        }
        open = 1;

        // Render the object in what remains (may GC)
        if (oty != ID_array)
            r.wantSpace();
        if ((lty == ID_program || lty == ID_block) &&
            r.editing() &&
            Settings.VerticalProgramRendering())
            r.wantCR();
        obj->render(r);
        if (!unnest && oty != ID_array)
            r.wantSpace();
    }

    // Add final space and closing separator
    if (close)
    {
        if (!unnest)
        {
            r.unindent();
            if (lty != ID_array)
                r.wantSpace();
        }
        if (crpgm)
            r.wantCR();
        r.put(close);
        if (!unnest && lty != ID_array)
            r.wantSpace();
    }
    if (crpgm)
        r.wantCR();

    return r.size();
}


list_p list::append(list_p a) const
// ----------------------------------------------------------------------------
//   Append the second list to this one
// ----------------------------------------------------------------------------
{
    text_g x = text_p(this);
    text_g y = text_p(a);
    return list_p(+(x + y));
}


list_p list::append(object_p o) const
// ----------------------------------------------------------------------------
//   Append object to list
// ----------------------------------------------------------------------------
{
    text_g x = text_p(this);
    text_g y = text::make(byte_p(o), o->size());
    return list_p(+(x + y));
}


bool list::expand_without_size(size_t *size) const
// ----------------------------------------------------------------------------
//   Expand items on the stack, but do not add the size
// ----------------------------------------------------------------------------
{
    size_t depth = rt.depth();
    for (object_p obj : *this)
    {
        if (!rt.push(obj))
        {
            rt.drop(rt.depth() - depth);
            return false;
        }
    }
    if (size)
        *size = rt.depth() - depth;
    return true;
}


bool list::expand() const
// ----------------------------------------------------------------------------
//   Expand items on the stack, adding the size
// ----------------------------------------------------------------------------
{
    size_t depth = rt.depth();
    if (expand_without_size())
        if (rt.push(integer::make(rt.depth() - depth)))
            return true;
    rt.drop(rt.depth() - depth);
    return false;
}


bool list::expand_deep(uint32_t which) const
// ----------------------------------------------------------------------------
//   Expand list content, expending inner expressions/programs/lists
// ----------------------------------------------------------------------------
{
    for (object_p obj : *this)
    {
        id ty = obj->type();
        if (uint(ty) < 32 && ((1U<<uint(ty)) & which))
        {
            if (!list_p(obj)->expand_deep(which))
                return false;
        }
        else if (!rt.push(obj))
        {
            return false;
        }
    }
    return true;
}


PARSE_BODY(list)
// ----------------------------------------------------------------------------
//    Try to parse this as an list
// ----------------------------------------------------------------------------
{
    return list_parse(ID_list, p, '{', '}');
}


EVAL_BODY(list)
// ----------------------------------------------------------------------------
//   Evaluate a list if `ListEvaluation` is set
// ----------------------------------------------------------------------------
{
    if (!program::running && o->type() == ID_list && Settings.ListAsProgram())
        return program_p(o)->run(false);
    return object::do_evaluate(o);
}


RENDER_BODY(list)
// ----------------------------------------------------------------------------
//   Render the list into the given list buffer
// ----------------------------------------------------------------------------
{
    return o->list_render(r, '{', '}');
}


using pixsize = grob::pixsize;


static pixsize row_height(size_t r, size_t rows, size_t cols)
// ----------------------------------------------------------------------------
//   Compute the height of a row (max of height of all grobs)
// ----------------------------------------------------------------------------
{
    size_t  nitems = rows * cols;
    pixsize rh     = 0;
    for (size_t c = 0; c < cols; c++)
    {
        size_t i       = r * cols + c;
        grob_p colitem = grob_p(rt.stack(nitems + ~i));
        ASSERT(colitem);
        pixsize h = colitem->height();
        if (rh < h)
            rh = h;
    }
    return rh;
}


static pixsize col_width(size_t c, size_t rows, size_t cols)
// ----------------------------------------------------------------------------
//   Compute the height of a row (max of height of all grobs)
// ----------------------------------------------------------------------------
{
    size_t  nitems = rows * cols;
    pixsize cw     = 0;
    for (size_t r = 0; r < rows; r++)
    {
        size_t  i       = r * cols + c;
        grob_p  colitem = grob_p(rt.stack(nitems + ~i));
        pixsize w       = colitem->width();
        if (cw < w)
            cw = w;
    }
    return cw;
}


grob_p list::graph(grapher &g, size_t rows, size_t cols, bool mat) const
// ----------------------------------------------------------------------------
//   Render a list or array exploded on stack
// ----------------------------------------------------------------------------
{
    list_g list = this;
    save<bool> sgraph(g.graph, false);

    // Convert all elements to graphical equivalent
    size_t nitems = rows * cols;
    for (size_t i = 0; i < nitems; i++)
    {
        object_p item = rt.stack(i);
        grob_g   grob = item->graph(g);
        if (!grob || grob->type() != ID_grob)
        {
            // Ran into a problem with one rendering, e.g. out of memory
            // Fallback to rendering as text.
            record(list_error, "Problem graphing %zu in %zu x %zu",
                   i, rows, cols);
            rt.drop(nitems);
            return object::do_graph(list, g);
        }
        rt.stack(i, grob);
    }

    // Compute the width, starting with 8 pixels for borders
    pixsize bw  = mat ? 4 : 2;
    pixsize sw  = 2;
    pixsize gap = 12;
    pixsize gw  = 4 * bw + (cols - 1) * gap + 2 * sw;
    pixsize gh  = 0;

    // Compute the width as the sum of the max width of all columns
    for (size_t c = 0; c < cols; c++)
    {
        size_t cw = col_width(c, rows, cols);
        gw += cw;
        if (gw > g.maxw)
        {
            rt.drop(nitems);
            return nullptr;
        }
    }

    // Compute the height as the sum of the max height of all rows
    for (size_t r = 0; r < rows; r++)
    {
        size_t rh = row_height(r, rows, cols);
        gh += rh;
        if (gh > g.maxh)
        {
            rt.drop(nitems);
            return nullptr;
        }
    }

    // Create the resulting graph
    grob_g result = g.grob(gw, gh);
    if (!result)
    {
        rt.drop(nitems);
        return nullptr;
    }
    grob::surface rs = result->pixels();
    rs.fill(0, 0, gw, gh, g.background);

    // Position of items in enclosing border and adjust
    coord xl = 0;
    coord yt = 0;

    // Draw all items inside the matrix
    coord yi = yt;
    for (size_t r = 0; r < rows; r++)
    {
        pixsize rh = row_height(r, rows, cols);
        coord   xi = xl + 2 * bw + sw;
        for (size_t c = 0; c < cols; c++)
        {
            pixsize cw      = col_width(c, rows, cols);
            size_t  i       = r * cols + c;
            grob_p  colitem = grob_p(rt.stack(nitems + ~i));
            grob::surface is      = colitem->pixels();
            pixsize iw      = is.width();
            pixsize ih      = is.height();
            rs.copy(is, xi + (cw - iw)/2, yi + (rh - ih)/2);
            xi += cw + gap;
        }
        yi += rh;
    }

    rt.drop(nitems);
    return result;
}


GRAPH_BODY(list)
// ----------------------------------------------------------------------------
//  Render the list vertically if the option is given
// ----------------------------------------------------------------------------
{
    if (Settings.VerticalLists() && o->type() == ID_list)
    {
        size_t depth = rt.depth();
        list_g list = o;
        if (list->expand_without_size())
        {
            size_t rows = rt.depth() - depth;
            size_t cols = 1;
            grob_g result = list->graph(g, rows, cols, true);
            if (!result)
                return nullptr;
            grob::surface rs = result->pixels();
            pixsize       gw = rs.width();
            pixsize       gh = rs.height();

            // Dimensions of border and adjust
            coord   xl = 0;
            coord   xr = coord(gw) - 2;
            coord   yt = 0;
            coord   yb = coord(gh) - 4;
            pixsize bw = 2;

            // Add the borders
            for (coord y = 1; y < yb; y++)
            {
                coord d = y < yb/2 - 1 || y > yb/2 + 1;
                rs.fill(xl+d,    y, xl+d+bw,   y, g.foreground);
                rs.fill(xr-d-bw, y, xr-d,      y, g.foreground);
            }
            rs.fill(xl+bw,   yt, xl+2*bw, yt+1, g.foreground);
            rs.fill(xr-2*bw, yt, xr-bw,   yt+1, g.foreground);
            rs.fill(xl+bw,   yb, xl+2*bw, yb+1, g.foreground);
            rs.fill(xr-2*bw, yb, xr-bw,   yb+1, g.foreground);

            return result;
        }
    }

    return object::do_graph(o, g);
}


HELP_BODY(list)
// ----------------------------------------------------------------------------
//   Help topic for lists
// ----------------------------------------------------------------------------
{
    return utf8("Lists");
}



// ============================================================================
//
//   Command implementation
//
// ============================================================================

list_p to_list_object(uint depth)
// ----------------------------------------------------------------------------
//  Make a list from the stack as an object
// ----------------------------------------------------------------------------
{
    scribble scr;
    for (uint i = 0; i < depth; i++)
        if (object_g obj = rt.stack(depth - 1 - i))
            if (!rt.append(obj))
                return nullptr;

    if (list_p result = list::make(scr.scratch(), scr.growth()))
    {
        rt.drop(depth);
        return result;
    }
    return nullptr;
}


object::result to_list(uint depth)
// ----------------------------------------------------------------------------
//  Make a list on the stack
// ----------------------------------------------------------------------------
{
    if (object_g list = to_list_object(depth))
        if (rt.push(list))
            return object::OK;
    return object::ERROR;
}


COMMAND_BODY(ToList)
// ----------------------------------------------------------------------------
//   Convert elements to a list
// ----------------------------------------------------------------------------
{
    uint32_t depth = uint32_arg();
    if (!rt.error())
    {
        if (!rt.args(depth + 1))
            return ERROR;

        if (rt.pop())
            return to_list(depth);
    }
    return ERROR;
}


COMMAND_BODY(FromList)
// ----------------------------------------------------------------------------
//   Convert elements to a list
// ----------------------------------------------------------------------------
{
    object_p obj = rt.top();
    if (list_p li = obj->as<list>())
    {
        rt.drop();
        if (li->expand())
            return OK;
    }
    else
    {
        rt.type_error();
    }
    return ERROR;
}


COMMAND_BODY(Size)
// ----------------------------------------------------------------------------
//   Return the size of an object
// ----------------------------------------------------------------------------
//   Behaves differently from standard RPL for integers, equations and units
//   where it returns 1 and not some weirdo internal count
{
    object_p obj    = rt.top();
    id       oty    = obj->type();
    size_t   size   = 1;

    switch (oty)
    {
    case ID_list:
        size = list_p(obj)->items(); break;
    case ID_array:
        if (object_p result = array_p(obj)->dimensions())
            if (rt.top(result))
                return OK;
        break;
    case ID_text:
        size = text_p(obj)->utf8_characters(); break;
    case ID_grob:
    case ID_bitmap:
        if (grob_p gr = grob_p(obj))
        {
            grob::pixsize w = 0, h = 0;
            if (gr->pixels(&w, &h))
            {
                integer_g wo = rt.make<based_integer>(w);
                integer_g ho = rt.make<based_integer>(h);
                if (wo && ho && rt.top(+wo) && rt.push(+ho))
                    return OK;
            }
        }
        return ERROR;
    default:
        break;
    }

    integer_p szo = integer::make(size);
    if (szo && rt.top(szo))
        return OK;
    return ERROR;
}


static object::result get(bool increment)
// ----------------------------------------------------------------------------
//   Get element from structure, incrementing index or not
// ----------------------------------------------------------------------------
{
    // Check we have an object at level 2
    if (object_p items = rt.stack(1))
    {
        if (symbol_p name = items->as_quoted<symbol>())
        {
            items = directory::recall_all(name, true);
            if (!items)
                return object::ERROR;
        }

        object_p item = items->at(rt.stack(0));
        if (!item)
        {
            if (!rt.error())
                rt.index_error();
        }
        else if (increment)
        {
            rt.push(item);
            object_g index = rt.stack(1);
            bool wrap = items->next_index(&+index);
            if (index)
            {
                rt.stack(1, index);
                Settings.IndexWrapped(wrap);
                return object::OK;
            }
        }
        else if (rt.pop() && rt.top(item))
        {
            return object::OK;
        }
    }
    return object::ERROR;
}


COMMAND_BODY(Get)
// ----------------------------------------------------------------------------
//   Get an element in a list
// ----------------------------------------------------------------------------
{
    return get(false);
}


COMMAND_BODY(GetI)
// ----------------------------------------------------------------------------
//   Get an element in a list and increment the index
// ----------------------------------------------------------------------------
{
    return get(true);
}


static object::result put(bool increment)
// ----------------------------------------------------------------------------
//   Put element in structure, incrementing index or not
// ----------------------------------------------------------------------------
{
    // Check that we have an object at level 2
    if (object_p items = rt.stack(2))
    {
        symbol_p name = items->as_quoted<symbol>();
        if (name)
        {
            items = directory::recall_all(name, true);
            if (!items)
                return object::ERROR;
        }

        if (object_g result = items->at(rt.stack(1), rt.top()))
        {
            if (increment)
            {
                object_g index = rt.stack(1);
                bool wrap = result->next_index(&+index);
                if (index)
                {
                    rt.stack(1, +index);
                    Settings.IndexWrapped(wrap);
                }
            }
            if (name)
            {
                name = rt.stack(2)->as_quoted<symbol>();
                if (directory::update(name, result))
                {
                    rt.drop(increment ? 1 : 3);
                    return object::OK;
                }
            }
            else
            {
                if (rt.drop(increment ? 1 : 2))
                    if (rt.stack(increment ? 1 : 0, result))
                        return object::OK;
            }
        }

        if(!rt.error())
            rt.index_error();
    }
    return object::ERROR;
}


COMMAND_BODY(Put)
// ----------------------------------------------------------------------------
//   Put an element in a list
// ----------------------------------------------------------------------------
{
    return put(false);
}


COMMAND_BODY(PutI)
// ----------------------------------------------------------------------------
//   Put an element in a list, incrementing the index
// ----------------------------------------------------------------------------
{
    return put(true);
}


COMMAND_BODY(Head)
// ----------------------------------------------------------------------------
//   Return first element in a list
// ----------------------------------------------------------------------------
{
    object_p obj = rt.top();
    id ty = obj->type();
    if (ty == ID_list || ty == ID_array)
    {
        if (object_p hd = list_p(obj)->head())
        {
            if (rt.top(hd))
                return OK;
        }
        else
        {
            rt.dimension_error();
        }
    }
    else
    {
        rt.type_error();
    }
    return ERROR;
}


COMMAND_BODY(Tail)
// ----------------------------------------------------------------------------
//   Return all but first element in a list
// ----------------------------------------------------------------------------
{
    object_p obj = rt.top();
    id ty = obj->type();
    if (ty == ID_list || ty == ID_array)
    {
        if (object_p tl = list_p(obj)->tail())
        {
            if (rt.top(tl))
                return OK;
        }
        else
        {
            // On HP50, HEAD returns an error, but TAIL returns an empty list
            // rt.dimension_error();
            return OK;
        }
    }
    else
    {
        rt.type_error();
    }
    return ERROR;
}


static object::result map_reduce_filter(object_p (list::*cmd)(object_p) const)
// ----------------------------------------------------------------------------
//   Shared code for map, reduce and filter
// ----------------------------------------------------------------------------
{
    size_t depth = rt.depth();
    object_p   obj = rt.stack(1);
    object_g   prg = rt.top();
    object::id ty  = obj->type();
    if (ty == object::ID_list || ty == object::ID_array)
    {
        object_p result = (list_p(obj)->*cmd)(prg);
        if (!result)
            goto error;
        if (rt.drop() && rt.top(result))
            return object::OK;
    }
    else
    {
        rt.type_error();
    }
error:
    if (rt.depth() > depth)
        rt.drop(rt.depth() - depth);
    return object::ERROR;
}


COMMAND_BODY(Map)
// ----------------------------------------------------------------------------
//   Apply unary function in level 1 to all elements in level 2
// ----------------------------------------------------------------------------
{
    return map_reduce_filter(&list::map_as_object);
}


COMMAND_BODY(Reduce)
// ----------------------------------------------------------------------------
//   Apply binary function in level 1 pairwise to combine elements in level 2
// ----------------------------------------------------------------------------
{
    return map_reduce_filter(&list::reduce);
}


COMMAND_BODY(Filter)
// ----------------------------------------------------------------------------
//   Filter the function in level 1 to all elements in level 2
// ----------------------------------------------------------------------------
{
    return map_reduce_filter(&list::filter_as_object);
}


static object::result list_reduce(object::id cmd)
// ----------------------------------------------------------------------------
//   Shared code for map, reduce and filter
// ----------------------------------------------------------------------------
{
    object_p   obj = rt.stack(0);
    object::id ty  = obj->type();
    if (ty == object::ID_list || ty == object::ID_array)
    {
        object_p result = list_p(obj)->reduce(command::static_object(cmd));
        if (result && rt.top(result))
            return object::OK;
    }
    else
    {
        rt.type_error();
    }
    return object::ERROR;
}


static object::result list_pair_map(object::id cmd)
// ----------------------------------------------------------------------------
//   Shared code for map, reduce and filter
// ----------------------------------------------------------------------------
{
    object_p   obj = rt.stack(0);
    object::id ty  = obj->type();
    if (ty == object::ID_list || ty == object::ID_array)
    {
        object_p cmdobj = command::static_object(cmd);
        object_p result = list_p(obj)->pair_map(cmdobj);
        if (result && rt.top(result))
            return object::OK;
    }
    else
    {
        rt.type_error();
    }
    return object::ERROR;
}


COMMAND_BODY(ListSum)
// ----------------------------------------------------------------------------
//   Return the sum of a list or array
// ----------------------------------------------------------------------------
{
    return list_reduce(ID_add);
}


COMMAND_BODY(ListProduct)
// ----------------------------------------------------------------------------
//   Return the product of a list or array
// ----------------------------------------------------------------------------
{
    return list_reduce(ID_mul);
}


COMMAND_BODY(ListDifferences)
// ----------------------------------------------------------------------------
//   Return the product of a list or array
// ----------------------------------------------------------------------------
{
    return list_pair_map(ID_sub);
}


object_p list::head() const
// ----------------------------------------------------------------------------
//   Return the first element in the list
// ----------------------------------------------------------------------------
{
    size_t   size  = 0;
    object_p first = objects(&size);
    if (!size)
        return nullptr;
    return first;
}


list_p list::tail() const
// ----------------------------------------------------------------------------
//   Return the tail elements of the list
// ----------------------------------------------------------------------------
{
    size_t   size  = 0;
    object_p first = objects(&size);
    if (!size)
        return nullptr;
    size_t osize = first->size();
    byte_p rest  = byte_p(first) + osize;
    return list::make(type(), rest, size - osize);
}


list_p list::map(object_p prgobj) const
// ----------------------------------------------------------------------------
//   Apply an RPL object (nominally a program) on all elements in the list
// ----------------------------------------------------------------------------
{
    id       ty    = type();
    object_g prg   = prgobj;
    size_t   depth = rt.depth();
    scribble scr;
    for (object_p obj : *this)
    {
        id oty = obj->type();
        if (oty == ID_array || oty == ID_list)
        {
            list_g sub = list_p(obj)->map(prg);
            obj = +sub;
        }
        else
        {
            if (!rt.push(obj))
                goto error;
            if (program::run(prg, true) != OK)
                goto error;
            if (rt.depth() != depth + 1)
            {
                rt.misbehaving_program_error();
                goto error;
            }
            obj = rt.pop();
        }
        if (!obj)
            goto error;

        if (!rt.append(obj))
            goto error;
    }

    return list::make(ty, scr.scratch(), scr.growth());

error:
    if (rt.depth() > depth)
        rt.drop(rt.depth() - depth);
    return nullptr;

}


object_p list::reduce(object_p prgobj) const
// ----------------------------------------------------------------------------
//   Apply an RPL object (nominally a program) on pairs of list elements
// ----------------------------------------------------------------------------
{
    object_g prg    = prgobj;
    size_t   depth  = rt.depth();
    object_g result = nullptr;
    for (object_p obj : *this)
    {
        if (!rt.push(obj))
            goto error;
        if (!result)
        {
            result = obj;
        }
        else
        {
            if (program::run(prg, true) != OK)
                goto error;
            if (rt.depth() != depth + 1)
                rt.misbehaving_program_error();
            result = rt.top();
        }
        if (rt.error())
            goto error;
    }
    if (rt.depth() > depth)
        rt.drop(rt.depth() - depth);
    return result;

error:
    if (rt.depth() > depth)
        rt.drop(rt.depth() - depth);
    return nullptr;
}


list_p list::filter(object_p prgobj) const
// ----------------------------------------------------------------------------
//   Apply an RPL object (nominally a program) to filter elements in a list
// ----------------------------------------------------------------------------
{
    id       ty    = type();
    object_g prg   = prgobj;
    size_t   depth = rt.depth();
    scribble scr;
    for (object_g obj : *this)
    {
        id   oty  = obj->type();
        bool keep = false;
        if (oty == ID_array || oty == ID_list)
        {
            object_g sub = list_p(+obj)->filter(prg);
            obj = +sub;
            keep = true;
        }
        else
        {
            if (!rt.push(obj))
                goto error;
            if (program::run(prg, true) != OK)
                goto error;
            if (rt.depth() != depth + 1)
            {
                rt.misbehaving_program_error();
                goto error;
            }
            object_p test = rt.pop();
            keep = test->as_truth(true);
            if (rt.error())
                goto error;
        }

        if (keep && !rt.append(obj))
            goto error;
    }

    return list::make(ty, scr.scratch(), scr.growth());

error:
    if (rt.depth() > depth)
        rt.drop(rt.depth() - depth);
    return nullptr;
}


list_p list::pair_map(object_p prgobj) const
// ----------------------------------------------------------------------------
//   Apply an RPL object (nominally a program) to combine successive elements
// ----------------------------------------------------------------------------
{
    id       ty    = type();
    object_g prg   = prgobj;
    size_t   depth = rt.depth();
    object_g prev;
    scribble scr;
    for (object_g obj : *this)
    {
        if (prev)
        {
            if (!rt.push(obj) || !rt.push(prev))
                goto error;
            if (program::run(prg, true) != OK)
                goto error;
            if (rt.depth() != depth + 1)
            {
                rt.misbehaving_program_error();
                goto error;
            }
            object_p item   = rt.pop();
            if (!rt.append(item))
                goto error;
        }
        prev = obj;
    }

    return list::make(ty, scr.scratch(), scr.growth());

error:
    if (rt.depth() > depth)
        rt.drop(rt.depth() - depth);
    return nullptr;
}


list_p list::map(algebraic_fn fn) const
// ----------------------------------------------------------------------------
//   Apply an algebraic function on all elements in the list
// ----------------------------------------------------------------------------
{
    id ty = type();
    scribble scr;
    for (object_p obj : *this)
    {
        id oty = obj->type();
        if (oty == ID_array || oty == ID_list)
        {
            list_g sub = list_p(obj)->map(fn);
            obj = +sub;
        }
        else
        {
            algebraic_g a = obj->as_algebraic();
            if (!a)
            {
                rt.type_error();
                return nullptr;
            }

            a = fn(a);
            if (!a)
                return nullptr;
            obj = +a;
        }

        if (!rt.append(obj))
            return nullptr;
    }

    return list::make(ty, scr.scratch(), scr.growth());
}


list_p list::map(arithmetic_fn fn, algebraic_r y) const
// ----------------------------------------------------------------------------
//   Right-apply an arithmtic function on all elements in the list
// ----------------------------------------------------------------------------
{
    id ty = type();
    scribble scr;
    for (object_p obj : *this)
    {
        id oty = obj->type();
        if (oty == ID_array || oty == ID_list)
        {
            list_g sub = list_p(obj)->map(fn, y);
            obj = +sub;
        }
        else
        {
            algebraic_g a = obj->as_algebraic();
            if (!a)
            {
                rt.type_error();
                return nullptr;
            }

            a = fn(a, y);
            if (!a)
                return nullptr;
            obj = +a;
        }

        if (!rt.append(obj))
            return nullptr;
    }

    return list::make(ty, scr.scratch(), scr.growth());
}


list_p list::map(algebraic_r x, arithmetic_fn fn) const
// ----------------------------------------------------------------------------
//   Left-apply an arithmtic function on all elements in the list
// ----------------------------------------------------------------------------
{
    id ty = type();
    scribble scr;
    for (object_p obj : *this)
    {
        id oty = obj->type();
        if (oty == ID_array || oty == ID_list)
        {
            list_g sub = list_p(obj)->map(x, fn);
            obj = +sub;
            if (!obj)
                return nullptr;
        }
        else
        {
            algebraic_g a = obj->as_algebraic();
            if (!a)
            {
                rt.type_error();
                return nullptr;
            }

            a = fn(x, a);
            if (!a)
                return nullptr;
            obj = +a;
        }

        if (!rt.append(obj))
            return nullptr;
    }

    return list::make(ty, scr.scratch(), scr.growth());
}



// ============================================================================
//
//   Sorting
//
// ============================================================================

int memory_compare(object_p *xp, object_p *yp)
// ----------------------------------------------------------------------------
//   Compare using type IDs and memory content
// ----------------------------------------------------------------------------
{
    object_p x = *xp;
    object_p y = *yp;
    return x->compare_to(y);
}


int value_compare(object_p *xp, object_p *yp)
// ----------------------------------------------------------------------------
//   Sort items according to value
// ----------------------------------------------------------------------------
{
    object_p x = *xp;
    object_p y = *yp;
    object::id xty = x->type();
    object::id yty = y->type();
    if ((object::is_algebraic(xty) && object::is_algebraic(yty)) ||
        (xty == object::ID_array   && yty == object::ID_array) ||
        (xty == object::ID_list    && yty == object::ID_list))
    {
        algebraic_g xa     = algebraic_p(x);
        algebraic_g ya     = algebraic_p(y);
        int         result = 0;
        if (comparison::compare(&result, xa, ya))
            return result;
    }
    return x->compare_to(y);
}


static int value_compare_reverse(object_p *xp, object_p *yp)
// ----------------------------------------------------------------------------
//   Sort item according in decreasing value order
// ----------------------------------------------------------------------------
{
    return -value_compare(xp, yp);
}


static int memory_compare_reverse(object_p *xp, object_p *yp)
// ----------------------------------------------------------------------------
//   Sort item according in decreasing value order
// ----------------------------------------------------------------------------
{
    return -memory_compare(xp, yp);
}


static object::result do_sort(int (*compare)(object_p *x, object_p *y))
// ----------------------------------------------------------------------------
//   RPL command for a sort
// ----------------------------------------------------------------------------
{
    typedef int (*qsort_fn)(const void *, const void*);

    if  (object_p obj = rt.stack(0))
    {
        object::id oty = obj->type();
        if (oty == object::ID_list || oty == object::ID_array)
        {
            size_t   depth = rt.depth();
            list_g   items = list_p(obj);
            size_t   count;
            scribble scr;
            qsort_fn cmp = qsort_fn(compare);

            for (object_p item : *items)
                if (!rt.push(item))
                    goto err;
            count = rt.depth() - depth;
            if (cmp)
                qsort(rt.stack_base(), count, sizeof(object_p), cmp);

            for (uint i = 0; i < count; i++)
                if (object_g obj = rt.stack(i))
                    if (!rt.append(obj))
                        goto err;
            rt.drop(count);
            items = list::make(oty, scr.scratch(), scr.growth());
            if (items && rt.top(+items))
                return object::OK;

        err:
            rt.drop(rt.depth() - depth);
            return object::ERROR;
        }
        else
        {
            rt.type_error();
        }
    }
    return object::ERROR;
}


COMMAND_BODY(Sort)
// ----------------------------------------------------------------------------
//   Sort contents of a list according to value
// ----------------------------------------------------------------------------
{
    return do_sort(value_compare);
}


COMMAND_BODY(QuickSort)
// ----------------------------------------------------------------------------
//   Sort contents of a list using memory comparisons
// ----------------------------------------------------------------------------
{
    return do_sort(memory_compare);
}


COMMAND_BODY(ReverseSort)
// ----------------------------------------------------------------------------
//   Sort contents of a list according to value
// ----------------------------------------------------------------------------
{
    return do_sort(value_compare_reverse);
}


COMMAND_BODY(ReverseQuickSort)
// ----------------------------------------------------------------------------
//   Sort contents of a list using memory comparisons
// ----------------------------------------------------------------------------
{
    return do_sort(memory_compare_reverse);
}


COMMAND_BODY(ReverseList)
// ----------------------------------------------------------------------------
//   Reverse a list
// ----------------------------------------------------------------------------
{
    return do_sort(nullptr);
}



// ============================================================================
//
//   Enumeration of names in an expression, a program or a list
//
// ============================================================================

bool list::names_enumerate(size_t depth, bool units) const
// ----------------------------------------------------------------------------
//   Enumerate the names in the object
// ----------------------------------------------------------------------------
{
    size_t vars;
    switch(type())
    {
    case ID_expression:
    case ID_funcall:
    case ID_program:
    case ID_list:
    case ID_block:
    case ID_array:
        for (auto obj : *this)
        {
            switch(obj->type())
            {
            case ID_symbol:
                if (!names_insert(depth, symbol_p(obj), nullptr))
                    return false;
                break;
            case ID_expression:
            case ID_funcall:
            case ID_polynomial:
                if (!expression_p(obj)->names_enumerate(depth, units))
                    return false;
                break;
            case ID_unit:
                if (algebraic_p value = unit_p(obj)->value())
                    if (symbol_p sym = value->as_quoted<symbol>())
                        if (!names_insert(depth, sym,
                                          units ? unit_p(obj) : nullptr))
                            return false;
                break;
            default:
                break;
            }
        }
        break;
    case ID_polynomial:
    {
        polynomial_g p = polynomial_p(this);
        vars = p->variables();
        for (size_t v = 0; v < vars; v++)
            if (!names_insert(depth, p->variable(v), nullptr))
                return false;
        break;
    }
    default:
        break;
    }
    return true;
}


bool list::names_insert(size_t depth, symbol_p sym, unit_p uobj)
// ----------------------------------------------------------------------------
//   Insert a name according to LNAME sorting order
// ----------------------------------------------------------------------------
{
    size_t existing = rt.depth() - depth;
    size_t len      = 0;
    utf8   name     = sym->value(&len);
    int    cmp      = -1;
    size_t level;
    for (level = 0; level < existing; level++)
    {
        object_p oobj  = rt.stack(level);
        unit_p   ounit = oobj->as<unit>();
        symbol_p other = ounit ? symbol_p(ounit->value()) : symbol_p(oobj);
        size_t   olen  = 0;
        utf8     oname = other->value(&olen);
        cmp = olen - len;
        if (!cmp)
            cmp = symbol::compare(name, oname, olen);
        if (cmp <= 0)
            break;
    }
    if (cmp)
    {
        symbol_g s = sym;
        unit_g   u = uobj;
        if (!rt.push(+sym))
            return false;
        for (size_t i = 0; i < level; i++)
            if (!rt.stack(i, rt.stack(i+1)))
                return false;
        if (!rt.stack(level, uobj ? object_p(+u) : object_p(+s)))
            return false;
    }
    return true;
}


list_p list::names(bool units, id type) const
// ----------------------------------------------------------------------------
//   Return the list of names in an expression
// ----------------------------------------------------------------------------
//   The sorting order defined for `LNAME` in HP50G manual is:
//   - Longer variables first
//   - For same size, alphabetic order
{
    size_t   depth = rt.depth();
    scribble scr;

    if (!names_enumerate(depth, units))
        goto error;

    // Copy the items to the list
    while (rt.depth() > depth)
    {
        object_g obj = rt.pop();
        if (!rt.append(obj))
            goto error;
    }

    return list::make(type, scr.scratch(), scr.growth());

error:
    if (size_t now = rt.depth())
        if (now > depth)
            rt.drop(now - depth);
    return nullptr;
}


static list_p list_variables(object_p obj, object::id type, bool units)
// ----------------------------------------------------------------------------
//   List all variables in an expressions
// ----------------------------------------------------------------------------
{
    object::id objty = obj->type();

    if (objty == object::ID_expression ||
        objty == object::ID_funcall    ||
        objty == object::ID_polynomial ||
        objty == object::ID_program    ||
        objty == object::ID_block      ||
        objty == object::ID_list       ||
        objty == object::ID_array)
        return list_p(obj)->names(units, type);

    rt.type_error();
    return nullptr;
}


COMMAND_BODY(XVars)
// ----------------------------------------------------------------------------
//   List all variable names in an array, leave expression in place
// ----------------------------------------------------------------------------
{
    if (object_p obj = rt.top())
        if (list_p result = list_variables(obj, ID_list, true))
            if (rt.top(result))
                return OK;
    return ERROR;
}


COMMAND_BODY(LName)
// ----------------------------------------------------------------------------
//   List variable names in expression as an array, leave expression in place
// ----------------------------------------------------------------------------
{
    if (object_p obj = rt.top())
        if (list_p result = list_variables(obj, ID_array, false))
            if (rt.push(result))
                return OK;
    return ERROR;
}



// ============================================================================
//
//   DoList and DoSubs commands
//
// ============================================================================

static size_t program_arity(object_p prg)
// ----------------------------------------------------------------------------
//   Compute the program arity for a given program
// ----------------------------------------------------------------------------
{
    if (!prg)
        return -1;
    object::id id = prg->type();
    if (id == object::ID_list || id == object::ID_program)
    {
        if (list_p(prg)->items() != 1)
        {
            rt.argument_count_error();
            return 0;
        }
        prg = list_p(prg)->at(0);
    }
    if (locals_p locs = prg->as<locals>())
        return locs->variables();

    return prg->arity();
}


COMMAND_BODY(DoList)
// ----------------------------------------------------------------------------
//   Apply commands to lists
// ----------------------------------------------------------------------------
{
    if (object_g prg = rt.stack(0))
    {
        if (object_p cntobj = rt.stack(1))
        {
            size_t count = 0;
            size_t base  = 0;
            if (cntobj->is_real())
            {
                count = cntobj->as_uint32(0, true);
                base = 2;
            }
            else if (cntobj->type() == ID_list)
            {
                base = 1;
                count = program_arity(prg);
            }
            else
            {
                rt.type_error();
            }
            if (rt.error())
                return ERROR;
            if (!rt.args(base + count))
                return ERROR;

            // Check consistency of list dimensions
            size_t length = 0;
            id lty = ID_list;
            for (size_t d = 0; d < count; d++)
            {
                object_p obj = rt.stack(base + d);
                if (!obj)
                    return ERROR;
                id ty = obj->type();
                if (ty != ID_list && ty != ID_array)
                    return ERROR;
                list_p lst = list_p(obj);
                if (!d)
                {
                    length = lst->items();
                    lty = ty;
                }
                else if (length != lst->items())
                {
                    rt.dimension_error();
                    return ERROR;
                }
                else if (ty != lty)
                {
                    rt.type_error();
                    return ERROR;
                }
            }

            // Run program on all the lists
            scribble scr;
            size_t depth = rt.depth();
            for (size_t i = 0; i < length; i++)
            {
                size_t offs = base + count - 1;
                for (size_t d = count; d --> 0; )
                {
                    list_p lst = list_p(rt.stack(offs));
                    if (!rt.push(lst->at(i)))
                        return ERROR;
                }

                program::run(prg, true);

                size_t after = rt.depth();
                if (after < depth)
                {
                    rt.missing_argument_error();
                    return ERROR;
                }

                size_t added = after - depth;
                for (size_t d = added; d --> 0; )
                {
                    object_g obj = rt.stack(d);
                    if (!rt.append(obj))
                        return ERROR;
                }
                rt.drop(added);
            }

            list_p result =  list::make(lty, scr.scratch(), scr.growth());
            if (rt.drop(base + count) && rt.push(result))
                return OK;
        }
    }

    return ERROR;
}


static size_t nsub = 0;
static size_t endsub = 0;


COMMAND_BODY(DoSubs)
// ----------------------------------------------------------------------------
//   Apply command to sublist
// ----------------------------------------------------------------------------
{
    if (object_g prg = rt.stack(0))
    {
        if (object_p cntobj = rt.stack(1))
        {
            size_t count = 0;
            size_t base  = 0;
            if (cntobj->is_real())
            {
                count = cntobj->as_uint32(0, true);
                base = 2;
            }
            else if (cntobj->type() == ID_list)
            {
                base = 1;
                count = program_arity(prg);
            }
            else
            {
                rt.type_error();
            }
            if (rt.error())
                return ERROR;
            if (!rt.args(base))
                return ERROR;

            // Check consistency of list dimensions
            object_p obj = rt.stack(base);
            if (!obj)
                return ERROR;
            id lty = obj->type();
            if (lty != ID_list && lty != ID_array)
                return ERROR;
            list_g lst = list_p(obj);
            size_t length = lst->items();
            if (count > length)
            {
                rt.dimension_error();
                return ERROR;
            }

            save<size_t> sns(nsub, 0);
            save<size_t> ses(endsub, length - count + 1);
            rt.drop(base + 1);

            // Run program on all subs in the list
            scribble scr;
            size_t depth = rt.depth();
            for (size_t i = 0; i < endsub; i++)
            {
                nsub = i + 1;
                for (size_t d = 0; d < count; d++)
                    if (!rt.push(lst->at(i + d)))
                        return ERROR;

                program::run(prg, true);

                size_t after = rt.depth();
                if (after < depth)
                {
                    rt.missing_argument_error();
                    return ERROR;
                }

                size_t added = after - depth;
                for (size_t d = added; d --> 0; )
                {
                    object_g obj = rt.stack(d);
                    if (!rt.append(obj))
                        return ERROR;
                }
                rt.drop(added);
            }

            list_p result =  list::make(lty, scr.scratch(), scr.growth());
            if (rt.push(result))
                return OK;
        }
    }

    return ERROR;
}


COMMAND_BODY(NSub)
// ----------------------------------------------------------------------------
//   Return current index in a DoSubs command
// ----------------------------------------------------------------------------
{
    if (!endsub)
        rt.undefined_local_name_error();
    else if (integer_p obj = integer::make(nsub))
        if (rt.push(obj))
            return OK;
    return ERROR;
}


COMMAND_BODY(EndSub)
// ----------------------------------------------------------------------------
//   Return last index in a DoSubs command
// ----------------------------------------------------------------------------
{
    if (!endsub)
        rt.undefined_local_name_error();
    else if (integer_p obj = integer::make(endsub))
        if (rt.push(obj))
            return OK;
    return ERROR;
}
