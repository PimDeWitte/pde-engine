#!/usr/bin/env python3
"""
A set of mathematical operations for use in expression generation.

This module defines standard unary, binary, and special operations
as top-level functions to ensure they are pickleable for multiprocessing.
"""

from sympy import sqrt, exp, log, Rational

def op_add(x, y):
    """Adds two expressions."""
    return x + y

def op_sub(x, y):
    """Subtracts two expressions."""
    return x - y

def op_mul(x, y):
    """Multiplies two expressions."""
    return x * y

def op_div(x, y):
    """Divides two expressions."""
    return x / y

def op_geom_sum(x, y):
    """Computes the geometric sum x / (1 - y)."""
    return x / (1 - y)

def op_neg(x):
    """Negates an expression."""
    return -x

def op_inv(x):
    """Computes the inverse of an expression."""
    return 1 / x

def op_sqrt(x):
    """Computes the square root of an expression."""
    return sqrt(x)

def op_square(x):
    """Squares an expression."""
    return x**2

def op_pow_3_2(x):
    """Raises an expression to the power of 3/2."""
    return x**Rational(3, 2)

def op_pow_neg_3_2(x):
    """Raises an expression to the power of -3/2."""
    return x**Rational(-3, 2)

def op_exp(x):
    """Computes the exponential of an expression."""
    return exp(x)

def op_exp_neg(x):
    """Computes the exponential of a negated expression."""
    return exp(-x)

def op_sqrt_shift_neg(x, y):
    """Computes sqrt((x - 1)**2 + y**2)."""
    return sqrt((x - 1)**2 + y**2)

def op_sqrt_shift_pos(x, y):
    """Computes sqrt((x + 1)**2 + y**2)."""
    return sqrt((x + 1)**2 + y**2)

def op_exp_mul(x, y):
    """Computes x * exp(y)."""
    return x * exp(y)

def op_log_mul(x, y):
    """Computes x * log(y)."""
    return x * log(y)

# Operation dictionaries for easy import and use
UNARY_OPS = {
    'neg': op_neg,
    'inv': op_inv,
    'sqrt': op_sqrt,
    'square': op_square,
    'pow_3_2': op_pow_3_2,
    'pow_neg_3_2': op_pow_neg_3_2,
    'exp': op_exp,
    'exp_neg': op_exp_neg,
}

BINARY_OPS = {
    'add': op_add,
    'sub': op_sub,
    'mul': op_mul,
    'div': op_div,
    'geom_sum': op_geom_sum,
}

SPECIAL_OPS = {
    'sqrt_shift_neg': op_sqrt_shift_neg,
    'sqrt_shift_pos': op_sqrt_shift_pos,
    'exp_mul': op_exp_mul,
    'log_mul': op_log_mul,
}

ALL_BINARY_OPS = {**BINARY_OPS, **SPECIAL_OPS}
