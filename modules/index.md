title: Bridging Physics and Computer Science

#Design

This program aims at being readable, fast and versatile using
Modern Fortran features.

#Units

The most general situation (i.e. short-range and dipolar interactions) requires
2 fundamental constants and 2 bases units, e.g.:

- \( k_\text{B} \): Boltzmann constant
- \( \varepsilon \): permittivity
- \( u_l \): unit of length
- \( u_E \): unit of energy.
If dipolar interactions aren't used, \( \varepsilon \) is useless.

#Constructors and Destructors

Module procedures and Final subroutines may improve
the implementations of constructors and destructors.
However, current compilers may not handle them correctly.

#Null Objects

Null Objects are used to deactivate a class behaviour without
impairing the readability nor the stability of the code.
However, they currently raise warnings for
``unused dummy arguments''.
A new intent may solve this problem (e.g.
    [intent(none)](http://fortranwiki.org/fortran/show/INTENT%28NONE%29)).

#References
- [SMAC](http://www.lps.ens.fr/~krauth/index.php/SMAC)
- [MolSim](https://www.elsevier.com/books/understanding-molecular-simulation/frenkel/978-0-12-267351-1)
