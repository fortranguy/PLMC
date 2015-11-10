title: Bridging Physics and Computer Science

#Design

This program aims at being fast and versatile using Modern
Fortran features.

#Base units

The most general situation (i.e. short-range and long-range
interactions) requires 3 bases units, e.g.:

\( u_l \): unit of length

\( u_E \): unit of energy

\( u_\epsilon \): unit of permittivity.

If long-range interactions aren't used, the first 2 will
suffice.

#Constructors and Destructors

Module procedures and Final subroutines may improve
the implementations of constructors and destructors.
However, current compilers may not handle them correctly.

#Null Objects

Null Objects are used to deactivate a feature without
impairing the readability nor the stability of the code.
However, they currently raise warnings for
``unused dummy arguments''.
A new intent may solve this problem (e.g. [intent(none)](http://fortranwiki.org/fortran/show/INTENT%28NONE%29)).
