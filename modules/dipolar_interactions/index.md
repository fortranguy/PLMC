title: Dipolar interactions

#Method In Bulk
Dipolar Ewald Summation (DES) is implemented to handle dipolar interactions in 3D periodicity.

#Slab Correction
In 2D periodicity, Dipolar Layer Correction (DLC) is implemented.
In [[classes_dlc_visitor:Abstract_visit_move]], [[classes_dlc_visitor:Abstract_visit_rotation]],
[[classes_dlc_visitor:Abstract_visit_exchange]] and [[classes_dlc_visitor:Abstract_visit_switch]],
the same pattern is used for the difference of energy:
\[
    \Delta U = \sum_{\vec{k}_{1:2}} w(\vec{k}_{1:2}) \Re\left[
        S_+(\vec{k}_{1:2}) \Delta S_-^\ast(\vec{k}_{1:2}) +
        S_-^\ast(\vec{k}_{1:2}) \Delta S_+(\vec{k}_{1:2}) +
        \Delta S_+(\vec{k}_{1:2}) \Delta S_-^\ast(\vec{k}_{1:2})
    \right].
\]
