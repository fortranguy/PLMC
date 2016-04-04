title: Environment

#FloorPenetration and WallsPotential
[[class_floor_penetration:Abstract_meet]] tells whether a position is above the floor
or below (i.e. overlap). It also gives the shortest vector from floor. Such vector will be useful
to calculate the potential between a particle and the wall.
[[class_walls_potential:Abstract_visit]] visits 2 [[Abstract_Floor_Penetration]] :
one below (floor) and one above (ceiling) which are mirror symmetric with respect to
the \( (x, y, 0) \) plan. The name may be confusing: walls are horizontal not vertical.

#CenteredBlockPenetration
[[Centered_Block_Penetration]] is a flat floor with a rounded block at the center, cf.
modules/environment/centered_block_penetration.tex which shows the right half.
When using [[class_floor_penetration:Block_meet]], if a position is in a blue area,
shortestVectorFromFloor's origin will be on a rounder corner. Otherwise (i.e. white area),
it will be on a flat portion.
