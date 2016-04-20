title: Environment

#FloorPenetration and WallsPotential
[[classes_floor_penetration:Abstract_meet]] tells whether a position is above the floor
or below (i.e. overlap). It also gives the shortest vector from floor. Such vector will be useful
to calculate the potential between a particle and the wall.
[[classes_walls_potential:Abstract_visit]] visits 2 [[Abstract_Floor_Penetration]] :
one below (floor) and one above (ceiling) which are mirror symmetric with respect to
the \( (x, y, 0) \) plan. The name may be confusing: walls are horizontal not vertical.
