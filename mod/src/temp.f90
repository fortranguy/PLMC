GC
xNew(1:2) = Lsize(1:2) * xRand(1:2)
        xNew(3) = (Height - this%get_sigma()) * xRand(3) + this%get_sigma()/2._DP
        this%Volume (GC)
        
class HS:
data_name = "Particles.Hard Spheres.number of Widom particles"
        call data_json%get(data_name, this%widom_num_particles, found)
        call test_data_found(data_name, found)
