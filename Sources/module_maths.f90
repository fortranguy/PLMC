module module_maths

implicit none

private
public gcd, lcm

contains

    integer recursive function gcd(a, b) result(great_common_divisor)
        
        integer, intent(in) :: a, b
        
        if (a == b) great_common_divisor = a
        if (a > b) great_common_divisor = gcd(a - b, b)
        if (a < b) great_common_divisor = gcd(a, b - a)
    
    end function gcd
    
    integer function lcm(a, b)
    
        integer, intent(in) :: a, b
        
        lcm = a * b / gcd(a, b)
    
    end function lcm

end module module_maths
