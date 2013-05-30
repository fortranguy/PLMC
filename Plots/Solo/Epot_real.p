epot1(alpha, r) = erfc(alpha*r)/r**3 + 2.*alpha/sqrt(pi) * exp(-alpha**2*r**2) / r**2

epot2(alpha, r) = 3.*erfc(alpha*r)/r**5 + 2.*alpha/sqrt(pi) * (2*alpha**2+3./r**2) * exp(-alpha**2*r**2) / r**2
