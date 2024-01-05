function circle(x)
    rr = 0.5*1.0 #chord 1.0
    return sqrt(rr^2-(x-rr)^2)
end

function NACA00(x; t=0.12)
    return 5*t*(0.2969*x^0.5-0.1260*x-0.3516*x^2+0.2843*x^3-0.1015*x^4)
end


"""
    CST_NACA0012()
It gives the CST of the NACA0012 airfoil
"""
function CST_NACA0012()
wu = [ 0.17136426990915907,0.13943231948894877,0.14906429056622653]
wl = -1 .* copy(wu)
return CST(wu,wl)
end
