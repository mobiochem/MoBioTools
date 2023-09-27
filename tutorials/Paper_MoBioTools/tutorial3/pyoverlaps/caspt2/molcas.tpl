&header
title = p-AZO molecule
group = C1
basis set = cc-pVDZ spherical all
RICD
&end

&seward
&end

&chgspin
0,1,0,1,0,1
&end

&externchg
&end

&guess
jobiph
&end

&caspt2
IPEA= 0.0
SHIFT= 0.0
IMAGinary= 0.2
maxiter=100
multistate = 10 1 2 3 4 5 6 7 8 9 10
PROP
&end

&rassi
Nr of JobIph = 1 10; 1 2 3 4 5 6 7 8 9 10
HEFF
&end
