(* Created with the Wolfram Language : www.wolfram.com *)
{Coupling[gL, {}, 0] -> Coupling[gL, {}, 0], 
 Coupling[gs, {}, 0] -> Coupling[gs, {}, 0], 
 Coupling[gY, {}, 0] -> Coupling[gY, {}, 0], 
 Coupling[Yd, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
  Coupling[cHqd, {Index[i1, Flavor], Index[i2, Flavor]}, 0], 
 Coupling[Ye, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
  Coupling[cHle, {Index[i1, Flavor], Index[i2, Flavor]}, 0], 
 Coupling[Yu, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
  Coupling[cHqu, {Index[i1, Flavor], Index[i2, Flavor]}, 0], 
 Coupling[\[Lambda], {}, 0] -> Coupling[cH4, {}, 0], 
 Coupling[\[Mu]2, {}, 2] -> Coupling[cH2, {}, 2], 
 Coupling[cllHH, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
  (Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]])/
    Coupling[MNR, {Index[d$$1, Flavor]}, 0] - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[gL, {}, 0]^2)/(4*Coupling[MNR, {Index[d$$1, Flavor]}, 0]) - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[gY, {}, 0]^2)/(4*Coupling[MNR, {Index[d$$1, Flavor]}, 0]) + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[i1, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0])/
    (2*Coupling[MNR, {Index[d$$3, Flavor]}, 0]) - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[i1, Flavor]}, 0]]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0])/
    (8*Coupling[MNR, {Index[d$$2, Flavor]}, 0]) - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0])/
    (8*Coupling[MNR, {Index[d$$2, Flavor]}, 0]) - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[i1, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {1, 1}])/
    Coupling[MNR, {Index[d$$3, Flavor]}, 0] - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[i1, Flavor]}, 0]]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {1, 1}])/
    (4*Coupling[MNR, {Index[d$$2, Flavor]}, 0]) - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {1, 1}])/
    (4*Coupling[MNR, {Index[d$$2, Flavor]}, 0]) + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[cH4, {}, 0]*Coupling[MNR, {Index[d$$1, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {1, 2}])/2 + 
   (3*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[gL, {}, 0]^2*Coupling[MNR, {Index[d$$1, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {1, 2}])/8 + 
   (3*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[gY, {}, 0]^2*Coupling[MNR, {Index[d$$1, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {1, 2}])/8, 
 Coupling[cdB, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 0, 
 Coupling[cdd, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 0, 
 Coupling[cdG, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 0, 
 Coupling[cdH, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
  (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[cH4, {}, 0]*Coupling[cHqd, {Index[i1, Flavor], 
       Index[i2, Flavor]}, 0]*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0])/(3*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^
      2) + (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 
       0]]*Coupling[cHqd, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[gL, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0])/(36*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^
      2) - (hbar*Bar[Coupling[cHle, {Index[d$$1, Flavor], 
        Index[d$$2, Flavor]}, 0]]*Bar[Coupling[Yn, {Index[d$$3, Flavor], 
        Index[d$$4, Flavor]}, 0]]*Coupling[cHle, {Index[d$$4, Flavor], 
       Index[d$$2, Flavor]}, 0]*Coupling[cHqd, {Index[i1, Flavor], 
       Index[i2, Flavor]}, 0]*Coupling[Yn, {Index[d$$3, Flavor], 
       Index[d$$1, Flavor]}, 0])/(4*Coupling[MNR, {Index[d$$3, Flavor]}, 0]^
      2) + (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 
       0]]*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[cHqd, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$4, Flavor], Index[d$$3, Flavor]}, 0])/
    (4*Coupling[MNR, {Index[d$$1, Flavor]}, 0]*
     Coupling[MNR, {Index[d$$4, Flavor]}, 0]) - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[cHqd, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[gL, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, 
      {1, 2}])/6 + 
   (hbar*Bar[Coupling[cHle, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
     Coupling[cHle, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]*
     Coupling[cHqd, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$3, Flavor]}, 0]}, {1, 2}])/2 + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
     Coupling[cHqd, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0], Coupling[MNR, 
        {Index[d$$3, Flavor]}, 0]}, {1, 1, 1}])/2 + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[cHqd, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[MNR, {Index[d$$1, Flavor]}, 0]*Coupling[MNR, 
      {Index[d$$4, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$2, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$3, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0], 
       Coupling[MNR, {Index[d$$4, Flavor]}, 0]}, {1, 1, 2}])/2 - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[cHqd, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[MNR, {Index[d$$1, Flavor]}, 0]*Coupling[MNR, 
      {Index[d$$4, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$2, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$3, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0], 
       Coupling[MNR, {Index[d$$4, Flavor]}, 0]}, {2, 1, 1}])/4 + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[cHqd, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[MNR, {Index[d$$1, Flavor]}, 0]*Coupling[MNR, 
      {Index[d$$4, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$2, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$3, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$4, Flavor]}, 0], 
       Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {2, 1, 1}])/4, 
 Coupling[cduq, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 0, 
 Coupling[cduu, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 0, 
 Coupling[cdW, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 0, 
 Coupling[ceB, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
  (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
    Coupling[cHle, {Index[d$$2, Flavor], Index[i2, Flavor]}, 0]*
    Coupling[gY, {}, 0]*Coupling[Yn, {Index[d$$1, Flavor], 
      Index[i1, Flavor]}, 0])/(24*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2), 
 Coupling[ced, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 0, 
 Coupling[cee, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 0, 
 Coupling[ceH, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
  (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[cH4, {}, 0]*Coupling[cHle, {Index[i1, Flavor], 
       Index[i2, Flavor]}, 0]*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0])/(3*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^
      2) + (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 
       0]]*Coupling[cHle, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[gL, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0])/(36*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^
      2) - (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 
       0]]*Coupling[cH4, {}, 0]*Coupling[cHle, {Index[d$$2, Flavor], 
       Index[i2, Flavor]}, 0]*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[i1, Flavor]}, 0])/(4*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^
      2) - (hbar*Bar[Coupling[cHle, {Index[d$$1, Flavor], 
        Index[d$$2, Flavor]}, 0]]*Bar[Coupling[Yn, {Index[d$$3, Flavor], 
        Index[d$$4, Flavor]}, 0]]*Coupling[cHle, {Index[d$$4, Flavor], 
       Index[i2, Flavor]}, 0]*Coupling[cHle, {Index[i1, Flavor], 
       Index[d$$2, Flavor]}, 0]*Coupling[Yn, {Index[d$$3, Flavor], 
       Index[d$$1, Flavor]}, 0])/(3*Coupling[MNR, {Index[d$$3, Flavor]}, 0]^
      2) - (hbar*Bar[Coupling[cHle, {Index[d$$1, Flavor], 
        Index[d$$2, Flavor]}, 0]]*Bar[Coupling[Yn, {Index[d$$3, Flavor], 
        Index[d$$4, Flavor]}, 0]]*Coupling[cHle, {Index[d$$4, Flavor], 
       Index[d$$2, Flavor]}, 0]*Coupling[cHle, {Index[i1, Flavor], 
       Index[i2, Flavor]}, 0]*Coupling[Yn, {Index[d$$3, Flavor], 
       Index[d$$1, Flavor]}, 0])/(4*Coupling[MNR, {Index[d$$3, Flavor]}, 0]^
      2) + (hbar*Bar[Coupling[cHle, {Index[d$$1, Flavor], 
        Index[d$$2, Flavor]}, 0]]*Bar[Coupling[Yn, {Index[d$$3, Flavor], 
        Index[d$$4, Flavor]}, 0]]*Coupling[cHle, {Index[d$$1, Flavor], 
       Index[i2, Flavor]}, 0]*Coupling[cHle, {Index[d$$4, Flavor], 
       Index[d$$2, Flavor]}, 0]*Coupling[Yn, {Index[d$$3, Flavor], 
       Index[i1, Flavor]}, 0])/(8*Coupling[MNR, {Index[d$$3, Flavor]}, 0]^
      2) + (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 
       0]]*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[cHle, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$4, Flavor], Index[d$$3, Flavor]}, 0])/
    (4*Coupling[MNR, {Index[d$$1, Flavor]}, 0]*
     Coupling[MNR, {Index[d$$4, Flavor]}, 0]) - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[cHle, {Index[d$$2, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$4, Flavor], Index[d$$3, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$4, Flavor], Index[i1, Flavor]}, 0])/
    (8*Coupling[MNR, {Index[d$$1, Flavor]}, 0]*
     Coupling[MNR, {Index[d$$4, Flavor]}, 0]) - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[cHle, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[gL, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, 
      {1, 2}])/6 - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[cH4, {}, 0]*Coupling[cHle, {Index[d$$2, Flavor], 
       Index[i2, Flavor]}, 0]*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[i1, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, 
      {1, 2}])/2 + 
   (hbar*Bar[Coupling[cHle, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
     Coupling[cHle, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]*
     Coupling[cHle, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$3, Flavor]}, 0]}, {1, 2}])/2 + 
   (hbar*Bar[Coupling[cHle, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
     Coupling[cHle, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[cHle, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$3, Flavor], Index[i1, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$3, Flavor]}, 0]}, {1, 2}])/4 + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
     Coupling[cHle, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0], Coupling[MNR, 
        {Index[d$$3, Flavor]}, 0]}, {1, 1, 1}])/2 - 
   (3*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
     Coupling[cHle, {Index[d$$2, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$3, Flavor], Index[i1, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0], Coupling[MNR, 
        {Index[d$$3, Flavor]}, 0]}, {1, 1, 1}])/4 + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
     Coupling[cHle, {Index[d$$4, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0], Coupling[MNR, 
        {Index[d$$3, Flavor]}, 0]}, {2, 1, 0}])/4 - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
     Coupling[cHle, {Index[d$$2, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$3, Flavor], Index[i1, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0], Coupling[MNR, 
        {Index[d$$3, Flavor]}, 0]}, {2, 1, 0}])/4 + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[cHle, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[MNR, {Index[d$$1, Flavor]}, 0]*Coupling[MNR, 
      {Index[d$$4, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$2, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$3, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0], 
       Coupling[MNR, {Index[d$$4, Flavor]}, 0]}, {1, 1, 2}])/2 - 
   (3*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[cHle, {Index[d$$2, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[MNR, {Index[d$$1, Flavor]}, 0]*Coupling[MNR, 
      {Index[d$$4, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$3, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[i1, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0], 
       Coupling[MNR, {Index[d$$4, Flavor]}, 0]}, {1, 1, 2}])/4 - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[cHle, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[MNR, {Index[d$$1, Flavor]}, 0]*Coupling[MNR, 
      {Index[d$$4, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$2, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$3, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0], 
       Coupling[MNR, {Index[d$$4, Flavor]}, 0]}, {2, 1, 1}])/4 - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[cHle, {Index[d$$3, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[MNR, {Index[d$$1, Flavor]}, 0]*Coupling[MNR, 
      {Index[d$$4, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$2, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[i1, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0], 
       Coupling[MNR, {Index[d$$4, Flavor]}, 0]}, {2, 1, 1}])/4 + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[cHle, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[MNR, {Index[d$$1, Flavor]}, 0]*Coupling[MNR, 
      {Index[d$$4, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$2, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$3, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$4, Flavor]}, 0], 
       Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {2, 1, 1}])/4 + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[cHle, {Index[d$$3, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[MNR, {Index[d$$1, Flavor]}, 0]*Coupling[MNR, 
      {Index[d$$4, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$2, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[i1, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$4, Flavor]}, 0], 
       Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {2, 1, 1}])/4, 
 Coupling[ceu, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 0, 
 Coupling[ceW, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
  (5*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
    Coupling[cHle, {Index[d$$2, Flavor], Index[i2, Flavor]}, 0]*
    Coupling[gL, {}, 0]*Coupling[Yn, {Index[d$$1, Flavor], 
      Index[i1, Flavor]}, 0])/(24*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2), 
 Coupling[cG, {}, 0] -> 0, Coupling[cGt, {}, 0] -> 0, 
 Coupling[cH, {}, 0] -> 
  (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[cH4, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0])/(3*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^
      2) + (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 
       0]]*Coupling[cH4, {}, 0]*Coupling[gL, {}, 0]^2*
     Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0])/
    (18*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) - 
   (hbar*Bar[Coupling[cHle, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
     Coupling[cH4, {}, 0]*Coupling[cHle, {Index[d$$4, Flavor], 
       Index[d$$2, Flavor]}, 0]*Coupling[Yn, {Index[d$$3, Flavor], 
       Index[d$$1, Flavor]}, 0])/(2*Coupling[MNR, {Index[d$$3, Flavor]}, 0]^
      2) + (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 
       0]]*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[cH4, {}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$2, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$3, Flavor]}, 0])/(2*Coupling[MNR, {Index[d$$1, Flavor]}, 0]*
     Coupling[MNR, {Index[d$$4, Flavor]}, 0]) - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[cH4, {}, 0]*Coupling[gL, {}, 0]^2*
     Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {1, 2}])/3 + 
   hbar*Bar[Coupling[cHle, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
    Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
    Coupling[cH4, {}, 0]*Coupling[cHle, {Index[d$$4, Flavor], 
      Index[d$$2, Flavor]}, 0]*Coupling[Yn, {Index[d$$3, Flavor], 
      Index[d$$1, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$3, Flavor]}, 0]}, 
     {1, 2}] + hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0]]*
    Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
    Coupling[cH4, {}, 0]*Coupling[Yn, {Index[d$$1, Flavor], 
      Index[d$$4, Flavor]}, 0]*Coupling[Yn, {Index[d$$3, Flavor], 
      Index[d$$2, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0], 
      Coupling[MNR, {Index[d$$3, Flavor]}, 0]}, {1, 1, 1}] + 
   hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
    Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
    Coupling[cH4, {}, 0]*Coupling[MNR, {Index[d$$1, Flavor]}, 0]*
    Coupling[MNR, {Index[d$$4, Flavor]}, 0]*
    Coupling[Yn, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]*
    Coupling[Yn, {Index[d$$4, Flavor], Index[d$$3, Flavor]}, 0]*
    LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0], 
      Coupling[MNR, {Index[d$$4, Flavor]}, 0]}, {1, 1, 2}] - 
   (2*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$5, Flavor], Index[d$$6, Flavor]}, 0]]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$3, Flavor], Index[d$$6, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$5, Flavor], Index[d$$2, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0], Coupling[MNR, 
        {Index[d$$3, Flavor]}, 0], Coupling[MNR, {Index[d$$5, Flavor]}, 0]}, 
      {1, 1, 1, 0}])/3 - 
   2*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
    Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
    Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[d$$5, Flavor]}, 0]]*
    Coupling[MNR, {Index[d$$3, Flavor]}, 0]*Coupling[MNR, 
     {Index[d$$6, Flavor]}, 0]*Coupling[Yn, {Index[d$$1, Flavor], 
      Index[d$$4, Flavor]}, 0]*Coupling[Yn, {Index[d$$6, Flavor], 
      Index[d$$2, Flavor]}, 0]*Coupling[Yn, {Index[d$$6, Flavor], 
      Index[d$$5, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0], 
      Coupling[MNR, {Index[d$$3, Flavor]}, 0], Coupling[MNR, 
       {Index[d$$6, Flavor]}, 0]}, {1, 1, 1, 1}], 
 Coupling[cHB, {}, 0] -> 
  (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
    Coupling[gY, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
      Index[d$$2, Flavor]}, 0])/(24*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^
     2), Coupling[cHBox, {}, 0] -> 
  (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[gL, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0])/(24*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^
      2) + (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 
       0]]*Coupling[gY, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0])/(72*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^
      2) - (hbar*Bar[Coupling[cHle, {Index[d$$1, Flavor], 
        Index[d$$2, Flavor]}, 0]]*Bar[Coupling[Yn, {Index[d$$3, Flavor], 
        Index[d$$4, Flavor]}, 0]]*Coupling[cHle, {Index[d$$4, Flavor], 
       Index[d$$2, Flavor]}, 0]*Coupling[Yn, {Index[d$$3, Flavor], 
       Index[d$$1, Flavor]}, 0])/(4*Coupling[MNR, {Index[d$$3, Flavor]}, 0]^
      2) - (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 
       0]]*Coupling[gL, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, 
      {1, 2}])/4 - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[gY, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, 
      {1, 2}])/12 + 
   (hbar*Bar[Coupling[cHle, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
     Coupling[cHle, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$3, Flavor]}, 0]}, {1, 2}])/2 + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0], Coupling[MNR, 
        {Index[d$$3, Flavor]}, 0]}, {1, 1, 1}])/2 - 
   hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
    Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
    Coupling[Yn, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]*
    Coupling[Yn, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]*
    LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0], 
      Coupling[MNR, {Index[d$$3, Flavor]}, 0]}, {2, 1, 0}] + 
   hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
    Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
    Coupling[Yn, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]*
    Coupling[Yn, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]*
    LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0], 
      Coupling[MNR, {Index[d$$3, Flavor]}, 0]}, {3, 1, -1}] + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[MNR, {Index[d$$1, Flavor]}, 0]*Coupling[MNR, 
      {Index[d$$4, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$2, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$3, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0], 
       Coupling[MNR, {Index[d$$4, Flavor]}, 0]}, {1, 1, 2}])/2 - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[MNR, {Index[d$$1, Flavor]}, 0]*Coupling[MNR, 
      {Index[d$$4, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$2, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$3, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0], 
       Coupling[MNR, {Index[d$$4, Flavor]}, 0]}, {2, 1, 1}])/4 + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[MNR, {Index[d$$1, Flavor]}, 0]*Coupling[MNR, 
      {Index[d$$4, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$2, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$3, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0], 
       Coupling[MNR, {Index[d$$4, Flavor]}, 0]}, {3, 1, 0}])/2 - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[MNR, {Index[d$$1, Flavor]}, 0]*Coupling[MNR, 
      {Index[d$$4, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$2, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$3, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$4, Flavor]}, 0], 
       Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {2, 1, 1}])/4 + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[MNR, {Index[d$$1, Flavor]}, 0]*Coupling[MNR, 
      {Index[d$$4, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$2, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$3, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$4, Flavor]}, 0], 
       Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {3, 1, 0}])/2, 
 Coupling[cHBt, {}, 0] -> 0, 
 Coupling[cHd, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
  -1/108*(hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 
        0]]*Coupling[gY, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
        Index[d$$2, Flavor]}, 0]*Delta[Index[i1, Flavor], Index[i2, Flavor]])/
     Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2 + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[gY, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0]*Delta[Index[i1, Flavor], Index[i2, Flavor]]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {1, 2}])/18, 
 Coupling[cHD, {}, 0] -> 
  (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[gY, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0])/(18*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^
      2) + (hbar*Bar[Coupling[cHle, {Index[d$$1, Flavor], 
        Index[d$$2, Flavor]}, 0]]*Bar[Coupling[Yn, {Index[d$$3, Flavor], 
        Index[d$$4, Flavor]}, 0]]*Coupling[cHle, {Index[d$$4, Flavor], 
       Index[d$$2, Flavor]}, 0]*Coupling[Yn, {Index[d$$3, Flavor], 
       Index[d$$1, Flavor]}, 0])/(2*Coupling[MNR, {Index[d$$3, Flavor]}, 0]^
      2) - (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 
       0]]*Coupling[gY, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, 
      {1, 2}])/3 - hbar*Bar[Coupling[cHle, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0]]*
    Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
    Coupling[cHle, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]*
    Coupling[Yn, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
    LF[{Coupling[MNR, {Index[d$$3, Flavor]}, 0]}, {1, 2}] + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0], Coupling[MNR, 
        {Index[d$$3, Flavor]}, 0]}, {1, 1, 1}])/2 + 
   hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
    Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
    Coupling[MNR, {Index[d$$1, Flavor]}, 0]*Coupling[MNR, 
     {Index[d$$4, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
      Index[d$$2, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
      Index[d$$3, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0], 
      Coupling[MNR, {Index[d$$4, Flavor]}, 0]}, {1, 1, 2}], 
 Coupling[cHe, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
  (hbar*Bar[Coupling[cHle, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[cHle, {Index[d$$3, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0])/
    (24*Coupling[MNR, {Index[d$$2, Flavor]}, 0]^2) - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[gY, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0]*Delta[Index[i1, Flavor], Index[i2, Flavor]])/
    (36*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[gY, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0]*Delta[Index[i1, Flavor], Index[i2, Flavor]]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {1, 2}])/6 - 
   (hbar*Bar[Coupling[cHle, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[cHle, {Index[d$$3, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$2, Flavor]}, 0]}, {1, 2}])/4, 
 Coupling[cHG, {}, 0] -> 0, Coupling[cHGt, {}, 0] -> 0, 
 Coupling[cHl1, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
  (Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0])/
    (4*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) + 
   (15*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[gL, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[i1, Flavor]}, 0])/(32*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^
      2) + (55*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], 
        Index[i2, Flavor]}, 0]]*Coupling[gY, {}, 0]^2*
     Coupling[Yn, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0])/
    (288*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0])/
    (8*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) - 
   (hbar*Bar[Coupling[cHle, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[cHle, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[i1, Flavor]}, 0])/
    (16*Coupling[MNR, {Index[d$$2, Flavor]}, 0]^2) - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[i1, Flavor]}, 0])/
    (32*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) - 
   (hbar*Bar[Coupling[cHle, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[cHle, {Index[i1, Flavor], Index[d$$2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0])/
    (16*Coupling[MNR, {Index[d$$3, Flavor]}, 0]^2) - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0])/
    (32*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[Yn, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$3, Flavor], Index[i1, Flavor]}, 0])/
    (4*Coupling[MNR, {Index[d$$1, Flavor]}, 0]*
     Coupling[MNR, {Index[d$$3, Flavor]}, 0]) - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[gY, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0]*Delta[Index[i1, Flavor], Index[i2, Flavor]])/
    (72*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) + 
   (9*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[gL, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[i1, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, 
      {1, 2}])/16 + 
   (11*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[gY, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[i1, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, 
      {1, 2}])/48 - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[MNR, {Index[d$$1, Flavor]}, 0]*Coupling[Yn, 
      {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$3, Flavor], Index[i1, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {1, 2}])/
    (8*Coupling[MNR, {Index[d$$3, Flavor]}, 0]) + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[gY, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0]*Delta[Index[i1, Flavor], Index[i2, Flavor]]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {1, 2}])/12 - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$2, Flavor]}, 0]}, {1, 1}])/
    (4*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[i1, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$2, Flavor]}, 0]}, {1, 1}])/
    (16*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) - 
   (hbar*Bar[Coupling[cHle, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[cHle, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[i1, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$2, Flavor]}, 0]}, {1, 2}])/8 - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$3, Flavor]}, 0]}, {1, 1}])/
    (16*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) - 
   (hbar*Bar[Coupling[cHle, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[cHle, {Index[i1, Flavor], Index[d$$2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$3, Flavor]}, 0]}, {1, 2}])/8 - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[MNR, {Index[d$$3, Flavor]}, 0]*Coupling[Yn, 
      {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$3, Flavor], Index[i1, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$3, Flavor]}, 0]}, {1, 2}])/
    (8*Coupling[MNR, {Index[d$$1, Flavor]}, 0]) + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0], Coupling[MNR, 
        {Index[d$$3, Flavor]}, 0]}, {1, 1, 1}])/2 + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[MNR, {Index[d$$1, Flavor]}, 0]*Coupling[MNR, 
      {Index[d$$3, Flavor]}, 0]*Coupling[Yn, {Index[d$$3, Flavor], 
       Index[d$$2, Flavor]}, 0]*Coupling[Yn, {Index[d$$3, Flavor], 
       Index[i1, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0], 
       Coupling[MNR, {Index[d$$3, Flavor]}, 0]}, {1, 1, 2}])/4, 
 Coupling[cHl3, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
  -1/4*(Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
      Coupling[Yn, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0])/
     Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2 + 
   (35*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[gL, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[i1, Flavor]}, 0])/(288*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^
      2) - (5*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 
       0]]*Coupling[gY, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[i1, Flavor]}, 0])/(32*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^
      2) - (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 
       0]]*Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0])/
    (8*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) - 
   (hbar*Bar[Coupling[cHle, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[cHle, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[i1, Flavor]}, 0])/
    (16*Coupling[MNR, {Index[d$$2, Flavor]}, 0]^2) + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[i1, Flavor]}, 0])/
    (32*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) - 
   (hbar*Bar[Coupling[cHle, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[cHle, {Index[i1, Flavor], Index[d$$2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0])/
    (16*Coupling[MNR, {Index[d$$3, Flavor]}, 0]^2) + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0])/
    (32*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[Yn, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$3, Flavor], Index[i1, Flavor]}, 0])/
    (8*Coupling[MNR, {Index[d$$1, Flavor]}, 0]*
     Coupling[MNR, {Index[d$$3, Flavor]}, 0]) + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[gL, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0]*Delta[Index[i1, Flavor], Index[i2, Flavor]])/
    (72*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) + 
   (7*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[gL, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[i1, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, 
      {1, 2}])/48 - 
   (3*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[gY, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[i1, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, 
      {1, 2}])/16 + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[MNR, {Index[d$$1, Flavor]}, 0]*Coupling[Yn, 
      {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$3, Flavor], Index[i1, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {1, 2}])/
    (8*Coupling[MNR, {Index[d$$3, Flavor]}, 0]) - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[gL, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0]*Delta[Index[i1, Flavor], Index[i2, Flavor]]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {1, 2}])/12 + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$2, Flavor]}, 0]}, {1, 1}])/
    (4*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[i1, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$2, Flavor]}, 0]}, {1, 1}])/
    (16*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) - 
   (hbar*Bar[Coupling[cHle, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[cHle, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[i1, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$2, Flavor]}, 0]}, {1, 2}])/8 + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$3, Flavor]}, 0]}, {1, 1}])/
    (16*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) - 
   (hbar*Bar[Coupling[cHle, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[cHle, {Index[i1, Flavor], Index[d$$2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$3, Flavor]}, 0]}, {1, 2}])/8 + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[MNR, {Index[d$$3, Flavor]}, 0]*Coupling[Yn, 
      {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$3, Flavor], Index[i1, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$3, Flavor]}, 0]}, {1, 2}])/
    (8*Coupling[MNR, {Index[d$$1, Flavor]}, 0]), 
 Coupling[cHq1, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
  (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[gY, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0]*Delta[Index[i1, Flavor], Index[i2, Flavor]])/
    (216*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[gY, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0]*Delta[Index[i1, Flavor], Index[i2, Flavor]]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {1, 2}])/36, 
 Coupling[cHq3, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
  (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[gL, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0]*Delta[Index[i1, Flavor], Index[i2, Flavor]])/
    (72*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[gL, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0]*Delta[Index[i1, Flavor], Index[i2, Flavor]]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {1, 2}])/12, 
 Coupling[cHu, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
  (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[gY, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0]*Delta[Index[i1, Flavor], Index[i2, Flavor]])/
    (54*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[gY, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0]*Delta[Index[i1, Flavor], Index[i2, Flavor]]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {1, 2}])/9, 
 Coupling[cHud, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 0, 
 Coupling[cHW, {}, 0] -> 
  (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
    Coupling[gL, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
      Index[d$$2, Flavor]}, 0])/(24*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^
     2), Coupling[cHWB, {}, 0] -> 
  (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
    Coupling[gL, {}, 0]*Coupling[gY, {}, 0]*
    Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0])/
   (12*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2), 
 Coupling[cHWt, {}, 0] -> 0, Coupling[cHWtB, {}, 0] -> 0, 
 Coupling[cld, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 
  (hbar*Bar[Coupling[cHqd, {Index[d$$1, Flavor], Index[i3, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[cHqd, {Index[d$$1, Flavor], Index[i4, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[i1, Flavor]}, 0])/
    (8*Coupling[MNR, {Index[d$$2, Flavor]}, 0]^2) - 
   (5*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[gY, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[i1, Flavor]}, 0]*Delta[Index[i3, Flavor], Index[i4, Flavor]])/
    (216*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[gY, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[i1, Flavor]}, 0]*Delta[Index[i3, Flavor], Index[i4, Flavor]]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {1, 2}])/36 + 
   (hbar*Bar[Coupling[cHqd, {Index[d$$1, Flavor], Index[i3, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[cHqd, {Index[d$$1, Flavor], Index[i4, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[i1, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$2, Flavor]}, 0]}, {1, 2}])/4, 
 Coupling[cle, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 
  -1/6*(hbar*Bar[Coupling[cHle, {Index[i2, Flavor], Index[i3, Flavor]}, 0]]*
      Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
      Coupling[cHle, {Index[i1, Flavor], Index[i4, Flavor]}, 0]*
      Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0])/
     Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2 + 
   (hbar*Bar[Coupling[cHle, {Index[i2, Flavor], Index[i3, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[cHle, {Index[d$$2, Flavor], Index[i4, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0])/
    (4*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) + 
   (hbar*Bar[Coupling[cHle, {Index[d$$1, Flavor], Index[i3, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[cHle, {Index[i1, Flavor], Index[i4, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[d$$1, Flavor]}, 0])/
    (4*Coupling[MNR, {Index[d$$2, Flavor]}, 0]^2) + 
   (hbar*Bar[Coupling[cHle, {Index[d$$1, Flavor], Index[i3, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[cHle, {Index[d$$1, Flavor], Index[i4, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[i1, Flavor]}, 0])/
    (8*Coupling[MNR, {Index[d$$2, Flavor]}, 0]^2) - 
   (5*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[gY, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[i1, Flavor]}, 0]*Delta[Index[i3, Flavor], Index[i4, Flavor]])/
    (72*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[gY, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[i1, Flavor]}, 0]*Delta[Index[i3, Flavor], Index[i4, Flavor]]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {1, 2}])/12 + 
   (hbar*Bar[Coupling[cHle, {Index[d$$1, Flavor], Index[i3, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[cHle, {Index[d$$1, Flavor], Index[i4, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[i1, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$2, Flavor]}, 0]}, {1, 2}])/4, 
 Coupling[cledq, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 
  (hbar*Bar[Coupling[cHqd, {Index[i4, Flavor], Index[i3, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[cHle, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0])/
    (3*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) - 
   (hbar*Bar[Coupling[cHqd, {Index[i4, Flavor], Index[i3, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[cHle, {Index[d$$2, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0])/
    (2*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2), 
 Coupling[clequ1, {Index[i1_, Flavor], Index[i2_, Flavor], 
    Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 
  -1/3*(hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
      Coupling[cHle, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
      Coupling[cHqu, {Index[i3, Flavor], Index[i4, Flavor]}, 0]*
      Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0])/
     Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2 + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[cHle, {Index[d$$2, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[cHqu, {Index[i3, Flavor], Index[i4, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0])/
    (2*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2), 
 Coupling[clequ3, {Index[i1_, Flavor], Index[i2_, Flavor], 
    Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 0, 
 Coupling[cll, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 
  -1/16*(hbar*Bar[Coupling[cHle, {Index[i4, Flavor], Index[d$$1, Flavor]}, 
        0]]*Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[i2, Flavor]}, 0]]*
      Coupling[cHle, {Index[i3, Flavor], Index[d$$1, Flavor]}, 0]*
      Coupling[Yn, {Index[d$$2, Flavor], Index[i1, Flavor]}, 0])/
     Coupling[MNR, {Index[d$$2, Flavor]}, 0]^2 + 
   (hbar*Bar[Coupling[cHle, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[i4, Flavor]}, 0]]*
     Coupling[cHle, {Index[i3, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[i1, Flavor]}, 0])/
    (16*Coupling[MNR, {Index[d$$2, Flavor]}, 0]^2) + 
   (hbar*Bar[Coupling[cHle, {Index[i4, Flavor], Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[cHle, {Index[i1, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[i3, Flavor]}, 0])/
    (16*Coupling[MNR, {Index[d$$2, Flavor]}, 0]^2) - 
   (hbar*Bar[Coupling[cHle, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[i4, Flavor]}, 0]]*
     Coupling[cHle, {Index[i1, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[i3, Flavor]}, 0])/
    (16*Coupling[MNR, {Index[d$$2, Flavor]}, 0]^2) + 
   (5*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i4, Flavor]}, 0]]*
     Coupling[gL, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[i3, Flavor]}, 0]*Delta[Index[i1, Flavor], Index[i2, Flavor]])/
    (288*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) - 
   (5*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i4, Flavor]}, 0]]*
     Coupling[gY, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[i3, Flavor]}, 0]*Delta[Index[i1, Flavor], Index[i2, Flavor]])/
    (288*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) - 
   (5*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[gL, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[i3, Flavor]}, 0]*Delta[Index[i1, Flavor], Index[i4, Flavor]])/
    (144*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) - 
   (5*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i4, Flavor]}, 0]]*
     Coupling[gL, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[i1, Flavor]}, 0]*Delta[Index[i2, Flavor], Index[i3, Flavor]])/
    (144*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) + 
   (5*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[gL, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[i1, Flavor]}, 0]*Delta[Index[i3, Flavor], Index[i4, Flavor]])/
    (288*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) - 
   (5*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[gY, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[i1, Flavor]}, 0]*Delta[Index[i3, Flavor], Index[i4, Flavor]])/
    (288*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i4, Flavor]}, 0]]*
     Coupling[gL, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[i3, Flavor]}, 0]*Delta[Index[i1, Flavor], Index[i2, Flavor]]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {1, 2}])/48 - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i4, Flavor]}, 0]]*
     Coupling[gY, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[i3, Flavor]}, 0]*Delta[Index[i1, Flavor], Index[i2, Flavor]]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {1, 2}])/48 - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[gL, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[i3, Flavor]}, 0]*Delta[Index[i1, Flavor], Index[i4, Flavor]]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {1, 2}])/24 - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i4, Flavor]}, 0]]*
     Coupling[gL, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[i1, Flavor]}, 0]*Delta[Index[i2, Flavor], Index[i3, Flavor]]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {1, 2}])/24 + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[gL, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[i1, Flavor]}, 0]*Delta[Index[i3, Flavor], Index[i4, Flavor]]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {1, 2}])/48 - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[gY, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[i1, Flavor]}, 0]*Delta[Index[i3, Flavor], Index[i4, Flavor]]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {1, 2}])/48 - 
   (hbar*Bar[Coupling[cHle, {Index[i4, Flavor], Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[cHle, {Index[i3, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[i1, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$2, Flavor]}, 0]}, {1, 2}])/8 + 
   (hbar*Bar[Coupling[cHle, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[i4, Flavor]}, 0]]*
     Coupling[cHle, {Index[i3, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[i1, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$2, Flavor]}, 0]}, {1, 2}])/8 + 
   (hbar*Bar[Coupling[cHle, {Index[i4, Flavor], Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[cHle, {Index[i1, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[i3, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$2, Flavor]}, 0]}, {1, 2}])/8 - 
   (hbar*Bar[Coupling[cHle, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[i4, Flavor]}, 0]]*
     Coupling[cHle, {Index[i1, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[i3, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$2, Flavor]}, 0]}, {1, 2}])/8 + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[i4, Flavor]}, 0]]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[i3, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[i1, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0], Coupling[MNR, 
        {Index[d$$2, Flavor]}, 0]}, {1, 1, 1}])/8 + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i4, Flavor]}, 0]]*
     Coupling[MNR, {Index[d$$1, Flavor]}, 0]*Coupling[MNR, 
      {Index[d$$2, Flavor]}, 0]*Coupling[Yn, {Index[d$$2, Flavor], 
       Index[i1, Flavor]}, 0]*Coupling[Yn, {Index[d$$2, Flavor], 
       Index[i3, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0], 
       Coupling[MNR, {Index[d$$2, Flavor]}, 0]}, {1, 1, 2}])/4, 
 Coupling[clq1, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 
  -1/16*(hbar*Bar[Coupling[cHqd, {Index[i4, Flavor], Index[d$$1, Flavor]}, 
        0]]*Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[i2, Flavor]}, 0]]*
      Coupling[cHqd, {Index[i3, Flavor], Index[d$$1, Flavor]}, 0]*
      Coupling[Yn, {Index[d$$2, Flavor], Index[i1, Flavor]}, 0])/
     Coupling[MNR, {Index[d$$2, Flavor]}, 0]^2 + 
   (hbar*Bar[Coupling[cHqu, {Index[i4, Flavor], Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[cHqu, {Index[i3, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[i1, Flavor]}, 0])/
    (16*Coupling[MNR, {Index[d$$2, Flavor]}, 0]^2) + 
   (5*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[gY, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[i1, Flavor]}, 0]*Delta[Index[i3, Flavor], Index[i4, Flavor]])/
    (432*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[gY, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[i1, Flavor]}, 0]*Delta[Index[i3, Flavor], Index[i4, Flavor]]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {1, 2}])/72 - 
   (hbar*Bar[Coupling[cHqd, {Index[i4, Flavor], Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[cHqd, {Index[i3, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[i1, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$2, Flavor]}, 0]}, {1, 2}])/8 + 
   (hbar*Bar[Coupling[cHqu, {Index[i4, Flavor], Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[cHqu, {Index[i3, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[i1, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$2, Flavor]}, 0]}, {1, 2}])/8, 
 Coupling[clq3, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 
  (hbar*Bar[Coupling[cHqd, {Index[i4, Flavor], Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[cHqd, {Index[i3, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[i1, Flavor]}, 0])/
    (16*Coupling[MNR, {Index[d$$2, Flavor]}, 0]^2) + 
   (hbar*Bar[Coupling[cHqu, {Index[i4, Flavor], Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[cHqu, {Index[i3, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[i1, Flavor]}, 0])/
    (16*Coupling[MNR, {Index[d$$2, Flavor]}, 0]^2) - 
   (5*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[gL, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[i1, Flavor]}, 0]*Delta[Index[i3, Flavor], Index[i4, Flavor]])/
    (144*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[gL, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[i1, Flavor]}, 0]*Delta[Index[i3, Flavor], Index[i4, Flavor]]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {1, 2}])/24 + 
   (hbar*Bar[Coupling[cHqd, {Index[i4, Flavor], Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[cHqd, {Index[i3, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[i1, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$2, Flavor]}, 0]}, {1, 2}])/8 + 
   (hbar*Bar[Coupling[cHqu, {Index[i4, Flavor], Index[d$$1, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[cHqu, {Index[i3, Flavor], Index[d$$1, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[i1, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$2, Flavor]}, 0]}, {1, 2}])/8, 
 Coupling[clu, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 
  -1/8*(hbar*Bar[Coupling[cHqu, {Index[d$$1, Flavor], Index[i3, Flavor]}, 0]]*
      Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[i2, Flavor]}, 0]]*
      Coupling[cHqu, {Index[d$$1, Flavor], Index[i4, Flavor]}, 0]*
      Coupling[Yn, {Index[d$$2, Flavor], Index[i1, Flavor]}, 0])/
     Coupling[MNR, {Index[d$$2, Flavor]}, 0]^2 + 
   (5*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[gY, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[i1, Flavor]}, 0]*Delta[Index[i3, Flavor], Index[i4, Flavor]])/
    (108*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[gY, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[i1, Flavor]}, 0]*Delta[Index[i3, Flavor], Index[i4, Flavor]]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {1, 2}])/18 - 
   (hbar*Bar[Coupling[cHqu, {Index[d$$1, Flavor], Index[i3, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$2, Flavor], Index[i2, Flavor]}, 0]]*
     Coupling[cHqu, {Index[d$$1, Flavor], Index[i4, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$2, Flavor], Index[i1, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$2, Flavor]}, 0]}, {1, 2}])/4, 
 Coupling[cqd1, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 
  -1/18*(hbar*Bar[Coupling[cHqd, {Index[i2, Flavor], Index[i3, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[cHqd, {Index[i1, Flavor], Index[i4, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0])/
    Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2, 
 Coupling[cqd8, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 
  -1/3*(hbar*Bar[Coupling[cHqd, {Index[i2, Flavor], Index[i3, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[cHqd, {Index[i1, Flavor], Index[i4, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0])/
    Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2, 
 Coupling[cqe, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 0, 
 Coupling[cqq1, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 0, 
 Coupling[cqq3, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 0, 
 Coupling[cqqq, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 0, 
 Coupling[cqqu, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 0, 
 Coupling[cqu1, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 
  -1/18*(hbar*Bar[Coupling[cHqu, {Index[i2, Flavor], Index[i3, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[cHqu, {Index[i1, Flavor], Index[i4, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0])/
    Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2, 
 Coupling[cqu8, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 
  -1/3*(hbar*Bar[Coupling[cHqu, {Index[i2, Flavor], Index[i3, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[cHqu, {Index[i1, Flavor], Index[i4, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0])/
    Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2, 
 Coupling[cquqd1, {Index[i1_, Flavor], Index[i2_, Flavor], 
    Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 
  (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
    Coupling[cHqd, {Index[i3, Flavor], Index[i4, Flavor]}, 0]*
    Coupling[cHqu, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
    Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0])/
   (3*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2), 
 Coupling[cquqd8, {Index[i1_, Flavor], Index[i2_, Flavor], 
    Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 0, 
 Coupling[cuB, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 0, 
 Coupling[cud1, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 0, 
 Coupling[cud8, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 0, 
 Coupling[cuG, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 0, 
 Coupling[cuH, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
  (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[cH4, {}, 0]*Coupling[cHqu, {Index[i1, Flavor], 
       Index[i2, Flavor]}, 0]*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0])/(3*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^
      2) + (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 
       0]]*Coupling[cHqu, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[gL, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0])/(36*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^
      2) - (hbar*Bar[Coupling[cHle, {Index[d$$1, Flavor], 
        Index[d$$2, Flavor]}, 0]]*Bar[Coupling[Yn, {Index[d$$3, Flavor], 
        Index[d$$4, Flavor]}, 0]]*Coupling[cHle, {Index[d$$4, Flavor], 
       Index[d$$2, Flavor]}, 0]*Coupling[cHqu, {Index[i1, Flavor], 
       Index[i2, Flavor]}, 0]*Coupling[Yn, {Index[d$$3, Flavor], 
       Index[d$$1, Flavor]}, 0])/(4*Coupling[MNR, {Index[d$$3, Flavor]}, 0]^
      2) + (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 
       0]]*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[cHqu, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$4, Flavor], Index[d$$3, Flavor]}, 0])/
    (4*Coupling[MNR, {Index[d$$1, Flavor]}, 0]*
     Coupling[MNR, {Index[d$$4, Flavor]}, 0]) - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[cHqu, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[gL, {}, 0]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, 
      {1, 2}])/6 + 
   (hbar*Bar[Coupling[cHle, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
     Coupling[cHle, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]*
     Coupling[cHqu, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$3, Flavor], Index[d$$1, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$3, Flavor]}, 0]}, {1, 2}])/2 + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
     Coupling[cHqu, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0], Coupling[MNR, 
        {Index[d$$3, Flavor]}, 0]}, {1, 1, 1}])/2 + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[cHqu, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[MNR, {Index[d$$1, Flavor]}, 0]*Coupling[MNR, 
      {Index[d$$4, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$2, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$3, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0], 
       Coupling[MNR, {Index[d$$4, Flavor]}, 0]}, {1, 1, 2}])/2 + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[cHqu, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[MNR, {Index[d$$1, Flavor]}, 0]*Coupling[MNR, 
      {Index[d$$4, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$2, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$3, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0], 
       Coupling[MNR, {Index[d$$4, Flavor]}, 0]}, {2, 1, 1}])/4 - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[cHqu, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[MNR, {Index[d$$1, Flavor]}, 0]*Coupling[MNR, 
      {Index[d$$4, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$2, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$3, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$4, Flavor]}, 0], 
       Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {2, 1, 1}])/4, 
 Coupling[cuu, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 0, 
 Coupling[cuW, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 0, 
 Coupling[cW, {}, 0] -> 0, Coupling[cWt, {}, 0] -> 0, 
 Coupling[cHqd, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
  Coupling[Yd, {Index[i1, Flavor], Index[i2, Flavor]}, 0] + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[Yd, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0])/4 + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[CH2, {}, 2]*Coupling[Yd, {Index[i1, Flavor], 
       Index[i2, Flavor]}, 0]*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0])/(3*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^
      2) - (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 
       0]]*Coupling[Yd, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {1, 1}])/2, 
 Coupling[cHqu, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
  Coupling[Yu, {Index[i1, Flavor], Index[i2, Flavor]}, 0] + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]*
     Coupling[Yu, {Index[i1, Flavor], Index[i2, Flavor]}, 0])/4 + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[CH2, {}, 2]*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0]*Coupling[Yu, {Index[i1, Flavor], 
       Index[i2, Flavor]}, 0])/(3*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^
      2) - (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 
       0]]*Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]*
     Coupling[Yu, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {1, 1}])/2, 
 Coupling[cHle, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
  Coupling[Ye, {Index[i1, Flavor], Index[i2, Flavor]}, 0] + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[Ye, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0])/4 + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[CH2, {}, 2]*Coupling[Ye, {Index[i1, Flavor], 
       Index[i2, Flavor]}, 0]*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0])/(3*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^
      2) - (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 
       0]]*Coupling[Ye, {Index[d$$2, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0])/8 - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[CH2, {}, 2]*Coupling[Ye, {Index[d$$2, Flavor], 
       Index[i2, Flavor]}, 0]*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[i1, Flavor]}, 0])/(4*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^
      2) - (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 
       0]]*Coupling[Ye, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {1, 1}])/2 + 
   (3*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[Ye, {Index[d$$2, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[Yn, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {1, 1}])/4 - 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[CH2, {}, 2]*Coupling[Ye, {Index[d$$2, Flavor], 
       Index[i2, Flavor]}, 0]*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[i1, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, 
      {1, 2}])/2, Coupling[cH2, {}, 2] -> Coupling[CH2, {}, 2] + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[CH2, {}, 2]*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0])/2 + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[CH2, {}, 2]^2*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0])/(3*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^
      2) - hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 
      0]]*Coupling[CH2, {}, 2]*Coupling[Yn, {Index[d$$1, Flavor], 
      Index[d$$2, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, 
     {1, 1}], Coupling[cH4, {}, 0] -> 
  (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[CH2, {}, 2]*Coupling[gL, {}, 0]^2*
     Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0])/
    (9*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) - 
   (hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
     Coupling[CH2, {}, 2]*Coupling[Ye, {Index[d$$4, Flavor], 
       Index[d$$2, Flavor]}, 0]*Coupling[Yn, {Index[d$$3, Flavor], 
       Index[d$$1, Flavor]}, 0])/Coupling[MNR, {Index[d$$3, Flavor]}, 0]^2 + 
   (hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
     Coupling[CH2, {}, 2]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$2, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
       Index[d$$3, Flavor]}, 0])/(Coupling[MNR, {Index[d$$1, Flavor]}, 0]*
     Coupling[MNR, {Index[d$$4, Flavor]}, 0]) + Coupling[\[Lambda], {}, 0] + 
   hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
    Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]*
    Coupling[\[Lambda], {}, 0] + 
   (4*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[CH2, {}, 2]*Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0]*Coupling[\[Lambda], {}, 0])/
    (3*Coupling[MNR, {Index[d$$1, Flavor]}, 0]^2) - 
   2*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
    Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]*
    Coupling[\[Lambda], {}, 0]*LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, 
     {1, 1}] - (2*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], 
        Index[d$$2, Flavor]}, 0]]*Coupling[CH2, {}, 2]*Coupling[gL, {}, 0]^2*
     Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]*
     LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {1, 2}])/3 + 
   2*hbar*Bar[Coupling[Ye, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
    Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
    Coupling[CH2, {}, 2]*Coupling[Ye, {Index[d$$4, Flavor], 
      Index[d$$2, Flavor]}, 0]*Coupling[Yn, {Index[d$$3, Flavor], 
      Index[d$$1, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$3, Flavor]}, 0]}, 
     {1, 2}] + 2*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], 
       Index[d$$2, Flavor]}, 0]]*
    Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
    Coupling[Yn, {Index[d$$1, Flavor], Index[d$$4, Flavor]}, 0]*
    Coupling[Yn, {Index[d$$3, Flavor], Index[d$$2, Flavor]}, 0]*
    LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0], 
      Coupling[MNR, {Index[d$$3, Flavor]}, 0]}, {1, 1, 0}] + 
   2*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
    Bar[Coupling[Yn, {Index[d$$3, Flavor], Index[d$$4, Flavor]}, 0]]*
    Coupling[CH2, {}, 2]*Coupling[Yn, {Index[d$$1, Flavor], 
      Index[d$$4, Flavor]}, 0]*Coupling[Yn, {Index[d$$3, Flavor], 
      Index[d$$2, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0], 
      Coupling[MNR, {Index[d$$3, Flavor]}, 0]}, {1, 1, 1}] + 
   2*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
    Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
    Coupling[MNR, {Index[d$$1, Flavor]}, 0]*Coupling[MNR, 
     {Index[d$$4, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
      Index[d$$2, Flavor]}, 0]*Coupling[Yn, {Index[d$$4, Flavor], 
      Index[d$$3, Flavor]}, 0]*LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0], 
      Coupling[MNR, {Index[d$$4, Flavor]}, 0]}, {1, 1, 1}] + 
   2*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
    Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$3, Flavor]}, 0]]*
    Coupling[CH2, {}, 2]*Coupling[MNR, {Index[d$$1, Flavor]}, 0]*
    Coupling[MNR, {Index[d$$4, Flavor]}, 0]*
    Coupling[Yn, {Index[d$$4, Flavor], Index[d$$2, Flavor]}, 0]*
    Coupling[Yn, {Index[d$$4, Flavor], Index[d$$3, Flavor]}, 0]*
    LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0], 
      Coupling[MNR, {Index[d$$4, Flavor]}, 0]}, {1, 1, 2}], 
 Coupling[CH2, {}, 2] -> Coupling[\[Mu]2, {}, 2] - 
   2*hbar*Bar[Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
    Coupling[Yn, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]*
    LF[{Coupling[MNR, {Index[d$$1, Flavor]}, 0]}, {1, 0}]}
