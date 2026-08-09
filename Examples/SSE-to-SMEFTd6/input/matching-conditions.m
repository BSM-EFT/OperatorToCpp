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
 Coupling[\[Lambda], {}, 0] -> -2*Coupling[cH4, {}, 0], 
 Coupling[\[Mu]2, {}, 2] -> Coupling[cH2, {}, 2], 
 Coupling[cllHH, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 0, 
 Coupling[cdB, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 0, 
 Coupling[cdd, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 0, 
 Coupling[cdG, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 0, 
 Coupling[cdH, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
  (5*hbar*Bar[Coupling[cHqd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[cHqd, {Index[d$$1, Flavor], 
       Index[i2, Flavor]}, 0]*Coupling[cHqd, {Index[i1, Flavor], 
       Index[d$$2, Flavor]}, 0])/(4*Coupling[M, {}, 0]^4) - 
   (29*hbar*Coupling[A, {}, 0]^2*Coupling[cH4, {}, 0]*
     Coupling[cHqd, {Index[i1, Flavor], Index[i2, Flavor]}, 0])/
    (6*Coupling[M, {}, 0]^4) - (31*hbar*Coupling[A, {}, 0]^2*
     Coupling[cHqd, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[gL, {}, 0]^2)/(36*Coupling[M, {}, 0]^4) - 
   (3*hbar*Coupling[A, {}, 0]^2*Coupling[cHqd, {Index[i1, Flavor], 
       Index[i2, Flavor]}, 0]*Coupling[\[Kappa], {}, 0])/
    (4*Coupling[M, {}, 0]^4) + (hbar*Coupling[A, {}, 0]^3*
     Coupling[cHqd, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[\[Mu], {}, 0])/(4*Coupling[M, {}, 0]^6) + 
   (3*hbar*Bar[Coupling[cHqd, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[cHqd, {Index[d$$1, Flavor], 
       Index[i2, Flavor]}, 0]*Coupling[cHqd, {Index[i1, Flavor], 
       Index[d$$2, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
    (2*Coupling[M, {}, 0]^4) - (hbar*Coupling[A, {}, 0]^2*
     Coupling[cH4, {}, 0]*Coupling[cHqd, {Index[i1, Flavor], 
       Index[i2, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
    Coupling[M, {}, 0]^4 - (5*hbar*Coupling[A, {}, 0]^2*
     Coupling[cHqd, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[gL, {}, 0]^2*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
    (6*Coupling[M, {}, 0]^4), 
 Coupling[cduq, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 0, 
 Coupling[cduu, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 0, 
 Coupling[cdW, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 0, 
 Coupling[ceB, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 0, 
 Coupling[ced, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 0, 
 Coupling[cee, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 0, 
 Coupling[ceH, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
  (5*hbar*Bar[Coupling[cHle, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[cHle, {Index[d$$1, Flavor], 
       Index[i2, Flavor]}, 0]*Coupling[cHle, {Index[i1, Flavor], 
       Index[d$$2, Flavor]}, 0])/(4*Coupling[M, {}, 0]^4) - 
   (29*hbar*Coupling[A, {}, 0]^2*Coupling[cH4, {}, 0]*
     Coupling[cHle, {Index[i1, Flavor], Index[i2, Flavor]}, 0])/
    (6*Coupling[M, {}, 0]^4) - (31*hbar*Coupling[A, {}, 0]^2*
     Coupling[cHle, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[gL, {}, 0]^2)/(36*Coupling[M, {}, 0]^4) - 
   (3*hbar*Coupling[A, {}, 0]^2*Coupling[cHle, {Index[i1, Flavor], 
       Index[i2, Flavor]}, 0]*Coupling[\[Kappa], {}, 0])/
    (4*Coupling[M, {}, 0]^4) + (hbar*Coupling[A, {}, 0]^3*
     Coupling[cHle, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[\[Mu], {}, 0])/(4*Coupling[M, {}, 0]^6) + 
   (3*hbar*Bar[Coupling[cHle, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[cHle, {Index[d$$1, Flavor], 
       Index[i2, Flavor]}, 0]*Coupling[cHle, {Index[i1, Flavor], 
       Index[d$$2, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
    (2*Coupling[M, {}, 0]^4) - (hbar*Coupling[A, {}, 0]^2*
     Coupling[cH4, {}, 0]*Coupling[cHle, {Index[i1, Flavor], 
       Index[i2, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
    Coupling[M, {}, 0]^4 - (5*hbar*Coupling[A, {}, 0]^2*
     Coupling[cHle, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[gL, {}, 0]^2*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
    (6*Coupling[M, {}, 0]^4), 
 Coupling[ceu, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 0, 
 Coupling[ceW, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 0, 
 Coupling[cG, {}, 0] -> 0, Coupling[cGt, {}, 0] -> 0, 
 Coupling[cH, {}, 0] -> (-2*hbar*Coupling[A, {}, 0]^6)/
    (3*Coupling[M, {}, 0]^8) + (12*hbar*Coupling[A, {}, 0]^4*
     Coupling[cH4, {}, 0])/Coupling[M, {}, 0]^6 + 
   (164*hbar*Coupling[A, {}, 0]^2*Coupling[cH4, {}, 0]^2)/
    (3*Coupling[M, {}, 0]^4) + (31*hbar*Coupling[A, {}, 0]^2*
     Coupling[cH4, {}, 0]*Coupling[gL, {}, 0]^2)/(9*Coupling[M, {}, 0]^4) + 
   (29*hbar*Coupling[A, {}, 0]^4*Coupling[\[Kappa], {}, 0])/
    (4*Coupling[M, {}, 0]^6) - 
   (Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0])/
    (2*Coupling[M, {}, 0]^4) + (27*hbar*Coupling[A, {}, 0]^2*
     Coupling[cH4, {}, 0]*Coupling[\[Kappa], {}, 0])/Coupling[M, {}, 0]^4 + 
   (11*hbar*Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]^2)/
    (2*Coupling[M, {}, 0]^4) - (hbar*Coupling[\[Kappa], {}, 0]^3)/
    (12*Coupling[M, {}, 0]^2) - (hbar*Coupling[A, {}, 0]^4*
     Coupling[\[Lambda]\[Phi], {}, 0])/(2*Coupling[M, {}, 0]^6) - 
   (hbar*Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
     Coupling[\[Lambda]\[Phi], {}, 0])/(2*Coupling[M, {}, 0]^4) - 
   (7*hbar*Coupling[A, {}, 0]^5*Coupling[\[Mu], {}, 0])/
    (4*Coupling[M, {}, 0]^8) + (Coupling[A, {}, 0]^3*Coupling[\[Mu], {}, 0])/
    (6*Coupling[M, {}, 0]^6) - (13*hbar*Coupling[A, {}, 0]^3*
     Coupling[cH4, {}, 0]*Coupling[\[Mu], {}, 0])/Coupling[M, {}, 0]^6 - 
   (5*hbar*Coupling[A, {}, 0]^3*Coupling[\[Kappa], {}, 0]*
     Coupling[\[Mu], {}, 0])/Coupling[M, {}, 0]^6 - 
   (hbar*Coupling[A, {}, 0]*Coupling[\[Kappa], {}, 0]^2*
     Coupling[\[Mu], {}, 0])/(4*Coupling[M, {}, 0]^4) + 
   (hbar*Coupling[A, {}, 0]^3*Coupling[\[Lambda]\[Phi], {}, 0]*
     Coupling[\[Mu], {}, 0])/(3*Coupling[M, {}, 0]^6) + 
   (hbar*Coupling[A, {}, 0]^4*Coupling[\[Mu], {}, 0]^2)/
    Coupling[M, {}, 0]^8 + (hbar*Coupling[A, {}, 0]^2*
     Coupling[\[Kappa], {}, 0]*Coupling[\[Mu], {}, 0]^2)/
    (2*Coupling[M, {}, 0]^6) - (hbar*Coupling[A, {}, 0]^3*
     Coupling[\[Mu], {}, 0]^3)/(6*Coupling[M, {}, 0]^8) + 
   (40*hbar*Coupling[A, {}, 0]^2*Coupling[cH4, {}, 0]^2*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4 + 
   (10*hbar*Coupling[A, {}, 0]^2*Coupling[cH4, {}, 0]*Coupling[gL, {}, 0]^2*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(3*Coupling[M, {}, 0]^4) + 
   (3*hbar*Coupling[A, {}, 0]^4*Coupling[\[Kappa], {}, 0]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(2*Coupling[M, {}, 0]^6) + 
   (18*hbar*Coupling[A, {}, 0]^2*Coupling[cH4, {}, 0]*
     Coupling[\[Kappa], {}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
    Coupling[M, {}, 0]^4 + (3*hbar*Coupling[A, {}, 0]^2*
     Coupling[\[Kappa], {}, 0]^2*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
    Coupling[M, {}, 0]^4 - (hbar*Coupling[A, {}, 0]^2*
     Coupling[\[Kappa], {}, 0]*Coupling[\[Lambda]\[Phi], {}, 0]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(4*Coupling[M, {}, 0]^4) - 
   (hbar*Coupling[A, {}, 0]^5*Coupling[\[Mu], {}, 0]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(2*Coupling[M, {}, 0]^8) - 
   (6*hbar*Coupling[A, {}, 0]^3*Coupling[cH4, {}, 0]*Coupling[\[Mu], {}, 0]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^6 - 
   (hbar*Coupling[A, {}, 0]^3*Coupling[\[Kappa], {}, 0]*
     Coupling[\[Mu], {}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
    Coupling[M, {}, 0]^6 + (hbar*Coupling[A, {}, 0]^3*
     Coupling[\[Lambda]\[Phi], {}, 0]*Coupling[\[Mu], {}, 0]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(12*Coupling[M, {}, 0]^6), 
 Coupling[cHB, {}, 0] -> (hbar*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2)/
   (12*Coupling[M, {}, 0]^4), Coupling[cHBox, {}, 0] -> 
  (19*hbar*Coupling[A, {}, 0]^4)/(12*Coupling[M, {}, 0]^6) - 
   Coupling[A, {}, 0]^2/(2*Coupling[M, {}, 0]^4) - 
   (31*hbar*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2)/
    (24*Coupling[M, {}, 0]^4) - (31*hbar*Coupling[A, {}, 0]^2*
     Coupling[gY, {}, 0]^2)/(72*Coupling[M, {}, 0]^4) + 
   (25*hbar*Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0])/
    (12*Coupling[M, {}, 0]^4) - (hbar*Coupling[\[Kappa], {}, 0]^2)/
    (24*Coupling[M, {}, 0]^2) - (hbar*Coupling[A, {}, 0]^2*
     Coupling[\[Lambda]\[Phi], {}, 0])/(2*Coupling[M, {}, 0]^4) - 
   (13*hbar*Coupling[A, {}, 0]^3*Coupling[\[Mu], {}, 0])/
    (12*Coupling[M, {}, 0]^6) - (5*hbar*Coupling[A, {}, 0]*
     Coupling[\[Kappa], {}, 0]*Coupling[\[Mu], {}, 0])/
    (12*Coupling[M, {}, 0]^4) + (11*hbar*Coupling[A, {}, 0]^2*
     Coupling[\[Mu], {}, 0]^2)/(24*Coupling[M, {}, 0]^6) + 
   (hbar*Coupling[A, {}, 0]^4*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
    Coupling[M, {}, 0]^6 - (5*hbar*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(4*Coupling[M, {}, 0]^4) - 
   (5*hbar*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(12*Coupling[M, {}, 0]^4) + 
   (2*hbar*Coupling[A, {}, 0]^2*Coupling[\[Kappa], {}, 0]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/Coupling[M, {}, 0]^4 - 
   (hbar*Coupling[A, {}, 0]^2*Coupling[\[Lambda]\[Phi], {}, 0]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(2*Coupling[M, {}, 0]^4), 
 Coupling[cHBt, {}, 0] -> 0, 
 Coupling[cHd, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
  (-5*hbar*Bar[Coupling[cHqd, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[cHqd, {Index[d$$1, Flavor], 
       Index[i2, Flavor]}, 0])/(8*Coupling[M, {}, 0]^4) + 
   (17*hbar*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
     Delta[Index[i1, Flavor], Index[i2, Flavor]])/
    (216*Coupling[M, {}, 0]^4) - 
   (hbar*Bar[Coupling[cHqd, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[cHqd, {Index[d$$1, Flavor], 
       Index[i2, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
    (4*Coupling[M, {}, 0]^4) + (hbar*Coupling[A, {}, 0]^2*
     Coupling[gY, {}, 0]^2*Delta[Index[i1, Flavor], Index[i2, Flavor]]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(36*Coupling[M, {}, 0]^4), 
 Coupling[cHD, {}, 0] -> 
  (-31*hbar*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2)/
    (18*Coupling[M, {}, 0]^4) - (5*hbar*Coupling[A, {}, 0]^2*
     Coupling[gY, {}, 0]^2*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
    (3*Coupling[M, {}, 0]^4), 
 Coupling[cHe, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
  (-5*hbar*Bar[Coupling[cHle, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[cHle, {Index[d$$1, Flavor], 
       Index[i2, Flavor]}, 0])/(8*Coupling[M, {}, 0]^4) + 
   (17*hbar*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
     Delta[Index[i1, Flavor], Index[i2, Flavor]])/(72*Coupling[M, {}, 0]^4) - 
   (hbar*Bar[Coupling[cHle, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[cHle, {Index[d$$1, Flavor], 
       Index[i2, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
    (4*Coupling[M, {}, 0]^4) + (hbar*Coupling[A, {}, 0]^2*
     Coupling[gY, {}, 0]^2*Delta[Index[i1, Flavor], Index[i2, Flavor]]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(12*Coupling[M, {}, 0]^4), 
 Coupling[cHG, {}, 0] -> 0, Coupling[cHGt, {}, 0] -> 0, 
 Coupling[cHl1, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
  (5*hbar*Bar[Coupling[cHle, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[cHle, {Index[i1, Flavor], 
       Index[d$$1, Flavor]}, 0])/(16*Coupling[M, {}, 0]^4) + 
   (17*hbar*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
     Delta[Index[i1, Flavor], Index[i2, Flavor]])/
    (144*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Coupling[cHle, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[cHle, {Index[i1, Flavor], 
       Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
    (8*Coupling[M, {}, 0]^4) + (hbar*Coupling[A, {}, 0]^2*
     Coupling[gY, {}, 0]^2*Delta[Index[i1, Flavor], Index[i2, Flavor]]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(24*Coupling[M, {}, 0]^4), 
 Coupling[cHl3, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
  (5*hbar*Bar[Coupling[cHle, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[cHle, {Index[i1, Flavor], 
       Index[d$$1, Flavor]}, 0])/(16*Coupling[M, {}, 0]^4) - 
   (17*hbar*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
     Delta[Index[i1, Flavor], Index[i2, Flavor]])/
    (144*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Coupling[cHle, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[cHle, {Index[i1, Flavor], 
       Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
    (8*Coupling[M, {}, 0]^4) - (hbar*Coupling[A, {}, 0]^2*
     Coupling[gL, {}, 0]^2*Delta[Index[i1, Flavor], Index[i2, Flavor]]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(24*Coupling[M, {}, 0]^4), 
 Coupling[cHq1, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
  (5*hbar*Bar[Coupling[cHqd, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[cHqd, {Index[i1, Flavor], 
       Index[d$$1, Flavor]}, 0])/(16*Coupling[M, {}, 0]^4) - 
   (5*hbar*Bar[Coupling[cHqu, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[cHqu, {Index[i1, Flavor], 
       Index[d$$1, Flavor]}, 0])/(16*Coupling[M, {}, 0]^4) - 
   (17*hbar*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
     Delta[Index[i1, Flavor], Index[i2, Flavor]])/
    (432*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Coupling[cHqd, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[cHqd, {Index[i1, Flavor], 
       Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
    (8*Coupling[M, {}, 0]^4) - 
   (hbar*Bar[Coupling[cHqu, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[cHqu, {Index[i1, Flavor], 
       Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
    (8*Coupling[M, {}, 0]^4) - (hbar*Coupling[A, {}, 0]^2*
     Coupling[gY, {}, 0]^2*Delta[Index[i1, Flavor], Index[i2, Flavor]]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(72*Coupling[M, {}, 0]^4), 
 Coupling[cHq3, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
  (5*hbar*Bar[Coupling[cHqd, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[cHqd, {Index[i1, Flavor], 
       Index[d$$1, Flavor]}, 0])/(16*Coupling[M, {}, 0]^4) + 
   (5*hbar*Bar[Coupling[cHqu, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[cHqu, {Index[i1, Flavor], 
       Index[d$$1, Flavor]}, 0])/(16*Coupling[M, {}, 0]^4) - 
   (17*hbar*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2*
     Delta[Index[i1, Flavor], Index[i2, Flavor]])/
    (144*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Coupling[cHqd, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[cHqd, {Index[i1, Flavor], 
       Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
    (8*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Coupling[cHqu, {Index[i2, Flavor], Index[d$$1, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[cHqu, {Index[i1, Flavor], 
       Index[d$$1, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
    (8*Coupling[M, {}, 0]^4) - (hbar*Coupling[A, {}, 0]^2*
     Coupling[gL, {}, 0]^2*Delta[Index[i1, Flavor], Index[i2, Flavor]]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(24*Coupling[M, {}, 0]^4), 
 Coupling[cHu, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
  (5*hbar*Bar[Coupling[cHqu, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[cHqu, {Index[d$$1, Flavor], 
       Index[i2, Flavor]}, 0])/(8*Coupling[M, {}, 0]^4) - 
   (17*hbar*Coupling[A, {}, 0]^2*Coupling[gY, {}, 0]^2*
     Delta[Index[i1, Flavor], Index[i2, Flavor]])/
    (108*Coupling[M, {}, 0]^4) + 
   (hbar*Bar[Coupling[cHqu, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[cHqu, {Index[d$$1, Flavor], 
       Index[i2, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
    (4*Coupling[M, {}, 0]^4) - (hbar*Coupling[A, {}, 0]^2*
     Coupling[gY, {}, 0]^2*Delta[Index[i1, Flavor], Index[i2, Flavor]]*
     Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/(18*Coupling[M, {}, 0]^4), 
 Coupling[cHud, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
  (-5*hbar*Bar[Coupling[cHqu, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[cHqd, {Index[d$$1, Flavor], 
       Index[i2, Flavor]}, 0])/(4*Coupling[M, {}, 0]^4) - 
   (hbar*Bar[Coupling[cHqu, {Index[d$$1, Flavor], Index[i1, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[cHqd, {Index[d$$1, Flavor], 
       Index[i2, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
    (2*Coupling[M, {}, 0]^4), Coupling[cHW, {}, 0] -> 
  (hbar*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]^2)/
   (12*Coupling[M, {}, 0]^4), Coupling[cHWB, {}, 0] -> 
  (hbar*Coupling[A, {}, 0]^2*Coupling[gL, {}, 0]*Coupling[gY, {}, 0])/
   (6*Coupling[M, {}, 0]^4), Coupling[cHWt, {}, 0] -> 0, 
 Coupling[cHWtB, {}, 0] -> 0, 
 Coupling[cld, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 0, 
 Coupling[cle, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 
  -1/12*(hbar*Bar[Coupling[cHle, {Index[i2, Flavor], Index[i3, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[cHle, {Index[i1, Flavor], 
       Index[i4, Flavor]}, 0])/Coupling[M, {}, 0]^4, 
 Coupling[cledq, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 
  (hbar*Bar[Coupling[cHqd, {Index[i4, Flavor], Index[i3, Flavor]}, 0]]*
    Coupling[A, {}, 0]^2*Coupling[cHle, {Index[i1, Flavor], 
      Index[i2, Flavor]}, 0])/(6*Coupling[M, {}, 0]^4), 
 Coupling[clequ1, {Index[i1_, Flavor], Index[i2_, Flavor], 
    Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 
  -1/6*(hbar*Coupling[A, {}, 0]^2*Coupling[cHle, {Index[i1, Flavor], 
       Index[i2, Flavor]}, 0]*Coupling[cHqu, {Index[i3, Flavor], 
       Index[i4, Flavor]}, 0])/Coupling[M, {}, 0]^4, 
 Coupling[clequ3, {Index[i1_, Flavor], Index[i2_, Flavor], 
    Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 0, 
 Coupling[cll, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 0, 
 Coupling[clq1, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 0, 
 Coupling[clq3, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 0, 
 Coupling[clu, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 0, 
 Coupling[cqd1, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 
  -1/36*(hbar*Bar[Coupling[cHqd, {Index[i2, Flavor], Index[i3, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[cHqd, {Index[i1, Flavor], 
       Index[i4, Flavor]}, 0])/Coupling[M, {}, 0]^4, 
 Coupling[cqd8, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 
  -1/6*(hbar*Bar[Coupling[cHqd, {Index[i2, Flavor], Index[i3, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[cHqd, {Index[i1, Flavor], 
       Index[i4, Flavor]}, 0])/Coupling[M, {}, 0]^4, 
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
  -1/36*(hbar*Bar[Coupling[cHqu, {Index[i2, Flavor], Index[i3, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[cHqu, {Index[i1, Flavor], 
       Index[i4, Flavor]}, 0])/Coupling[M, {}, 0]^4, 
 Coupling[cqu8, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 
  -1/6*(hbar*Bar[Coupling[cHqu, {Index[i2, Flavor], Index[i3, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[cHqu, {Index[i1, Flavor], 
       Index[i4, Flavor]}, 0])/Coupling[M, {}, 0]^4, 
 Coupling[cquqd1, {Index[i1_, Flavor], Index[i2_, Flavor], 
    Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 
  (hbar*Coupling[A, {}, 0]^2*Coupling[cHqd, {Index[i3, Flavor], 
      Index[i4, Flavor]}, 0]*Coupling[cHqu, {Index[i1, Flavor], 
      Index[i2, Flavor]}, 0])/(6*Coupling[M, {}, 0]^4), 
 Coupling[cquqd8, {Index[i1_, Flavor], Index[i2_, Flavor], 
    Index[i3_, Flavor], Index[i4_, Flavor]}, 0] -> 0, 
 Coupling[cuB, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 0, 
 Coupling[cud1, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 0, 
 Coupling[cud8, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 0, 
 Coupling[cuG, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 0, 
 Coupling[cuH, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 
  (5*hbar*Bar[Coupling[cHqu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[cHqu, {Index[d$$1, Flavor], 
       Index[i2, Flavor]}, 0]*Coupling[cHqu, {Index[i1, Flavor], 
       Index[d$$2, Flavor]}, 0])/(4*Coupling[M, {}, 0]^4) - 
   (29*hbar*Coupling[A, {}, 0]^2*Coupling[cH4, {}, 0]*
     Coupling[cHqu, {Index[i1, Flavor], Index[i2, Flavor]}, 0])/
    (6*Coupling[M, {}, 0]^4) - (31*hbar*Coupling[A, {}, 0]^2*
     Coupling[cHqu, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[gL, {}, 0]^2)/(36*Coupling[M, {}, 0]^4) - 
   (3*hbar*Coupling[A, {}, 0]^2*Coupling[cHqu, {Index[i1, Flavor], 
       Index[i2, Flavor]}, 0]*Coupling[\[Kappa], {}, 0])/
    (4*Coupling[M, {}, 0]^4) + (hbar*Coupling[A, {}, 0]^3*
     Coupling[cHqu, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[\[Mu], {}, 0])/(4*Coupling[M, {}, 0]^6) + 
   (3*hbar*Bar[Coupling[cHqu, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]*
     Coupling[A, {}, 0]^2*Coupling[cHqu, {Index[d$$1, Flavor], 
       Index[i2, Flavor]}, 0]*Coupling[cHqu, {Index[i1, Flavor], 
       Index[d$$2, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
    (2*Coupling[M, {}, 0]^4) - (hbar*Coupling[A, {}, 0]^2*
     Coupling[cH4, {}, 0]*Coupling[cHqu, {Index[i1, Flavor], 
       Index[i2, Flavor]}, 0]*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
    Coupling[M, {}, 0]^4 - (5*hbar*Coupling[A, {}, 0]^2*
     Coupling[cHqu, {Index[i1, Flavor], Index[i2, Flavor]}, 0]*
     Coupling[gL, {}, 0]^2*Log[\[Mu]bar2/Coupling[M, {}, 0]^2])/
    (6*Coupling[M, {}, 0]^4), 
 Coupling[cuu, {Index[i1_, Flavor], Index[i2_, Flavor], Index[i3_, Flavor], 
    Index[i4_, Flavor]}, 0] -> 0, 
 Coupling[cuW, {Index[i1_, Flavor], Index[i2_, Flavor]}, 0] -> 0, 
 Coupling[cW, {}, 0] -> 0, Coupling[cWt, {}, 0] -> 0}
