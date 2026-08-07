(* ::Package:: *)

(* ::Title:: *)
(*SMEFT model file*)


(* ::Text:: *)
(*This file differs from the predefined "SMEFT_Warsaw.m" model fine, by a different sign convention for the covariant derivatives. *)
(*The coefficients of operators that contain an odd number of field-strength tensors differ by an overall minus sign with respect to the definition of these coefficients in https://arxiv.org/abs/1008.4884.*)


ParentModel["SM"]


DefineCoupling[cllHH, Indices->{Flavor, Flavor}, Symmetries-> {SymmetricPermutation[2, 1]},
	NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(llHH\)]\)"];


DefineCoupling[cG,    SelfConjugate-> True, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(G\)]\)"];
DefineCoupling[cGt,   SelfConjugate-> True, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), OverscriptBox[\(G\), \(~\)]]\)"];
DefineCoupling[cW,    SelfConjugate-> True, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(W\)]\)"];
DefineCoupling[cWt,   SelfConjugate-> True, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), OverscriptBox[\(W\), \(~\)]]\)"];
DefineCoupling[cHG,   SelfConjugate-> True, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(HG\)]\)"];
DefineCoupling[cHGt,  SelfConjugate-> True, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(H \*OverscriptBox[\(G\), \(~\)]\)]\)"];
DefineCoupling[cHW,   SelfConjugate-> True, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(HW\)]\)"];
DefineCoupling[cHWt,  SelfConjugate-> True, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(H \*OverscriptBox[\(W\), \(~\)]\)]\)"];
DefineCoupling[cHB,   SelfConjugate-> True, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(HB\)]\)"];
DefineCoupling[cHBt,  SelfConjugate-> True, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(H \*OverscriptBox[\(B\), \(~\)]\)]\)"];
DefineCoupling[cHWB,  SelfConjugate-> True, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(HWB\)]\)"];
DefineCoupling[cHWtB, SelfConjugate-> True, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(H \*OverscriptBox[\(W\), \(~\)] B\)]\)"];
DefineCoupling[cH,    SelfConjugate-> True, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(H\)]\)"];
DefineCoupling[cHBox, SelfConjugate-> True, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(H \[Square]\)]\)"];
DefineCoupling[cHD,   SelfConjugate-> True, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(HD\)]\)"];


DefineCoupling[ceH, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(eH\)]\)"];
DefineCoupling[cuH, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(uH\)]\)"];
DefineCoupling[cdH, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(dH\)]\)"];
DefineCoupling[ceW, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(eW\)]\)"];
DefineCoupling[ceB, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(eB\)]\)"];
DefineCoupling[cuG, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(uG\)]\)"];
DefineCoupling[cuW, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(uW\)]\)"];
DefineCoupling[cuB, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(uB\)]\)"];
DefineCoupling[cdG, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(dG\)]\)"];
DefineCoupling[cdW, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(dW\)]\)"];
DefineCoupling[cdB, Indices-> {Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(dB\)]\)"];


DefineCoupling[cHl1, Indices-> {Flavor, Flavor}, SelfConjugate-> {2, 1}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(Hl\), \((1)\)]\)"];
DefineCoupling[cHl3, Indices-> {Flavor, Flavor}, SelfConjugate-> {2, 1}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(Hl\), \((3)\)]\)"];
DefineCoupling[cHe,  Indices-> {Flavor, Flavor}, SelfConjugate-> {2, 1}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(He\)]\)"];
DefineCoupling[cHq1, Indices-> {Flavor, Flavor}, SelfConjugate-> {2, 1}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(Hq\), \((1)\)]\)"];
DefineCoupling[cHq3, Indices-> {Flavor, Flavor}, SelfConjugate-> {2, 1}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(Hq\), \((3)\)]\)"];
DefineCoupling[cHu,  Indices-> {Flavor, Flavor}, SelfConjugate-> {2, 1}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(Hu\)]\)"];
DefineCoupling[cHd,  Indices-> {Flavor, Flavor}, SelfConjugate-> {2, 1}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(Hd\)]\)"];
DefineCoupling[cHud, Indices-> {Flavor, Flavor}, NiceForm->"\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(Hud\)]\)"];


DefineCoupling[cll,  Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, Symmetries-> {SymmetricPermutation[3, 4, 1, 2]}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(ll\)]\)"];
DefineCoupling[cqq1, Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, Symmetries-> {SymmetricPermutation[3, 4, 1, 2]}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(qq\), \((1)\)]\)"];
DefineCoupling[cqq3, Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, Symmetries-> {SymmetricPermutation[3, 4, 1, 2]}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(qq\), \((3)\)]\)"];
DefineCoupling[clq1, Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(lq\), \((1)\)]\)"];
DefineCoupling[clq3, Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(lq\), \((3)\)]\)"];
DefineCoupling[cee,  Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, Symmetries-> {SymmetricPermutation[3, 4, 1, 2], SymmetricPermutation[1, 4, 3, 2]}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(ee\)]\)"];
DefineCoupling[cuu,  Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, Symmetries-> {SymmetricPermutation[3, 4, 1, 2]}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(uu\)]\)"];
DefineCoupling[cdd,  Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, Symmetries-> {SymmetricPermutation[3, 4, 1, 2]}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(dd\)]\)"];
DefineCoupling[ceu,  Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(eu\)]\)"];
DefineCoupling[ced,  Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(ed\)]\)"];
DefineCoupling[cud1, Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(ud\), \((1)\)]\)"];
DefineCoupling[cud8, Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(ud\), \((8)\)]\)"];
DefineCoupling[cle,  Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(le\)]\)"];
DefineCoupling[clu,  Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(lu\)]\)"];
DefineCoupling[cld,  Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(ld\)]\)"];
DefineCoupling[cqe,  Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(qe\)]\)"];
DefineCoupling[cqu1, Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(qu\), \((1)\)]\)"];
DefineCoupling[cqu8, Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(qu\), \((8)\)]\)"];
DefineCoupling[cqd1, Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(qd\), \((1)\)]\)"];
DefineCoupling[cqd8, Indices-> {Flavor, Flavor, Flavor, Flavor}, SelfConjugate-> {2, 1, 4, 3}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(qd\), \((8)\)]\)"];


DefineCoupling[cduq, Indices-> {Flavor, Flavor, Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(duq\)]\)"];
DefineCoupling[cqqu, Indices-> {Flavor, Flavor, Flavor, Flavor}, Symmetries-> {SymmetricPermutation[2, 1, 3, 4]}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(qqu\)]\)"];
DefineCoupling[cqqq, Indices-> {Flavor, Flavor, Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(qqq\)]\)"];
DefineCoupling[cduu, Indices-> {Flavor, Flavor, Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(duu\)]\)"];


DefineCoupling[cledq,  Indices-> {Flavor, Flavor, Flavor, Flavor}, NiceForm-> "\!\(\*SubscriptBox[\(\[ScriptCapitalC]\), \(ledq\)]\)"];
DefineCoupling[cquqd1, Indices-> {Flavor, Flavor, Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(quqd\), \((1)\)]\)"];
DefineCoupling[cquqd8, Indices-> {Flavor, Flavor, Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(quqd\), \((8)\)]\)"];
DefineCoupling[clequ1, Indices-> {Flavor, Flavor, Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(lequ\), \((1)\)]\)"];
DefineCoupling[clequ3, Indices-> {Flavor, Flavor, Flavor, Flavor}, NiceForm-> "\!\(\*SubsuperscriptBox[\(\[ScriptCapitalC]\), \(lequ\), \((3)\)]\)"];


Module[{p,r,s,t,i,j,k,m,n,o,J,K,L,\[Alpha],\[Beta],\[Delta],A,C,D,\[Nu],\[Rho],\[Theta],\[Eta],\[Kappa],
		\[Tau]SU2L,fSU2L,\[Epsilon]SU2L,TSU3c,fSU3c,\[Epsilon]SU3c,HermitianCD,
		Lagd5,Lagd60F,Lagd62F, Lagd62FH ,Lagd64F,Lagd64FH, Lagd6BNV},

	\[Tau]SU2L[Jadj_,ifund_,jfund_]:= 2 CG[gen[SU2L[fund]],{Jadj,ifund,jfund}];
	fSU2L[Iadj_,Jadj_,Kadj_]:=CG[fStruct[SU2L],{Iadj,Jadj,Kadj}];
	\[Epsilon]SU2L[ifund_,jfund_]:= CG[eps[SU2L],{ifund,jfund}];
	TSU3c[Aadj_,\[Alpha]fund_,\[Beta]fund_]:=CG[gen[SU3c[fund]],{Aadj,\[Alpha]fund,\[Beta]fund}];
	fSU3c[Aadj_,Badj_,Cadj_]:=CG[fStruct[SU3c],{Aadj,Badj,Cadj}];
	\[Epsilon]SU3c[\[Alpha]fund_,\[Beta]fund_,\[Gamma]fund_]:= CG[eps[SU3c],{\[Alpha]fund,\[Beta]fund,\[Gamma]fund}];
	HermitianCD[ind_,term1_,term2_]:= I term1 CD[ind,term2] - I CD[ind,term1] term2;

	Lagd5 = cllHH[p,r] Bar@\[Epsilon]SU2L[i,j]Bar@\[Epsilon]SU2L[k,m] H[i] H[k] Bar@CConj@l[j,p]**l[m,r];

	Lagd60F = cG[] fSU3c[A,C,D] FS[G, \[Nu], \[Rho], A] FS[G, \[Rho], \[Theta], C] FS[G, \[Theta], \[Nu], D]/gs[]^3 + 1/2 cGt[] LCTensor[\[Nu],\[Rho],\[Eta],\[Kappa]]fSU3c[A,C,D] FS[G, \[Eta], \[Kappa], A] FS[G, \[Rho], \[Theta], C] FS[G, \[Theta], \[Nu], D]/gs[]^3
		    + cW[] fSU2L[J,K,L] FS[W, \[Nu], \[Rho], J] FS[W, \[Rho], \[Theta], K] FS[W, \[Theta], \[Nu], L]/gL[]^3 + 1/2 cWt[] LCTensor[\[Nu],\[Rho],\[Eta],\[Kappa]]fSU2L[J,K,L] FS[W, \[Eta], \[Kappa], J] FS[W, \[Rho], \[Theta], K] FS[W, \[Theta], \[Nu], L]/gL[]^3
	        + cHG[] Bar@H[i] H[i] FS[G, \[Nu], \[Rho], A] FS[G, \[Nu], \[Rho], A]/gs[]^2 + 1/2 cHGt[] Bar@H[i] H[i] LCTensor[\[Nu],\[Rho],\[Eta],\[Kappa]] FS[G, \[Eta], \[Kappa], A] FS[G, \[Nu], \[Rho], A]/gs[]^2
	        + cHW[] Bar@H[i] H[i] FS[W, \[Nu], \[Rho], J] FS[W, \[Nu], \[Rho], J]/gL[]^2 + 1/2 cHWt[] Bar@H[i] H[i] LCTensor[\[Nu],\[Rho],\[Eta],\[Kappa]] FS[W, \[Eta], \[Kappa], J] FS[W, \[Nu], \[Rho], J]/gL[]^2
	        + cHB[] Bar@H[i] H[i] FS[B, \[Nu], \[Rho]] FS[B, \[Nu], \[Rho]]/gY[]^2 + 1/2 cHBt[] Bar@H[i] H[i] LCTensor[\[Nu],\[Rho],\[Eta],\[Kappa]] FS[B, \[Eta], \[Kappa]] FS[B, \[Nu], \[Rho]]/gY[]^2
	        + cHWB[] Bar@H[i] \[Tau]SU2L[J,i,j] H[j] FS[W, \[Nu], \[Rho], J] FS[B, \[Nu], \[Rho]]/(gL[]gY[]) + 1/2 cHWtB[] Bar@H[i] \[Tau]SU2L[J,i,j] H[j] LCTensor[\[Nu],\[Rho],\[Eta],\[Kappa]] FS[W, \[Eta], \[Kappa], J] FS[B, \[Nu], \[Rho]]/(gL[]gY[])
	        + cH[] Bar@H[i] H[i] Bar@H[j] H[j] Bar@H[k] H[k]
	        + cHBox[] Bar@H[i] H[i] CD[{\[Nu],\[Nu]},Bar@H[j] H[j]]
	        + cHD[] H[i] CD[\[Nu], Bar@H[i]] Bar@H[j] CD[\[Nu],H[j]];
	Lagd62F = ceH[p,r] Bar@H[i] H[i] H[j] Bar@l[j,p]**e[r]
			+ cuH[p,r] Bar@H[i] H[i] \[Epsilon]SU2L[j,k] Bar@H[k] Bar@q[\[Alpha],j,p]**u[\[Alpha],r]
			+ cdH[p,r] Bar@H[i] H[i] H[j] Bar@q[\[Alpha],j,p]**d[\[Alpha],r]
			+ ceW[p,r] Bar@l[i,p]**\[Sigma][\[Nu],\[Rho]]**e[r] \[Tau]SU2L[J,i,j] H[j] FS[W,\[Nu],\[Rho],J]/gL[]
			+ ceB[p,r] Bar@l[i,p]**\[Sigma][\[Nu],\[Rho]]**e[r] H[i] FS[B,\[Nu],\[Rho]]/gY[]
			+ cuG[p,r] Bar@q[\[Alpha],i,p]**\[Sigma][\[Nu],\[Rho]]**u[\[Beta],r]TSU3c[A,\[Alpha],\[Beta]] \[Epsilon]SU2L[i,j] Bar@H[j] FS[G,\[Nu],\[Rho],A]/gs[]
			+ cuW[p,r] Bar@q[\[Alpha],i,p]**\[Sigma][\[Nu],\[Rho]]**u[\[Alpha],r] \[Tau]SU2L[J,i,j]\[Epsilon]SU2L[j,k]Bar@H[k] FS[W,\[Nu],\[Rho],J]/gL[]
			+ cuB[p,r] Bar@q[\[Alpha],i,p]**\[Sigma][\[Nu],\[Rho]]**u[\[Alpha],r] \[Epsilon]SU2L[i,j] Bar@H[j] FS[B,\[Nu],\[Rho]]/gY[]
			+ cdG[p,r] Bar@q[\[Alpha],i,p]**\[Sigma][\[Nu],\[Rho]]**d[\[Beta],r] H[i] TSU3c[A,\[Alpha],\[Beta]]FS[G,\[Nu],\[Rho],A]/gs[]
			+ cdW[p,r] Bar@q[\[Alpha],i,p]**\[Sigma][\[Nu],\[Rho]]**d[\[Alpha],r] \[Tau]SU2L[J,i,j] H[j] FS[W,\[Nu],\[Rho],J]/gL[]
			+ cdB[p,r] Bar@q[\[Alpha],i,p]**\[Sigma][\[Nu],\[Rho]]**d[\[Alpha],r] H[i] FS[B,\[Nu],\[Rho]]/gY[];
	Lagd62FH = cHl1[p,r] HermitianCD[\[Nu],Bar@H[i],H[i]] Bar@l[j,p]**\[Gamma][\[Nu]]**l[j,r]
			+ cHl3[p,r] HermitianCD[\[Nu],Bar@H[i],\[Tau]SU2L[J,i,j] H[j]] \[Tau]SU2L[J,k,m] Bar@l[k,p]**\[Gamma][\[Nu]]**l[m,r]
			+ cHe[p,r] HermitianCD[\[Nu],Bar@H[i], H[i]] Bar@e[p]**\[Gamma][\[Nu]]**e[r]
			+ cHq1[p,r] HermitianCD[\[Nu],Bar@H[i],H[i]] Bar@q[\[Alpha],j,p]**\[Gamma][\[Nu]]**q[\[Alpha],j,r]
			+ cHq3[p,r] HermitianCD[\[Nu],Bar@H[i],\[Tau]SU2L[J,i,j] H[j]] \[Tau]SU2L[J,k,m] Bar@q[\[Alpha],k,p]**\[Gamma][\[Nu]]**q[\[Alpha],m,r]
			+ cHu[p,r] HermitianCD[\[Nu],Bar@H[i],H[i]] Bar@u[\[Alpha],p]**\[Gamma][\[Nu]]**u[\[Alpha],r]
			+ cHd[p,r] HermitianCD[\[Nu],Bar@H[i],H[i]] Bar@d[\[Alpha],p]**\[Gamma][\[Nu]]**d[\[Alpha],r]
			+ PlusHc[cHud[p,r] I Bar@\[Epsilon]SU2L[i,j] H[j] CD[\[Nu],H[i]] Bar@u[\[Alpha],p]**\[Gamma][\[Nu]]**d[\[Alpha],r]];
	Lagd64FH = cll[p,r,s,t] Bar@l[i,p]**\[Gamma][\[Nu]]**l[i,r] Bar@l[j,s]**\[Gamma][\[Nu]]**l[j,t]
			+ cqq1[p,r,s,t] Bar@q[\[Alpha],i,p]**\[Gamma][\[Nu]]**q[\[Alpha],i,r] Bar@q[\[Beta],j,s]**\[Gamma][\[Nu]]**q[\[Beta],j,t]
			+ cqq3[p,r,s,t] \[Tau]SU2L[J,i,j]Bar@q[\[Alpha],i,p]**\[Gamma][\[Nu]]**q[\[Alpha],j,r] \[Tau]SU2L[J,k,m]Bar@q[\[Beta],k,s]**\[Gamma][\[Nu]]**q[\[Beta],m,t]
			+ clq1[p,r,s,t] Bar@l[i,p]**\[Gamma][\[Nu]]**l[i,r] Bar@q[\[Alpha],j,s]**\[Gamma][\[Nu]]**q[\[Alpha],j,t]
			+ clq3[p,r,s,t] \[Tau]SU2L[J,i,j]Bar@l[i,p]**\[Gamma][\[Nu]]**l[j,r] \[Tau]SU2L[J,k,m]Bar@q[\[Alpha],k,s]**\[Gamma][\[Nu]]**q[\[Alpha],m,t]
			+ cee[p,r,s,t] Bar@e[p]**\[Gamma][\[Nu]]**e[r] Bar@e[s]**\[Gamma][\[Nu]]**e[t]
			+ cuu[p,r,s,t] Bar@u[\[Alpha],p]**\[Gamma][\[Nu]]**u[\[Alpha],r] Bar@u[\[Beta],s]**\[Gamma][\[Nu]]**u[\[Beta],t]
			+ cdd[p,r,s,t] Bar@d[\[Alpha],p]**\[Gamma][\[Nu]]**d[\[Alpha],r] Bar@d[\[Beta],s]**\[Gamma][\[Nu]]**d[\[Beta],t]
			+ ceu[p,r,s,t] Bar@e[p]**\[Gamma][\[Nu]]**e[r] Bar@u[\[Alpha],s]**\[Gamma][\[Nu]]**u[\[Alpha],t]
			+ ced[p,r,s,t] Bar@e[p]**\[Gamma][\[Nu]]**e[r] Bar@d[\[Alpha],s]**\[Gamma][\[Nu]]**d[\[Alpha],t]
			+ cud1[p,r,s,t] Bar@u[\[Alpha],p]**\[Gamma][\[Nu]]**u[\[Alpha],r] Bar@d[\[Beta],s]**\[Gamma][\[Nu]]**d[\[Beta],t]
			+ cud8[p,r,s,t] TSU3c[A,\[Alpha],\[Beta]] Bar@u[\[Alpha],p]**\[Gamma][\[Nu]]**u[\[Beta],r] TSU3c[A,\[Delta],\[Kappa]] Bar@d[\[Delta],s]**\[Gamma][\[Nu]]**d[\[Kappa],t]
			+ cle[p,r,s,t] Bar@l[i,p]**\[Gamma][\[Nu]]**l[i,r] Bar@e[s]**\[Gamma][\[Nu]]**e[t]
			+ clu[p,r,s,t] Bar@l[i,p]**\[Gamma][\[Nu]]**l[i,r] Bar@u[\[Alpha],s]**\[Gamma][\[Nu]]**u[\[Alpha],t]
			+ cld[p,r,s,t] Bar@l[i,p]**\[Gamma][\[Nu]]**l[i,r] Bar@d[\[Alpha],s]**\[Gamma][\[Nu]]**d[\[Alpha],t]
			+ cqe[p,r,s,t] Bar@q[\[Alpha],i,p]**\[Gamma][\[Nu]]**q[\[Alpha],i,r] Bar@e[s]**\[Gamma][\[Nu]]**e[t]
			+ cqu1[p,r,s,t] Bar@q[\[Alpha],i,p]**\[Gamma][\[Nu]]**q[\[Alpha],i,r] Bar@u[\[Beta],s]**\[Gamma][\[Nu]]**u[\[Beta],t]
			+ cqu8[p,r,s,t] TSU3c[A,\[Alpha],\[Beta]] Bar@q[\[Alpha],i,p]**\[Gamma][\[Nu]]**q[\[Beta],i,r] TSU3c[A,\[Delta],\[Kappa]] Bar@u[\[Delta],s]**\[Gamma][\[Nu]]**u[\[Kappa],t]
			+ cqd1[p,r,s,t] Bar@q[\[Alpha],i,p]**\[Gamma][\[Nu]]**q[\[Alpha],i,r] Bar@d[\[Beta],s]**\[Gamma][\[Nu]]**d[\[Beta],t]
			+ cqd8[p,r,s,t] TSU3c[A,\[Alpha],\[Beta]] Bar@q[\[Alpha],i,p]**\[Gamma][\[Nu]]**q[\[Beta],i,r] TSU3c[A,\[Delta],\[Kappa]] Bar@d[\[Delta],s]**\[Gamma][\[Nu]]**d[\[Kappa],t];
	Lagd64F = cledq[p,r,s,t] Bar@l[i,p]**e[r] Bar@d[\[Alpha],s]**q[\[Alpha],i,t]
			+ cquqd1[p,r,s,t] Bar@q[\[Alpha],i,p]**u[\[Alpha],r] \[Epsilon]SU2L[i,j] Bar@q[\[Beta],j,s]**d[\[Beta],t]
			+ cquqd8[p,r,s,t] TSU3c[A,\[Alpha],\[Beta]] Bar@q[\[Alpha],i,p]**u[\[Beta],r] \[Epsilon]SU2L[i,j] TSU3c[A,\[Delta],\[Kappa]] Bar@q[\[Delta],j,s]**d[\[Kappa],t]
			+ clequ1[p,r,s,t] Bar@l[i,p]**e[r] \[Epsilon]SU2L[i,j] Bar@q[\[Alpha],j,s]**u[\[Alpha],t]
			+ clequ3[p,r,s,t] Bar@l[i,p]**\[Sigma][\[Nu],\[Rho]]**e[r] \[Epsilon]SU2L[i,j] Bar@q[\[Alpha],j,s]**\[Sigma][\[Nu],\[Rho]]**u[\[Alpha],t];
	Lagd6BNV = cduq[p,r,s,t] Bar@\[Epsilon]SU3c[\[Alpha],\[Beta],\[Delta]]Bar@\[Epsilon]SU2L[i,j] CConj@Bar@d[\[Alpha],p]**u[\[Beta],r] CConj@Bar@q[\[Delta],i,s]**l[j,t]
			 + cqqu[p,r,s,t] Bar@\[Epsilon]SU3c[\[Alpha],\[Beta],\[Delta]]Bar@\[Epsilon]SU2L[i,j] CConj@Bar@q[\[Alpha],i,p]**q[\[Beta],j,r] CConj@Bar@u[\[Delta],s]**e[t]
			 + cqqq[p,r,s,t] Bar@\[Epsilon]SU3c[\[Alpha],\[Beta],\[Delta]]Bar@\[Epsilon]SU2L[i,j] Bar@\[Epsilon]SU2L[k,m] CConj@Bar@q[\[Alpha],i,p]**q[\[Beta],k,r] CConj@Bar@q[\[Delta],m,s]**l[j,t]
			 + cduu[p,r,s,t] Bar@\[Epsilon]SU3c[\[Alpha],\[Beta],\[Delta]] CConj@Bar@d[\[Alpha],p]**u[\[Beta],r] CConj@Bar@u[\[Delta],s]**e[t];

	PlusHc[Lagd5] + Lagd60F + Lagd62FH + Lagd64FH + PlusHc[Lagd62F + Lagd64F + Lagd6BNV]//RelabelIndices
]
