(* ::Package:: *)

(* ::Chapter:: *)
(*-- Package Context --*)


Needs["Matchete`"];
BeginPackage["EFTParser`", {"Global`", "Matchete`"}];


(* ::Subsection:: *)
(*Description of public functions*)


ReadEFTInfo::usage = "ReadEFTInfo[modelPath,EFTPath] receives the path (or name if already in the Matchete models-database) of the low-energy model and the corresponding EFT extension, and return an association containing the relevant metadata for each Wilson coefficient."
ExportWCInfo::usgae = "ExportWCInfo[EFT,Basis,WCInfo] admits the name of the EFT, the specific Basis and an association containing the attributes of each Wilson coefficient and stores the information within a JSON file, to be used at the Python frontend."


(* ::Subsection:: *)
(*External symbols*)


(* External symbols*)
Matchete`GetCouplings;
Matchete`SelfConjugate;
Matchete`Symmetries;
Matchete`Indices;
Matchete`ResetAll


Begin["`Private`"];


NotebookDirectory[];


(* ::Chapter:: *)
(*-- Function definitions --*)


(* ::Subsection:: *)
(*Function to read EFT.m file and generate the WCInfo dictionary*)


(* ::Subsubsection:: *)
(*Helper functions*)


(* Admits the path of the Low-energy model and the EFT, and returns a list of Wilson Coefficients *)
LoadCoeffs[modelPath_,EFTPath_]:=Module[{couplingsRen,couplingsAll,WilsonCoeffs},
	LoadModel[modelPath];
	couplingsRen =Keys[Matchete`GetCouplings[]];
	LoadModel[EFTPath];
	couplingsAll = Keys[Matchete`GetCouplings[]];
	WilsonCoeffs =Complement[couplingsAll,couplingsRen];
	Return[WilsonCoeffs];
];


(* Combine all index permutations *)
CollectPermutations[wc_]:=Module[{perms},
	perms=If[
		ListQ[wc[Matchete`SelfConjugate]],
		Join[{wc[Matchete`SelfConjugate]},Keys[wc[Matchete`Symmetries]]],
		Replace[Keys[wc[Matchete`Symmetries]],{ {{}}:>{} }]
	];
	Return[perms]
]



(* Read the WC data and store it as an association *)
CreateWCDict[wcs_]:=Module[{WCDict},
	WCDict = Association[];
	Do[
		wc=Matchete`GetCouplings[wcs[[i]]];
		WCDict[ToString[wcs[[i]]]]=Association[{
			"Nf"->Length[wc[Matchete`Indices]],
			"SelfConjugate"->If[wc[Matchete`SelfConjugate]===True,True,False],
			"Permutations"->CollectPermutations[wc]
		}],{i,1,Length[wcs]}
	];
	Return[SortBy[WCDict,#"Nf"&]];
]


(* ::Subsubsection:: *)
(*User-facing functions*)


ReadEFTInfo[modelPath_,EFTPath_]:=Module[{WCs,WCdict},
	Matchete`ResetAll[];
	WCs=LoadCoeffs[modelPath,EFTPath];
	WCdict=CreateWCDict[WCs];
	Return[WCdict];
]


(* ::Input::Initialization:: *)
ExportWCInfo[EFT_,Basis_,WCInfo_]:=Module[{filename,path},
filename=EFT<>"-"<>Basis<>".json";
path=FileNameJoin[{NotebookDirectory[],"py","wcInfo",filename}];
Export[path,WCInfo,"JSON","Compact"->2];
];


(* ::Chapter:: *)
(*-- Package Context --*)


End[];
EndPackage[];
