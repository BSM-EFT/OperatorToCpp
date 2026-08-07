(* ::Package:: *)

(* ::Subsection:: *)
(*Load Matchete to access the EvaluateLoopFunctions[ ] function*)


SetDirectory[NotebookDirectory[]];
Needs["Matchete`"];


(* ::Subsection:: *)
(*Load one-loop matched result from (globally specified) path*)


matchedResultLocal=Import[$MatchedResultsPath];


(* ::Subsection:: *)
(*Extract out all varieties of loop functions in the matched results*)


ExtractAllLoopFuncExps[exprList_]:=Module[{expr,expComb,expCombs},
	expCombs=List[];
	Do[
		expr=exprList[[i]];
		expr=expr/.LF[{x__},{y__}]:>Exponents[y];
		If[Length[List@@expr]>0,
			expComb=Select[List@@expr,Head[#]===Exponents&]/.{Exponents[y__]:>y};
			If[Length[expComb]>0,AppendTo[expCombs,expComb]],
			Continue[]]
		,{i,1,Length[exprList]}];
	expCombs=DeleteDuplicates[Sort[expCombs]];
	Return[expCombs];
]


ExpCombs=ExtractAllLoopFuncExps[Total[Values[matchedResultLocal]]]
Clear[matchedResultLocal]


(* ::Subsection:: *)
(*Collect new loop-functions which appear when some masses are degenerate*)


SelectiveSum[list1_,indices_]:=Module[{sum,newlist,idxRev},
	sum=0;
	Do[sum=sum+list1[[indices[[i]]]],{i,1,Length[indices]}];
	newlist=ReplacePart[list1,First[indices]->sum];
	idxRev=Reverse[indices];
	Do[newlist=Drop[newlist,{idxRev[[i]]}],{i,1,Length[indices]-1}];
	Return[newlist]
]

EvaluateAllSums[lst_]:=Module[{len,idxList,sums},
	len=Length[lst]-1;
	idxList={};
	Do[idxList=Join[idxList,Subsets[Range[len],{i}]],{i,2,len}];

	sums={};
	Do[AppendTo[sums,SelectiveSum[lst,idxList[[i]]]],{i,1,Length[idxList]}];
	Return[sums]
]

UpdateExpCombs[ExpCombs_]:=Module[{summedExps,newCombs,updatedCombs},
	newCombs={};
	Do[
	summedExps=EvaluateAllSums[ExpCombs[[i]]]//DeleteDuplicates;
	newCombs=Join[newCombs,summedExps]//Sort//DeleteDuplicates,
	{i,1,Length[ExpCombs]}
	];
	updatedCombs=Union[ExpCombs,newCombs]
]


Clear[oldCombs,newCombs,AllExpCombs];
oldCombs = ExpCombs;
newCombs = UpdateExpCombs[ExpCombs];

While[
	oldCombs!=newCombs,
	oldCombs=newCombs;
	newCombs=UpdateExpCombs[newCombs];
]
AllExpCombs = newCombs;
AllExpCombs//Length


(* ::Subsubsection:: *)
(*Replacement rules to be added to Operator_Export.m*)


$ReplLoopFuncExpns={};
Do[AppendTo[$ReplLoopFuncExpns,AllExpCombs[[i]]->i],{i,1,Length[AllExpCombs]}];


(*Export["LFRules.m",$ReplLoopFuncExpns];*)


(* ::Subsection:: *)
(*Helper Functions to convert analytical formulas into a form consistent with C++*)


GenerateMassArgList[numExps_]:=Module[{mArgList},
	mArgList=List[];
	Do[AppendTo[mArgList,ToExpression["m"<>ToString[i]]],{i,1,numExps-1}];
	Return[mArgList];
]

GenerateMassRepRules[numExps_]:=Module[{mRules},
	mRules={};
	Do[AppendTo[mRules,"m"<>ToString[i]->"m"<>"["<>ToString[i-1]<>"]"],{i,1,numExps-1}];
	Return[mRules];
]

ConvertFormulaToString[ExpComb_]:=Module[{str},
	str=ToString[CForm[EvaluateLoopFunctions[LF[GenerateMassArgList[Length[ExpComb]],ExpComb]]/.\[Mu]bar2->mubarsq]];
	str=StringReplace[str,{"Power" -> "pow", "Log"->"log",".*"->"*", ". "->".0 ", ".)"->".0)" , "+ -" -> "-"}];
	str=StringReplace[str,GenerateMassRepRules[Length[ExpComb]]];
	Return[str]
]

ConvertFunctionCallToString[ExpComb_,rule_]:=Module[{str},
	str=ToString[
	CForm[LF[GenerateMassArgList[Length[ExpComb]],ExpComb]/.rule/.$ReplLoopFuncExpns]];
	str=StringReplace[str,"List("~~x__~~"),":>"{"~~x~~"},"];
	str=StringReplace[str,")":>",mubarsq)"];
	Return[str]
]

CreateIfElseBlock[ExpComb_,massArgList_]:=Module[{pairs,repRules,rules,strList},
	(* Remember if there is only 1 mass then there is no if, no else-if and no else simply a return statement. *)
	pairs=Subsets[massArgList,{2}];
	repRules={a_,b_}:>b->a;
	rules=pairs/.repRules;

	strList={};
	(* Create the first if branch *)
	AppendTo[
		strList,
		StringReplace[StringJoin["if (rel_diff(",ToString[pairs[[1]][[1]]],",",ToString[pairs[[1]][[2]]],") <= 5e-3) { return ", ConvertFunctionCallToString[ExpComb,rules[[1]]],"; }"  ],GenerateMassRepRules[Length[ExpComb]]]
	];

	(* Create the else-if branches *)
	If[
		Length[pairs]>1,
		Do[AppendTo[strList,StringReplace[StringJoin["else if (rel_diff(",ToString[pairs[[i]][[1]]],",",ToString[pairs[[i]][[2]]],") <= 5e-3) { return ", ConvertFunctionCallToString[ExpComb,rules[[i]]],"; }"  ],GenerateMassRepRules[Length[ExpComb]]]],{i,2,Length[pairs]}]
	];

	(* Create the non-degenerate else branch *)
	AppendTo[strList,StringJoin["else { return LFh",ToString[(ExpComb/.$ReplLoopFuncExpns)],"(m, mubarsq); }"]];
	Return[strList]
]


(* ::Subsubsection:: *)
(*Dictionary to enable accessing individual loop functions using a C++ switch statement*)


dict1 = Association[];
Do[
	massArgs=GenerateMassArgList[Length[AllExpCombs[[i]]]];
	strList=If[
		Length[massArgs]==1,
		{StringJoin["return LFh",ToString[(AllExpCombs[[i]]/.$ReplLoopFuncExpns)],"(m, mubarsq); "]},
		CreateIfElseBlock[AllExpCombs[[i]],massArgs]
];
AppendTo[dict1,ToString[i]->strList],{i,1,Length[AllExpCombs]}]


(* ::Subsection:: *)
(*Create the C++ source and header files*)


(* ::Subsubsection:: *)
(*Header file - LF.h*)


declaration="std::complex<double> LF(const std::vector<std::complex<double> >& masses, int code, double mubarsq);";


CreateLFHeaderFile[]:=Module[{path,fullPath,line},
	path = FileNameJoin[{NotebookDirectory[],"include"}];
	If[!DirectoryQ[path], CreateDirectory[path]];
	fullPath = FileNameJoin[{path,"LF.h"}];
	line = OpenWrite[fullPath];
	WriteLine[line, "#pragma once"];
	WriteLine[line, "#include <vector>"];
	WriteLine[line, "#include <complex>"];
	WriteLine[line, ""];
	WriteLine[line, "/**"];
	WriteLine[line, " * Evaluate loop functions."];
	WriteLine[line, " *"];
	WriteLine[line, " * @param masses A vector of mass parameters involved in the loop function evaluation."];
	WriteLine[line, " * @param code A unique identifier specifying the specific loop function variety."];
	WriteLine[line, " * @param mubarsq The mass scale at which the evaluation is to occur."];
	WriteLine[line, " * @return The evaluated numerical value."];
	WriteLine[line, " */"];
	WriteLine[line, declaration];
	Close[line];
]


(* ::Subsubsection:: *)
(*Header file - LF_helper.h*)


CreateLFHelperHeaderFile[]:=Module[{path,fullPath,line,str},
	path = FileNameJoin[{NotebookDirectory[],"include"}];
	If[!DirectoryQ[path], CreateDirectory[path]];
	fullPath = FileNameJoin[{path,"LF_helper.h"}];
	line = OpenWrite[fullPath];
	WriteLine[line, "#pragma once"];
	WriteLine[line, "#include <vector>"];
	WriteLine[line, "#include <complex>"];
	WriteLine[line, ""];
	WriteLine[line, "double rel_diff(std::complex<double> a, std::complex<double> b);"];
	WriteLine[line, ""];
	Do[WriteLine[line, "std::complex<double> LFh"<>ToString[i]<>"(const std::vector<std::complex<double> >& m, double mubarsq);"],{i,1,Length[AllExpCombs]}];
	WriteLine[line, ""];
	Close[line];
]


(* ::Subsubsection:: *)
(*Source file - LF.cpp*)


declL1 = "#include \"LF.h\"";
declL2 = "#include \"LF_helper.h\"";
doc = "// elaborate definitions of loop-functions in terms of masses";
defL1 = "std::complex<double> LF(const std::vector<std::complex<double> >& m, int code, double mubarsq) {";
defL2 = "    switch(code) {";
defClosingBraceSwitch = "    }";
defClosingBrace = "}";


CreateLFSourceFile[]:=Module[{path,fullPath,line,ifElseBranches},
	path = FileNameJoin[{NotebookDirectory[],"lib"}];
	If[!DirectoryQ[path], CreateDirectory[path]];
	fullPath = FileNameJoin[{path,"LF.cpp"}];
	line = OpenWrite[fullPath];
	WriteLine[line, declL1];
	WriteLine[line, declL2];
	WriteLine[line, ""];
	WriteLine[line, doc];
	WriteLine[line, defL1];
	WriteLine[line, defL2];
	Do[
		WriteLine[line, "        case "<>Keys[dict1][[i]]<>": " ];
		ifElseBranches=Values[dict1][[i]];
		Do[
			WriteLine[line, "            "<>ifElseBranches[[j]]],
			{j,1,Length[ifElseBranches]}
		];
		WriteLine[line, "            break;"],
		{i,1,Length[dict1]}
	];
	WriteLine[line, defClosingBraceSwitch];
	WriteLine[line, "    return 0.0;"];
	WriteLine[line, defClosingBrace];
	Close[line];
]


(* ::Subsubsection:: *)
(*Source file - LF_helper.cpp*)


CreateLFHelperSourceFile[]:=Module[{path,fullPath,line,str},
	path = FileNameJoin[{NotebookDirectory[],"lib"}];
	If[!DirectoryQ[path], CreateDirectory[path]];
	fullPath = FileNameJoin[{path,"LF_helper.cpp"}];
	line = OpenWrite[fullPath];
	WriteLine[line, "#include \"LF_helper.h\""];
	WriteLine[line, "#include \"complex_math.h\""];
	WriteLine[line, "#include <vector>"];
	WriteLine[line, "#include <complex>"];
	WriteLine[line, "#include <algorithm>"];
	WriteLine[line, "#include <cmath>"];
	WriteLine[line, ""];
	WriteLine[line, "double rel_diff(std::complex<double> a, std::complex<double> b) {"];
	WriteLine[line, "    return std::abs(a-b) / std::min(std::abs(a),std::abs(b));"];
	WriteLine[line, "}"];
	WriteLine[line, ""];
	Do[
		WriteLine[line, "std::complex<double> LFh"<>ToString[i]<>"(const std::vector<std::complex<double> >& m, double mubarsq) {"];
		WriteLine[line, "    return "<>ConvertFormulaToString[AllExpCombs[[i]]]<>";"];
		WriteLine[line, "}\n"],{i,1,Length[AllExpCombs]}];
	WriteLine[line, ""];
	Close[line];
]


CreateLFHeaderFile[]
CreateLFSourceFile[]


CreateLFHelperHeaderFile[]
CreateLFHelperSourceFile[]


Remove["Matchete`*"]
