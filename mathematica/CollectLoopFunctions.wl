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


ReplLoopFuncExpns={};
Do[AppendTo[ReplLoopFuncExpns,AllExpCombs[[i]]->i],{i,1,Length[AllExpCombs]}];


Export["LFRules.m",ReplLoopFuncExpns];


(* ::Subsection:: *)
(*Helper Functions to convert analytical formulas into a form consistent with C++*)


GenerateMassArgList[numExps_]:=Module[{mArgList},
	mArgList=List[];
	Do[AppendTo[mArgList,ToExpression["m"<>ToString[i]]],{i,1,numExps-1}];
	Return[mArgList];
]

GenerateMassRepRules[numExps_]:=Module[{mRules},
	mRules={};
	Do[AppendTo[mRules,"m"<>ToString[i]->"masses"<>"["<>ToString[i-1]<>"]"],{i,1,numExps-1}];
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
	CForm[LF[GenerateMassArgList[Length[ExpComb]],ExpComb]/.rule/.ReplLoopFuncExpns]];
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
		StringReplace[StringJoin["if (abs(",ToString[pairs[[1]][[1]]]," - ",ToString[pairs[[1]][[2]]],") < 1e-6) { return ", ConvertFunctionCallToString[ExpComb,rules[[1]]],"; }"  ],GenerateMassRepRules[Length[ExpComb]]]
	];

	(* Create the else-if branches *)
	If[
		Length[pairs]>1,
		Do[AppendTo[strList,StringReplace[StringJoin["else if (abs(",ToString[pairs[[i]][[1]]]," - ",ToString[pairs[[i]][[2]]],") < 1e-6) { return ", ConvertFunctionCallToString[ExpComb,rules[[i]]],"; }"  ],GenerateMassRepRules[Length[ExpComb]]]],{i,2,Length[pairs]}]
	];

	(* Create the non-degenerate else branch *)
	AppendTo[strList,StringReplace[StringJoin["else { return ",ConvertFormulaToString[ExpComb], "; }"],GenerateMassRepRules[Length[ExpComb]]]];
	Return[strList]
]


(* ::Subsubsection:: *)
(*Dictionary to enable accessing individual loop functions using a C++ switch statement*)


dict1 = Association[];
Do[
	massArgs=GenerateMassArgList[Length[AllExpCombs[[i]]]];
	strList=If[
		Length[massArgs]==1,
		{StringJoin["return ",ConvertFormulaToString[AllExpCombs[[i]]], ";"]},
		CreateIfElseBlock[AllExpCombs[[i]],massArgs]
];
AppendTo[dict1,ToString[i]->strList],{i,1,Length[AllExpCombs]}]


(* ::Subsection:: *)
(*Create the C++ source and header files*)


(* ::Subsubsection:: *)
(*Header file*)


declaration="std::complex<double> LF(std::vector<std::complex<double> > masses, int code, double mubarsq);";


CreateLFHeaderFile[]:=Module[{path,fullPath,line},
	path = FileNameJoin[{ParentDirectory[NotebookDirectory[]],"include"}];
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
(*Source file*)


declL1 = "#include \"pch.h\"";
declL2 = "#include \"LF.h\"";
declL3 = "#include \"complex_math.h\"";
declL4 = "using std::vector;";
declL5 = "using std::complex;";
doc = "// elaborate definitions of loop-functions in terms of masses";
defL1 = "complex<double> LF(vector<complex<double> > masses, int code, double mubarsq) {";
defL2 = "    switch(code) {";
defClosingBraceSwitch = "    }";
defClosingBrace = "}";


CreateLFSourceFile[]:=Module[{path,fullPath,line,ifElseBranches},
	path = FileNameJoin[{ParentDirectory[NotebookDirectory[]],"lib"}];
	If[!DirectoryQ[path], CreateDirectory[path]];
	fullPath = FileNameJoin[{path,"LF.cpp"}];
	line = OpenWrite[fullPath];
	WriteLine[line, declL1];
	WriteLine[line, declL2];
	WriteLine[line, declL3];
	WriteLine[line, declL4];
	WriteLine[line, declL5];
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


CreateLFHeaderFile[]


CreateLFSourceFile[]


Remove["Matchete`*"]
