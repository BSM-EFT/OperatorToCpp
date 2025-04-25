(* ::Package:: *)

(* ::Chapter:: *)
(*-- Package Context --*)


BeginPackage["OperatorExport`", {"Global`"}];


(*Begin["`Private`"];*)


NotebookDirectory[]


(* ::Chapter:: *)
(*-- Function definitions --*)


(* ::Section:: *)
(*Clean up and simplification of results *)


(* ::Subsection:: *)
(*Create a dictionary mapping operator content and matched expressions*)


createMatchingDict[WarsawOutput_] := Module[{d6Dict, keyList, key1, key2},
	d6Dict= Association[];
	keyList = Keys[WarsawOutput];
	
	Do[
		key1=StringCases[ToString[keyList[[i]]], RegularExpression["(?<=Coupling\\[c)"]~~x:WordCharacter..:>"c"<>x];
		(* this assumption is crucial - Matchete produces Warsaw basis Wilson Coefficients and their names start with "c" *)
		If[Length[key1]==0,Continue[], AppendTo[d6Dict,key1[[1]]->Values[WarsawOutput][[i]]]],
	{i,Length[keyList]}];
	
	Do[
		key2=Keys[WarsawAll][[i]];
		If[KeyExistsQ[d6Dict,key2],
			Continue[],
			AppendTo[d6Dict,key2->0.0]],
		{i,Length[WarsawAll]}];
		
	Return[d6Dict];
]


(* ::Subsection:: *)
(*Simplify the output through suitable replacement of symbols*)


SimplifyOutput[output_]:=Module[{dict,newDict},
	dict = createMatchingDict[output];
	dict = dict /. {Index[i1,Flavor]->i1,Index[i2,Flavor]->i2,Index[i3,Flavor]->i3,Index[i4,Flavor]->i4 };
	dict = dict /. {Index[x_,__]:>ToExpression[StringJoin["r",StringPart[ToString[x],-1]]]};
	dict = dict /. Bar[Coupling[x_,{y___},z_]] :> Coupling[ToExpression[ToString[x]<>"c"],{y},z];
	dict = dict /. {Coupling[x_,{},__]:>x} /.{Coupling[x_,{a_},__]:>Mass[x,a]}/.{Coupling[x_,{p_,q_},__]:>TwoDim[x,p,q](*Which[StringContainsQ[ToString[x],"cY"],YukawaFn[x,p,q],True,TwoDim[x,p,q]]*)};
	dict = dict /. {FlavorSum[x_]:>1};
		
	newDict=Association[];
	Do[AppendTo[newDict,Keys[dict][[i]]->Values[dict][[i]]],{i,1,Length[dict]}];
	Return[newDict];
]


(* ::Section:: *)
(*Extraction of parameters *)


(* ::Subsection:: *)
(*Replacement list for complex parameters*)


ReplaceCCVars[list_] := Module[{rules},
	rules={};
	If[
		Length[list]==0,
		Return[rules],
		Do[AppendTo[rules,ToExpression[ToString[list[[i]]]<>"c"]->list[[i]]],{i,1,Length[list]}];
		Return[rules]
	]
]

ReplaceCCString[list_] := Module[{rules},
	rules={};
	If[
		Length[list]==0,
		Return[rules],
		Do[AppendTo[rules,ToString[list[[i]]]<>"c"->ToString[list[[i]]]],{i,1,Length[list]}];
		Return[rules]
	]
]


(* ::Subsection:: *)
(*Extract parameters from the matched expressions*)


ExtractParams[expression_,ComplexPars_] := Module[{paramList},
	paramList = Variables[expression];
	paramList = paramList /. LF[x__,__]:>x;
	paramList = paramList /. TwoDim[x_,__] :> Mat[x];
	(*paramList = paramList /. YukawaFn[x_,__] :> YF[x];*)
	paramList = paramList /. Mass[x_,_] :> Vec[x];
	paramList = paramList /. Log[x_] :> (ToExpression["mubarsq"]/x);
	paramList = paramList /. y__^2 :> y;
	paramList = paramList /. ReplaceCCVars[ComplexPars];
	paramList = paramList /. {Delta[__] :> 1, (*cH2 -> 1, cG2 -> 1, cB2 -> 1, cW2 -> 1,*) hbar->1};
	Return[DeleteCases[DeleteDuplicates[Flatten[paramList]],_Integer]];
]

ExtractAllParams[matchingDict_,ComplexPars_] := Module[{paramList,expr},
	paramList = List[];
	Do[
		expr = Values[matchingDict][[i]];
		AppendTo[paramList,ExtractParams[expr,ComplexPars]], {i,1,Length[matchingDict]}
	];

	Return[Sort[DeleteDuplicates[Flatten[paramList]]]];
]

SegregateParams[matchingDict_,ComplexPars_]:=Module[{extractedList},
	extractedList=ExtractAllParams[matchingDict,ComplexPars];
	Return[{
		Select[extractedList,Head[#]===Symbol&], 
		Select[extractedList,Head[#]===Vec&]/.Vec[x_]:>x,
		Select[extractedList,Head[#]===Mat&]/.Mat[x_]:>x(*,
		Select[extractedList,Head[#]===YF&]/.YF[x_]:>x*)}];
]


(* ::Section:: *)
(*Rewrite the matched expression in a C++ friendly form*)


(* ::Subsection:: *)
(*Create wrappers for functions of vector masses and matrix products*)


(* ::Subsubsection:: *)
(*Replace the loop function exponents*)


ReplLoopFuncExpns = {
{1,0}->1,{1,1}->2,{1,2}->3,{2,-1}->4,{2,0}->5,{2,1}->6,{3,-1}->7,{3,0}->8,{4,-2}->9,{4,-1}->10,
{4,0}->11,{5,-2}->12,{5,-1}->13,{5,0}->14,{6,-3}->15,{6,-2}->16,{6,-1}->17,{6,0}->18,{1,1,-1}->19,{1,1,0}->20,
{1,2,-1}->21,{1,2,0}->22,{1,3,-1}->23,{1,3,0}->24,{1,4,-1}->25,{1,4,0}->26,{1,5,0}->27,{2,1,-1}->28,{2,1,0}->29,{2,2,-2}->30,
{2,2,-1}->31,{2,2,0}->32,{2,3,-2}->33,{2,3,-1}->34,{2,3,0}->35,{2,4,-1}->36,{2,4,0}->37,{3,1,-2}->38,{3,1,-1}->39,{3,1,0}->40,
{3,2,-2}->41,{3,2,-1}->42,{3,2,0}->43,{3,3,-3}->44,{3,3,-2}->45,{3,3,-1}->46,{3,3,0}->47,{4,1,-2}->48,{4,1,-1}->49,{4,1,0}->50,
{4,2,-3}->51,{4,2,-2}->52,{4,2,-1}->53,{4,2,0}->54,{5,1,-3}->55,{5,1,-2}->56,{5,1,-1}->57,{5,1,0}->58,{1,1,1,-1}->59,{1,1,1,0}->60,
{1,1,2,-1}->61,{1,1,2,0}->62,{1,1,3,-1}->63,{1,1,3,0}->64,{1,1,4,0}->65,{1,2,1,-1}->66,{1,2,1,0}->67,{1,2,2,-1}->68,{1,2,2,0}->69,{1,2,3,0}->70,
{1,3,1,-1}->71,{1,3,1,0}->72,{1,3,2,0}->73,{1,4,1,0}->74,{2,1,1,-2}->75,{2,1,1,-1}->76,{2,1,1,0}->77,{2,1,2,-2}->78,{2,1,2,-1}->79,{2,1,2,0}->80,
{2,1,3,0}->81,{2,2,1,-2}->82,{2,2,1,-1}->83,{2,2,1,0}->84,{2,2,2,-1}->85,{2,2,2,0}->86,{2,3,1,-1}->87,{2,3,1,0}->88,{3,1,1,-2}->89,{3,1,1,-1}->90,
{3,1,1,0}->91,{3,1,2,-1}->92,{3,1,2,0}->93,{3,2,1,-3}->94,{3,2,1,-2}->95,{3,2,1,-1}->96,{3,2,1,0}->97,{4,1,1,-3}->98,{4,1,1,-2}->99,{4,1,1,-1}->100,
{4,1,1,0}->101,{1,1,1,1,-1}->102,{1,1,1,1,0}->103,{1,1,1,2,-1}->104,{1,1,1,2,0}->105,{1,1,1,3,0}->106,{1,1,2,1,-1}->107,{1,1,2,1,0}->108,{1,1,2,2,0}->109,{1,1,3,1,0}->110,
{1,2,1,1,-1}->111,{1,2,1,1,0}->112,{1,2,1,2,0}->113,{1,2,2,1,0}->114,{1,3,1,1,0}->115,{2,1,1,1,-2}->116,{2,1,1,1,-1}->117,{2,1,1,1,0}->118,{2,1,1,2,0}->119,{2,1,2,1,0}->120,
{2,2,1,1,-1}->121,{2,2,1,1,0}->122,{3,1,1,1,-1}->123,{3,1,1,1,0}->124,{1,1,1,1,1,-1}->125,{1,1,1,1,1,0}->126,{1,1,1,1,2,0}->127,{1,1,1,2,1,0}->128,{1,1,2,1,1,0}->129,{1,2,1,1,1,0}->130,
{2,1,1,1,1,0}->131,{1,1,1,1,1,1,0}->132
};

(* automate the creation of this list of rules as well *)


(* ::Subsubsection:: *)
(*Wrapper for loop functions + extraction of flavour indices*)


WrapLFinTF=LF[arg1_,arg2_]:>Module[{argsList,vecMasses,vecMassindices},
	argsList=List@@arg1;
	vecMasses=Select[argsList,MatchQ[#,_[_,__]]&];
	If[Length[vecMasses]==0,
		TF[LoopFunc[arg1,arg2,ToExpression["mubarsq"]],List[]],
		vecMassindices = Flatten[vecMasses[[All,2]]];
		TF[LoopFunc[arg1,arg2,ToExpression["mubarsq"]],List[Sequence@@vecMassindices]]
	]
];


(* ::Subsubsection:: *)
(*Wrapper for powers of masses + extraction of flavour indices*)


WrapMassPowinTF={Mass[x_,y_]^n_:>TF[MassPow[Mass[x,y],n],{y}],Mass[x_,y_]:>TF[MassPow[Mass[x,y],1],{y}]};


(* ::Subsubsection:: *)
(*Replace the dummy Mass[x,y] function call by simply the symbol*)


(* it's important to be careful of the fact that this substitution happens only inside the TF wrapper and not in the rest of the expression *)

ReplMassVecinTF[expr_]:=Module[{rest,rel},
	rest=DeleteCases[expr,_?(Head[#]===TF&)];
	rel=Select[List@@expr,Head[#]===TF&]/.{Mass[x_,_]:>x};
	If[Length[rel]==0,
		Return[rest],
		Return[rest*First[rel]]
	]
];


(* ::Subsubsection:: *)
(*Further extraction of free indices*)


ExtractFIdxfromTF =TF[arg1_,arg2_]:>Module[{freeIdxList},
	freeIdxList=Select[arg2,StringMatchQ[ToString[#],"i"~~___]&];
	TF[arg1,arg2,freeIdxList]
];


(* ::Subsubsection:: *)
(*Combining individual functions*)


LoopFuncWrapper[expr_]:=Module[{},
	If[Head[expr]===Symbol,Return[expr]];
	ReplMassVecinTF[expr/.ReplLoopFuncExpns/.WrapLFinTF]/.ExtractFIdxfromTF
];

MassPowWrapper[expr_]:=Module[{},
	If[Head[expr]===Symbol,Return[expr]];
	ReplMassVecinTF[expr/.WrapMassPowinTF]/.ExtractFIdxfromTF
];


(* ::Subsubsection:: *)
(*Wrapper for matrix products*)


MatProdWrapper[expr_]:=Module[{nonMat,lst,symbList,idxList,freeIdxList},
	If[Head[expr]===Symbol,Return[expr]];
	nonMat=DeleteCases[expr,_?((Head[#]===TwoDim || Head[#]===YukawaFn)&)];
	(*lst=Join[
		Select[List@@expr,Head[#]===TwoDim&]/.{TwoDim[x_,y___]:>{x,{y}}},
		Select[List@@expr,Head[#]===YukawaFn&]/.{YukawaFn[x_,y___]:>{MakePair[x,ToExpression["mubarsq"]],{y}}}
	];*)
	lst=Select[List@@expr,Head[#]===TwoDim&]/.{TwoDim[x_,y___]:>{x,{y}}};
	If[Length[lst]==0,
		Return[nonMat],
		symbList=lst[[;;,1]];
		idxList=lst[[;;,2]];
		freeIdxList=Select[Flatten[idxList],StringMatchQ[ToString[#],"i"~~___]&];
		Return[nonMat*TM[symbList,idxList,freeIdxList]]
	]
];



(* ::Subsection:: *)
(*Wrap all indical objects within a call to EinsSum[]*)


(* ::Subsubsection:: *)
(*Replace the index variables by integers*)


ReplIdxbyNum=Join[
	Table[ToExpression["r"<>ToString[j]]->j,{j,1,6}], 
	Table[ToExpression["i"<>ToString[j]]->10+j,{j,1,6}]
];


(* ::Subsubsection:: *)
(*Combine TF[] and TM[] wrappers and their index sequences within EinsSum[]*)


EinsSumCall[expr_]:=Module[{nonWrapped,lstTF,lstTM},
	If[Head[expr]===Symbol,Return[expr]];
	nonWrapped=DeleteCases[expr,_?((Head[#]===TF||Head[#]===TM)&)];
	lstTF=Select[List@@expr,Head[#]===TF&]/.TF[x_,y_,z_]:>{x,y,z};
	lstTM=Select[List@@expr,Head[#]===TM&]/.TM[a_,b_,c_]:>{a,b,c};
	Which[
		Length[lstTF]!=0&&Length[lstTM]!=0,
		Return[nonWrapped*EinsSum[Flatten[{lstTF[[;;,1]],First[lstTM][[1]]}],(Join[lstTF[[;;,2]],First[lstTM][[2]]]/.ReplIdxbyNum),Join[Flatten[lstTF[[;;,3]]],First[lstTM][[3]] ]]],
		Length[lstTF]!=0&&Length[lstTM]==0,
		Return[nonWrapped*EinsSum[Flatten[{lstTF[[;;,1]]}],(lstTF[[;;,2]]/.ReplIdxbyNum),Flatten[lstTF[[;;,3]]] ]],
		Length[lstTF]==0&&Length[lstTM]!=0,
		Return[nonWrapped*EinsSum[Flatten[{First[lstTM][[1]]}],(First[lstTM][[2]]/.ReplIdxbyNum),First[lstTM][[3]] ]],
		Length[lstTF]==0&&Length[lstTM]==0,
		Return[nonWrapped]
	];
];


(* ::Subsection:: *)
(*Convert expressions to string (for writing to files)*)


CreateString[expr_]:=Module[{nonTensor,Tensor,str},
	If[Head[expr]===Symbol,Return[ToString[expr]]];
	nonTensor=ToString[CForm[DeleteCases[expr,_?(Head[#]===EinsSum&)]]];
	Tensor=Select[List@@expr,Head[#]===EinsSum&];
	If[Length[Tensor]==0,
		Return[nonTensor],
		str=ToString[First[Tensor]];
		str=StringReplace[str,{"["->"(","]"->")"(*,"MakePair"->"std::make_pair"*)}];
		Return[StringJoin[{nonTensor,"*",str }]]
	];
];


ConvertSingleTerm[expr_,ComplexPars_]:=Module[{str,newStr},
	str=CreateString[
			EinsSumCall[
				MatProdWrapper[
					MassPowWrapper[
						LoopFuncWrapper[expr]]]]];
	newStr=StringReplace[str,{(*"cH2"->"cH2(mubarsq)","cG2"->"cG2(mubarsq)","cB2"->"cB2(mubarsq)","cW2"->"cW2(mubarsq)",*)"Power" -> "pow", "Log"->"log", "Sqrt"->"sqrt", "Delta("->"KronDelta(", "Pi"->"3.14159", "\""->"",".*"->"*", ". "->".0 ", ".)"->".0)" , "+ -" -> "- ",Whitespace->""}];
	newStr=StringReplace[newStr,ReplaceCCString[ComplexPars]];
	Return[newStr]
];


ConvertFullExpression[sumOfExprs_,ComplexPars_]:=Module[{exprList,str},
	Which[
		Head[sumOfExprs]===Times,
		Return[ConvertSingleTerm[sumOfExprs,ComplexPars]],
		Head[sumOfExprs]===Plus,
		exprList = List @@ sumOfExprs;
		str=ConvertSingleTerm[exprList[[1]],ComplexPars];
		Do[str=str<>" + "<>ConvertSingleTerm[exprList[[i]],ComplexPars],{i,2,Length[exprList]}];
		str=StringReplace[str,{"+ -"->"- "}];
		Return[str],
		Head[sumOfExprs]===Symbol,
		Return[ToString[sumOfExprs]],
		True,
		Return[ToString[0]]
	]
];


(* ::Section:: *)
(*Build C++ functions and write them to file*)


(* ::Subsection:: *)
(*Build the header (.h) file*)


(* ::Subsubsection:: *)
(*Builder for preprocessor directives, model class with data members, constructor and updater declarations*)


HeaderPreprocessorDirectives[line_] := Module[{},
	WriteLine[line, ""];
	WriteLine[line, "#include <vector>"];
	WriteLine[line, "#include <string>"];
	WriteLine[line, "#include <map>"];
	WriteLine[line, "#include <functional>"];
	WriteLine[line, "#define hbar 0.006332574"]
];

HeaderModelClass[className_,paramList_,line_] := Module[{args},
	WriteLine[line, ""];
	WriteLine[line, "class "<>className<>" {"];
	WriteLine[line, "    private:"];
	
	(* define and initialize 0-dimensional parameters *)
	If[Length[paramList[[1]]]!=0,
		Do[WriteLine[line, "        "<>"double "<>ToString[paramList[[1]][[i]]]<>" = 0.0;"],{i,1,Length[paramList[[1]]]}]
	];
	WriteLine[line, ""];
	
	(* define and initialize 1-dimensional parameters (masses) *)
	If[Length[paramList[[2]]]!=0,
		Do[WriteLine[line, "        "<>"std::vector<double> "<>ToString[paramList[[2]][[i]]]<>" = {0.0, 0.0, 0.0};"], {i,1,Length[paramList[[2]]]}];
	];
	WriteLine[line, ""];
	
	(* define and initialize 2-dimensional parameters *)
	If[Length[paramList[[3]]]!=0,
		Do[WriteLine[line, "        "<>"std::vector<std::vector<double> > "<>ToString[paramList[[3]][[i]]]<>" = {{0.0, 0.0, 0.0},{0.0, 0.0, 0.0},{0.0, 0.0, 0.0}};"],
		{i,1,Length[paramList[[3]]]}];
	];
	WriteLine[line, ""];
		
	WriteLine[line, "    public:"];
	(* declaration for the default constructor*)
	WriteLine[line, "        "<>className<>"() = default;"];
	WriteLine[line, ""];
	
	(* declaration for the overloaded constructor*)
	WriteLine[line, StringJoin["        ", className, "(std::map<std::string, double> params);"]];
	WriteLine[line, ""];

	(* declaration for the updater method *)
	WriteLine[line, StringJoin["        void updateParams", "(std::map<std::string, double> params);"]];
	WriteLine[line, ""];
	
	(* declaration for a method that prints the names of all parameters *)
	WriteLine[line, StringJoin["        void printParamNames", "();"]];
	WriteLine[line, ""];
	
	(* declaration for a method that prints the names and current values of all parameters *)
	WriteLine[line, StringJoin["        void printParams", "();"]];
	WriteLine[line, ""];
	
	(* declarations for WC methods (Warsaw) *)
	Do[WriteLine[line,"        double "<>Values[WarsawAll][[i]]<>";"],{i,1,Length[WarsawAll]}];
	
	WriteLine[line, ""];
	
	(* declaration and definition of std::function lambdas that bind to the Yukawa member functions of the model class *)
	(*If[Length[paramList[[4]]]!=0,
		Do[WriteLine[line, "        "<>"std::function<double(int, int, double)> "<>ToString[paramList[[4]][[i]]]<>"_ = [this](int i, int j, double m){ return this->"<>ToString[paramList[[4]][[i]]]<>"(i, j, m); };"],
		{i,1,Length[paramList[[4]]]}];
	];*)
	WriteLine[line, ""];
	
	WriteLine[line, "};"];
];


(* ::Subsubsection:: *)
(*Master builder*)


HeaderFileBuilder[modelName_,paramList_]:=Module[{path,line1},
	path = FileNameJoin[{NotebookDirectory[],"include"}];
	line1 = OpenWrite[path<>"/"<>modelName<>".h"];
	HeaderPreprocessorDirectives[line1];
	WriteLine[line1,""];
	HeaderModelClass[modelName,paramList,line1];
	Close[line1];
];



(* ::Subsection:: *)
(*Build the source (.cpp) file*)


(* ::Subsubsection:: *)
(*Builder for preprocessor directives, constructor and updater methods*)


BuildPreprocessorDirectives[modelName_,line_] := Module[{},
	WriteLine[line, ""];
	WriteLine[line, "#include <iostream>"];
	WriteLine[line, "#include <cmath>"];
	WriteLine[line, "#include <vector>"];
	WriteLine[line, "#include <string>"];
	WriteLine[line, "#include <map>"];
	WriteLine[line, "#include <tuple>"];
	WriteLine[line, "#include \"OperatorImport.h\""];
	WriteLine[line, "#include \""<>modelName<>".h\""];
];

BuildConstructor[className_, paramList_, line_]:=Module[{args},
	WriteLine[line, ""];
	WriteLine[line, StringJoin[className, "::", className, "(std::map<std::string, double> params) {"]];
	
	If[Length[paramList[[1]]]!=0,
		Do[WriteLine[line, "    if (params.contains(\""<>ToString[paramList[[1]][[i]]]<>"\")) this->"<>ToString[paramList[[1]][[i]]]<>" = params[\""<>ToString[paramList[[1]][[i]]]<>"\"];"], 
			{i,1,Length[paramList[[1]]]}]	
	];
	
	If[Length[paramList[[2]]]!=0,
		Do[
			WriteLine[line, "    if (params.contains(\""<>ToString[paramList[[2]][[i]]]<>"1\")) this->"<>ToString[paramList[[2]][[i]]]<>"[0] = params[\""<>ToString[paramList[[2]][[i]]]<>"1\"];"];
			WriteLine[line, "    if (params.contains(\""<>ToString[paramList[[2]][[i]]]<>"2\")) this->"<>ToString[paramList[[2]][[i]]]<>"[1] = params[\""<>ToString[paramList[[2]][[i]]]<>"2\"];"];
			WriteLine[line, "    if (params.contains(\""<>ToString[paramList[[2]][[i]]]<>"3\")) this->"<>ToString[paramList[[2]][[i]]]<>"[2] = params[\""<>ToString[paramList[[2]][[i]]]<>"3\"];"], 
		{i,1,Length[paramList[[2]]]}]	
	];
	
	If[Length[paramList[[3]]]!=0,
		Do[
			WriteLine[line, "    if (params.contains(\""<>ToString[paramList[[3]][[i]]]<>"11\")) this->"<>ToString[paramList[[3]][[i]]]<>"[0][0] = params[\""<>ToString[paramList[[3]][[i]]]<>"11\"];"];
			WriteLine[line, "    if (params.contains(\""<>ToString[paramList[[3]][[i]]]<>"12\")) this->"<>ToString[paramList[[3]][[i]]]<>"[0][1] = params[\""<>ToString[paramList[[3]][[i]]]<>"12\"];"];
			WriteLine[line, "    if (params.contains(\""<>ToString[paramList[[3]][[i]]]<>"13\")) this->"<>ToString[paramList[[3]][[i]]]<>"[0][2] = params[\""<>ToString[paramList[[3]][[i]]]<>"13\"];"];
			WriteLine[line, "    if (params.contains(\""<>ToString[paramList[[3]][[i]]]<>"21\")) this->"<>ToString[paramList[[3]][[i]]]<>"[1][0] = params[\""<>ToString[paramList[[3]][[i]]]<>"21\"];"];
			WriteLine[line, "    if (params.contains(\""<>ToString[paramList[[3]][[i]]]<>"22\")) this->"<>ToString[paramList[[3]][[i]]]<>"[1][1] = params[\""<>ToString[paramList[[3]][[i]]]<>"22\"];"];
			WriteLine[line, "    if (params.contains(\""<>ToString[paramList[[3]][[i]]]<>"23\")) this->"<>ToString[paramList[[3]][[i]]]<>"[1][2] = params[\""<>ToString[paramList[[3]][[i]]]<>"23\"];"];
			WriteLine[line, "    if (params.contains(\""<>ToString[paramList[[3]][[i]]]<>"31\")) this->"<>ToString[paramList[[3]][[i]]]<>"[2][0] = params[\""<>ToString[paramList[[3]][[i]]]<>"31\"];"];
			WriteLine[line, "    if (params.contains(\""<>ToString[paramList[[3]][[i]]]<>"32\")) this->"<>ToString[paramList[[3]][[i]]]<>"[2][1] = params[\""<>ToString[paramList[[3]][[i]]]<>"32\"];"];
			WriteLine[line, "    if (params.contains(\""<>ToString[paramList[[3]][[i]]]<>"33\")) this->"<>ToString[paramList[[3]][[i]]]<>"[2][2] = params[\""<>ToString[paramList[[3]][[i]]]<>"33\"];"], 
		{i,1,Length[paramList[[3]]]}]	
	];
	
	WriteLine[line, StringJoin["}"]];
];

BuildUpdater[className_, paramList_, line_]:=Module[{},
	WriteLine[line, ""];
	WriteLine[line, StringJoin["void ", className, "::updateParams", "(std::map<std::string, double> params) {"]];
	
	If[Length[paramList[[1]]]!=0,
		Do[WriteLine[line, "    if (params.contains(\""<>ToString[paramList[[1]][[i]]]<>"\")) this->"<>ToString[paramList[[1]][[i]]]<>" = params[\""<>ToString[paramList[[1]][[i]]]<>"\"];"], 
		{i,1,Length[paramList[[1]]]}]	
	];
	
	If[Length[paramList[[2]]]!=0,
		Do[
			WriteLine[line, "    if (params.contains(\""<>ToString[paramList[[2]][[i]]]<>"1\")) this->"<>ToString[paramList[[2]][[i]]]<>"[0] = params[\""<>ToString[paramList[[2]][[i]]]<>"1\"];"];
			WriteLine[line, "    if (params.contains(\""<>ToString[paramList[[2]][[i]]]<>"2\")) this->"<>ToString[paramList[[2]][[i]]]<>"[1] = params[\""<>ToString[paramList[[2]][[i]]]<>"2\"];"];
			WriteLine[line, "    if (params.contains(\""<>ToString[paramList[[2]][[i]]]<>"3\")) this->"<>ToString[paramList[[2]][[i]]]<>"[2] = params[\""<>ToString[paramList[[2]][[i]]]<>"3\"];"], 
		{i,1,Length[paramList[[2]]]}]	
	];
	
	If[Length[paramList[[3]]]!=0,
		Do[
			WriteLine[line, "    if (params.contains(\""<>ToString[paramList[[3]][[i]]]<>"11\")) this->"<>ToString[paramList[[3]][[i]]]<>"[0][0] = params[\""<>ToString[paramList[[3]][[i]]]<>"11\"];"];
			WriteLine[line, "    if (params.contains(\""<>ToString[paramList[[3]][[i]]]<>"12\")) this->"<>ToString[paramList[[3]][[i]]]<>"[0][1] = params[\""<>ToString[paramList[[3]][[i]]]<>"12\"];"];
			WriteLine[line, "    if (params.contains(\""<>ToString[paramList[[3]][[i]]]<>"13\")) this->"<>ToString[paramList[[3]][[i]]]<>"[0][2] = params[\""<>ToString[paramList[[3]][[i]]]<>"13\"];"];
			WriteLine[line, "    if (params.contains(\""<>ToString[paramList[[3]][[i]]]<>"21\")) this->"<>ToString[paramList[[3]][[i]]]<>"[1][0] = params[\""<>ToString[paramList[[3]][[i]]]<>"21\"];"];
			WriteLine[line, "    if (params.contains(\""<>ToString[paramList[[3]][[i]]]<>"22\")) this->"<>ToString[paramList[[3]][[i]]]<>"[1][1] = params[\""<>ToString[paramList[[3]][[i]]]<>"22\"];"];
			WriteLine[line, "    if (params.contains(\""<>ToString[paramList[[3]][[i]]]<>"23\")) this->"<>ToString[paramList[[3]][[i]]]<>"[1][2] = params[\""<>ToString[paramList[[3]][[i]]]<>"23\"];"];
			WriteLine[line, "    if (params.contains(\""<>ToString[paramList[[3]][[i]]]<>"31\")) this->"<>ToString[paramList[[3]][[i]]]<>"[2][0] = params[\""<>ToString[paramList[[3]][[i]]]<>"31\"];"];
			WriteLine[line, "    if (params.contains(\""<>ToString[paramList[[3]][[i]]]<>"32\")) this->"<>ToString[paramList[[3]][[i]]]<>"[2][1] = params[\""<>ToString[paramList[[3]][[i]]]<>"32\"];"];
			WriteLine[line, "    if (params.contains(\""<>ToString[paramList[[3]][[i]]]<>"33\")) this->"<>ToString[paramList[[3]][[i]]]<>"[2][2] = params[\""<>ToString[paramList[[3]][[i]]]<>"33\"];"], 
		{i,1,Length[paramList[[3]]]}]	
	];
	
	WriteLine[line, "}"];
];


(* ::Subsubsection::Closed:: *)
(*Builder for printer methods*)


BuildPrinter[className_, paramList_, line_]:=Module[{},
	WriteLine[line, ""];
	WriteLine[line, "void "<>className<>"::printParamNames"<>"()"<>"{"];
	
	If[Length[paramList[[1]]]!=0,
		WriteLine[line, StringJoin["    std::cout << ","\"\\n", "Scalar-valued parameters", "\"", " << \"\\n\";"]];
		WriteLine[line, StringJoin["    std::cout << ","\"", "-----------------", "\"", " << \"\\n\";"]];
		Do[WriteLine[line, StringJoin["    std::cout << ", "\"", ToString[paramList[[1]][[i]]], "\" << \"\\n\";"]],{i,1,Length[paramList[[1]]]}];
		WriteLine[line, ""]
	];
	
	If[Length[paramList[[2]]]!=0,
		WriteLine[line, StringJoin["    std::cout << ","\"\\n", "Vector-valued parameters", "\"", " << \"\\n\";"]];
		WriteLine[line, StringJoin["    std::cout << ","\"", "-----------------", "\"", " << \"\\n\";"]];
		Do[WriteLine[line, StringJoin["    std::cout << ", "\"", ToString[paramList[[2]][[i]]], "\" << \"\\n\";"]],{i,1,Length[paramList[[2]]]}];
		WriteLine[line, ""]
	];
	
	If[Length[paramList[[3]]]!=0,
		WriteLine[line, StringJoin["    std::cout << ","\"\\n", "Matrix-valued parameters", "\"", " << \"\\n\";"]];
		WriteLine[line, StringJoin["    std::cout << ","\"", "-----------------", "\"", " << \"\\n\";"]];
		Do[WriteLine[line, StringJoin["    std::cout << ", "\"", ToString[paramList[[3]][[i]]], "\" << \"\\n\";"]],{i,1,Length[paramList[[3]]]}];
		WriteLine[line, ""]
	];
	
	WriteLine[line, "}"];
	
	WriteLine[line, ""];
	WriteLine[line, "void "<>className<>"::printParams"<>"()"<>"{"];
	
	If[Length[paramList[[1]]]!=0,
		WriteLine[line, StringJoin["    std::cout << ","\"\\n", "Scalar-valued parameters", "\"", " << \"\\n\";"]];
		WriteLine[line, StringJoin["    std::cout << ","\"", "-----------------", "\"", " << \"\\n\";"]];
		Do[
			WriteLine[line, StringJoin["    std::cout << ", "\"", ToString[paramList[[1]][[i]]], ": \" << ", ToString[paramList[[1]][[i]]], " << \"\\n\";"]],
		{i,1,Length[paramList[[1]]]}];
		WriteLine[line, ""]
	];
	
	If[Length[paramList[[2]]]!=0,
		WriteLine[line, StringJoin["    std::cout << ","\"\\n", "Vector-valued parameters", "\"", " << \"\\n\";"]];
		WriteLine[line, StringJoin["    std::cout << ","\"", "-----------------", "\"", " << \"\\n\";"]];
		Do[
			WriteLine[line, StringJoin["    std::cout << ", "\"", ToString[paramList[[2]][[i]]], "1: \" << ", ToString[paramList[[2]][[i]]], "[0] << \"\\n\";"]];
			WriteLine[line, StringJoin["    std::cout << ", "\"", ToString[paramList[[2]][[i]]], "2: \" << ", ToString[paramList[[2]][[i]]], "[1] << \"\\n\";"]];
			WriteLine[line, StringJoin["    std::cout << ", "\"", ToString[paramList[[2]][[i]]], "3: \" << ", ToString[paramList[[2]][[i]]], "[2] << \"\\n\";"]],
		{i,1,Length[paramList[[2]]]}];
		WriteLine[line, ""]
	];
	
	If[Length[paramList[[3]]]!=0,
		WriteLine[line, StringJoin["    std::cout << ","\"\\n", "Matrix-valued parameters", "\"", " << \"\\n\";"]];
		WriteLine[line, StringJoin["    std::cout << ","\"", "-----------------", "\"", " << \"\\n\";"]];
		Do[
			WriteLine[line, StringJoin["    std::cout << ", "\"", ToString[paramList[[3]][[i]]], "11: \" << ", ToString[paramList[[3]][[i]]], "[0][0] << \"\\n\";"]];
			WriteLine[line, StringJoin["    std::cout << ", "\"", ToString[paramList[[3]][[i]]], "12: \" << ", ToString[paramList[[3]][[i]]], "[0][1] << \"\\n\";"]];
			WriteLine[line, StringJoin["    std::cout << ", "\"", ToString[paramList[[3]][[i]]], "13: \" << ", ToString[paramList[[3]][[i]]], "[0][2] << \"\\n\";"]];
			WriteLine[line, StringJoin["    std::cout << ", "\"", ToString[paramList[[3]][[i]]], "21: \" << ", ToString[paramList[[3]][[i]]], "[1][0] << \"\\n\";"]];
			WriteLine[line, StringJoin["    std::cout << ", "\"", ToString[paramList[[3]][[i]]], "22: \" << ", ToString[paramList[[3]][[i]]], "[1][1] << \"\\n\";"]];
			WriteLine[line, StringJoin["    std::cout << ", "\"", ToString[paramList[[3]][[i]]], "23: \" << ", ToString[paramList[[3]][[i]]], "[1][2] << \"\\n\";"]];
			WriteLine[line, StringJoin["    std::cout << ", "\"", ToString[paramList[[3]][[i]]], "31: \" << ", ToString[paramList[[3]][[i]]], "[2][0] << \"\\n\";"]];
			WriteLine[line, StringJoin["    std::cout << ", "\"", ToString[paramList[[3]][[i]]], "32: \" << ", ToString[paramList[[3]][[i]]], "[2][1] << \"\\n\";"]];
			WriteLine[line, StringJoin["    std::cout << ", "\"", ToString[paramList[[3]][[i]]], "33: \" << ", ToString[paramList[[3]][[i]]], "[2][2] << \"\\n\";"]],
		{i,1,Length[paramList[[3]]]}];
		WriteLine[line, ""]
	];
	
	WriteLine[line, "}"];
];


(* ::Subsubsection:: *)
(*Builder for a single WC function (Warsaw Basis)*)


BuildFunctionWarsaw[modelName_,WCname_,expr_,ComplexPars_,(*YFReplRule_,*)line_]:=Module[{returnExpr},
	returnExpr = ConvertFullExpression[expr,ComplexPars];
	(*If[
		StringContainsQ[WCname,"cY"],
		ConvertFullExpression[expr,ComplexPars],
		StringReplace[ConvertFullExpression[expr,ComplexPars],YFReplRule]
	];*)
	
	WriteLine[line,""];
	WriteLine[line, "double "<>modelName<>"::"<>WCname<> " {"];		
	WriteLine[line, "    return ("<>returnExpr<>");"];
	WriteLine[line, "}"];
]


(* ::Subsubsection:: *)
(*Builder for a single SMEFiT WC function*)


BuildFunctionSMEFiT[modelName_,WCname_,expr_,line_]:=Module[{},
	WriteLine[line," "];
	WriteLine[line, "double "<>modelName<>"::"<>WCname<> "(double mubarsq) {"];			
	WriteLine[line, "    return "<>expr<>";"];
	WriteLine[line, "}"];
]



(* ::Subsubsection:: *)
(*Master builder *)


ReplaceVarName[list_,str_] := Module[{rules},
	rules={};
	If[
		Length[list]==0,
		Return[rules],
		Do[AppendTo[rules,ToString[list[[i]]]->ToString[list[[i]]]<>str],{i,1,Length[list]}];
		Return[rules]
	]
];


SourceFileBuilder[modelName_, paramList_, ComplexPars_, matchingOutput_]:=Module[{keyList, exprList, path, (*YFReplRule,*) line1},
	keyList=Keys[matchingOutput];
	exprList=Values[matchingOutput];
	
	path = FileNameJoin[{NotebookDirectory[],"lib"}];
	line1 = OpenWrite[path<>"/"<>modelName<>".cpp"];
	
	BuildPreprocessorDirectives[modelName,line1];
	BuildConstructor[modelName,paramList,line1];
	BuildUpdater[modelName,paramList,line1];
	BuildPrinter[modelName,paramList,line1];
	
	(*YFReplRule = ReplaceVarName[paramList[[4]],"_"];*)
	
	Do[
		BuildFunctionWarsaw[modelName,WarsawAll[keyList[[k]]],exprList[[k]],ComplexPars(*,YFReplRule*),line1],
	{k,1,Length[matchingOutput]}];
	
	Close[line1];
]



(* ::Subsection:: *)
(*Create .cpp and .h files*)


CreateSourceAndHeader[modelName_,paramList_,WarsawOutput_]:=Module[{line1},
	HeaderFileBuilder[modelName,paramList];
	SourceFileBuilder[modelName, paramList, WarsawOutput]
];


(* ::Section:: *)
(*Dictionaries (WCs and Miscellaneous) *)


(* ::Subsection:: *)
(*A dictionary containing all Warsaw basis operator names and the corresponding C++ function prototype*)


(* In the future, the creation of such a dictionary can also be automated *)

(* Here, we are treating all coefficients to be real*)

WarsawAll = Association[
	(* extra coefficients *)
	(*"cH2" -> "cH2(double mubarsq)",
	"cG2" -> "cG2(double mubarsq)",
	"cW2" -> "cW2(double mubarsq)",
	"cB2" -> "cB2(double mubarsq)",
	
	"cYle" -> "cYle(int i1, int i2, double mubarsq)",
	"cYqu" -> "cYqu(int i1, int i2, double mubarsq)",
	"cYqd" -> "cYqd(int i1, int i2, double mubarsq)",*)
	
	(* operator coefficients *)
	"cllHH" -> "cllHH(int i1, int i2, double mubarsq)",
	
	"cG" -> "cG(double mubarsq)",
	"cW" -> "cW(double mubarsq)",
	"cGt" -> "cGt(double mubarsq)",
	"cWt" -> "cWt(double mubarsq)",
	
	"cH" -> "cH(double mubarsq)",
	"cHBox" -> "cHBox(double mubarsq)",
	"cHD" -> "cHD(double mubarsq)",
	
	"cHG" -> "cHG(double mubarsq)",
	"cHW" -> "cHW(double mubarsq)",
	"cHB" -> "cHB(double mubarsq)",
	"cHWB" -> "cHWB(double mubarsq)",
	"cHGt" -> "cHGt(double mubarsq)",
	"cHWt" -> "cHWt(double mubarsq)",
	"cHBt" -> "cHBt(double mubarsq)",
	"cHWtB" -> "cHWtB(double mubarsq)",
	
	"ceH" -> "ceH(int i1, int i2, double mubarsq)",
	"cuH" -> "cuH(int i1, int i2, double mubarsq)",
	"cdH" -> "cdH(int i1, int i2, double mubarsq)",
	
	"ceW" -> "ceW(int i1, int i2, double mubarsq)",
	"ceB" -> "ceB(int i1, int i2, double mubarsq)",
	"cuG" -> "cuG(int i1, int i2, double mubarsq)",
	"cuW" -> "cuW(int i1, int i2, double mubarsq)",
	"cuB" -> "cuB(int i1, int i2, double mubarsq)",
	"cdG" -> "cdG(int i1, int i2, double mubarsq)",
	"cdW" -> "cdW(int i1, int i2, double mubarsq)",
	"cdB" -> "cdB(int i1, int i2, double mubarsq)",
	
	"cHl1" -> "cHl1(int i1, int i2, double mubarsq)",
	"cHl3" -> "cHl3(int i1, int i2, double mubarsq)",
	"cHe" -> "cHe(int i1, int i2, double mubarsq)",
	"cHq1" -> "cHq1(int i1, int i2, double mubarsq)",
	"cHq3" -> "cHq3(int i1, int i2, double mubarsq)",
	"cHu" -> "cHu(int i1, int i2, double mubarsq)",
	"cHd" -> "cHd(int i1, int i2, double mubarsq)",
	"cHud" -> "cHud(int i1, int i2, double mubarsq)",
	
	"cll" -> "cll(int i1, int i2, int i3, int i4, double mubarsq)",
	"cqq1" -> "cqq1(int i1, int i2, int i3, int i4, double mubarsq)",
	"cqq3" -> "cqq3(int i1, int i2, int i3, int i4, double mubarsq)",
	"clq1" -> "clq1(int i1, int i2, int i3, int i4, double mubarsq)",
	"clq3" -> "clq3(int i1, int i2, int i3, int i4, double mubarsq)",
	"cee" -> "cee(int i1, int i2, int i3, int i4, double mubarsq)",
	"cuu" -> "cuu(int i1, int i2, int i3, int i4, double mubarsq)",
	"cdd" -> "cdd(int i1, int i2, int i3, int i4, double mubarsq)",
	"ceu" -> "ceu(int i1, int i2, int i3, int i4, double mubarsq)",
	"ced" -> "ced(int i1, int i2, int i3, int i4, double mubarsq)",
	"cud1" -> "cud1(int i1, int i2, int i3, int i4, double mubarsq)",
	"cud8" -> "cud8(int i1, int i2, int i3, int i4, double mubarsq)",
	"cle" -> "cle(int i1, int i2, int i3, int i4, double mubarsq)",
	"clu" -> "clu(int i1, int i2, int i3, int i4, double mubarsq)",
	"cld" -> "cld(int i1, int i2, int i3, int i4, double mubarsq)",
	"cqe" -> "cqe(int i1, int i2, int i3, int i4, double mubarsq)",
	"cqu1" -> "cqu1(int i1, int i2, int i3, int i4, double mubarsq)",
	"cqu8" -> "cqu8(int i1, int i2, int i3, int i4, double mubarsq)",
	"cqd1" -> "cqd1(int i1, int i2, int i3, int i4, double mubarsq)",
	"cqd8" -> "cqd8(int i1, int i2, int i3, int i4, double mubarsq)",
	
	"cledq" -> "cledq(int i1, int i2, int i3, int i4, double mubarsq)",
	"cquqd1" -> "cquqd1(int i1, int i2, int i3, int i4, double mubarsq)",
	"cquqd8" -> "cquqd8(int i1, int i2, int i3, int i4, double mubarsq)",
	"clequ1" -> "clequ1(int i1, int i2, int i3, int i4, double mubarsq)",
	"clequ3" -> "clequ3(int i1, int i2, int i3, int i4, double mubarsq)" ,
	
	"cduq" -> "cduq(int i1, int i2, int i3, int i4, double mubarsq)",
	"cqqu" -> "cqqu(int i1, int i2, int i3, int i4, double mubarsq)",
	"cqqq" -> "cqqq(int i1, int i2, int i3, int i4, double mubarsq)",
	"cduu" -> "cduu(int i1, int i2, int i3, int i4, double mubarsq)"
];


(* ::Chapter:: *)
(*-- Package Context --*)


(*End[];*)
EndPackage[];
