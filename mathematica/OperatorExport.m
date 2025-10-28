(* ::Package:: *)

(* ::Chapter:: *)
(*-- Package Context --*)


BeginPackage["OperatorExport`", {"Global`"}];


(* ::Subsection:: *)
(*Description of public functions*)


SimplifyOutput::usage = "SimplifyOutput[matchedResult] creates a dictionary containing only Wilson coefficients and the corresponding matching conditions, while unwraping the Matchete API and converting the expressions to a lighter intermediate form.";
SegregateParams::usage = "SegregateParams[SimplifiedOutput, ComlexParams] collects all unique parameters within the simpplified matching output and segregates them based on their dimensionality. Keeps track of the complex parameters to ensure that their conjugates are not counted as distinct parameters.";
HeaderFileBuilder::usage = "HeaderFileBuilder[ModelName,ModelParams] creates a ModelName.h file that defines a ModelName class with the model parameters as member variables and Warsaw basis Wilson coefficients as methods. Also, declares methods for initalizing, updating and printing the parameters of a ModelName object.";
SourceFileBuilder::usage = "SourceFileBuilder[ModelName,ModelParams,ComplexParams,SimplifiedOutput] creates a ModelName.cpp file with implementations for each Wilson coefficient method corresponding to the ModelName class. The method bodies are built by further manipulating the SimplifiedOutput.";
GeneratePythonDeclarations::usage = "GeneratePythonDeclarations[ModelName] creates a ModelName.pyi file with declarations for the model class and each Wilson coefficient method.";


(* ::Subsection:: *)
(*External symbols*)


(* External symbols*)
Index;
Flavor;
FlavorSum;
i1; i2; i3; i4; i5; i6;
r1; r2; r3; r4; r5; r6;
Coupling;
Bar;
LF;
hbar;
Delta;
LoopFunc;
MassPow;
EinsSum;


(* ::Subsection:: *)
(*Package and author information*)


Print[
 Style[
  Column[{
    Row[{
      "OperatorExport  v0.1 \[LongDash] by ",
      "Suraj Prakash",
      " (", 
      Style["suraj.prakash@ific.uv.es", "Hyperlink", FontColor -> Blue],
      ")"
    }],
    Row[{
      "Affiliation: IFIC (Universitat de Valencia - CSIC)"
    }],
    Row[{
      "GitHub: ", 
      Style["https://github.com/BSM-EFT/OperatorToCpp", "Hyperlink", FontColor -> Blue]
    }],
    Style["Mathematica component of the OperatorToC++ code.", FontColor->Gray]
  }],
  "Text"
 ]
];



Begin["`Private`"];


NotebookDirectory[];


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
	dict = dict /. {Coupling[x_,{},__]:>x} /.{Coupling[x_,{a_},__]:>Mass[x,a]}/.{Coupling[x_,{p_,q_},__]:>TwoDim[x,p,q]};
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
	paramList = paramList /. Mass[x_,_] :> Vec[x];
	paramList = paramList /. Log[x_] :> (ToExpression["mubarsq"]/x);
	paramList = paramList /. y__^2 :> y;
	paramList = paramList /. ReplaceCCVars[ComplexPars];
	paramList = paramList /. {Delta[__] :> 1,hbar->1};
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
		Select[extractedList,Head[#]===Mat&]/.Mat[x_]:>x}];
]


(* ::Section:: *)
(*Rewrite the matched expression in a C++ friendly form*)


(* ::Subsection:: *)
(*Create wrappers for functions of vector masses and matrix products*)


(* ::Subsubsection:: *)
(*Replace the loop function exponents*)


RulesFile=FileNameJoin[{NotebookDirectory[],"./LFRules.m"}];
ReplLoopFuncExpns=Import[RulesFile];


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
	nonMat=DeleteCases[expr,_?(Head[#]===TwoDim&)];
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
		str=StringReplace[str,{"["->"(","]"->")"}];
		Return[StringJoin[{nonTensor,"*",str }]]
	];
];


ConvertSingleTerm[expr_,ComplexPars_]:=Module[{str,newStr},
	str=CreateString[
			EinsSumCall[
				MatProdWrapper[
					MassPowWrapper[
						LoopFuncWrapper[expr]]]]];
	newStr=StringReplace[str,{"Power" -> "pow", "Log"->"log", "Sqrt"->"sqrt", "Delta("->"KronDelta(", "Pi"->"3.14159", "\""->"",".*"->"*", ". "->".0 ", ".)"->".0)" , "+ -" -> "- ",Whitespace->""}];
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
	WriteLine[line, "#pragma once"];
	WriteLine[line, "#include <vector>"];
	WriteLine[line, "#include <string>"];
	WriteLine[line, "#include <unordered_map>"];
];

HeaderModelClass[className_,paramList_,line_] := Module[{args},
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
	WriteLine[line, StringJoin["        ", className, "(std::unordered_map<std::string, double> params);"]];
	WriteLine[line, ""];

	(* declaration for the updater method *)
	WriteLine[line, StringJoin["        void updateParams", "(std::unordered_map<std::string, double> params);"]];
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
	
	WriteLine[line, "};"];
	WriteLine[line, ""];
];


(* ::Subsubsection:: *)
(*Master builder*)


HeaderFileBuilder[modelName_,paramList_]:=Module[{path,line1,line2},
	path = FileNameJoin[{ParentDirectory[NotebookDirectory[]],"include"}];
	
	line1 = OpenWrite[path<>"/"<>modelName<>".h"];
	HeaderPreprocessorDirectives[line1];
	WriteLine[line1,""];
	HeaderModelClass[modelName,paramList,line1];
	Close[line1];
	
	line2 = OpenWrite[path<>"/"<>"modelName.h"];
	WriteLine[line2, "#pragma once"];
	WriteLine[line2, "#include \""<>modelName<>".h\""];
	WriteLine[line2, "using Model = "<>modelName<>";"];
	WriteLine[line2, "const char py_class[] = \""<>modelName<>"\";"];
	Close[line2];	
];



(* ::Subsection:: *)
(*Build the source (.cpp) file*)


(* ::Subsubsection:: *)
(*Builder for preprocessor directives, constructor and updater methods*)


BuildPreprocessorDirectives[modelName_,line_] := Module[{},
	WriteLine[line, "#include \"pch.h\""];
	WriteLine[line, "#include \"OperatorImport.h\""];
	WriteLine[line, "#include \""<>modelName<>".h\""];
];

BuildConstructor[className_, paramList_, line_]:=Module[{args},
	WriteLine[line, ""];
	WriteLine[line, StringJoin[className, "::", className, "(std::unordered_map<std::string, double> params) {"]];
	
	If[Length[paramList[[1]]]!=0,
		Do[WriteLine[line, "    if (params.find(\""<>ToString[paramList[[1]][[i]]]<>"\") != params.end()) this->"<>ToString[paramList[[1]][[i]]]<>" = params[\""<>ToString[paramList[[1]][[i]]]<>"\"];"], 
			{i,1,Length[paramList[[1]]]}]	
	];
	
	If[Length[paramList[[2]]]!=0,
		Do[
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[2]][[i]]]<>"1\") != params.end()) this->"<>ToString[paramList[[2]][[i]]]<>"[0] = params[\""<>ToString[paramList[[2]][[i]]]<>"1\"];"];
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[2]][[i]]]<>"2\") != params.end()) this->"<>ToString[paramList[[2]][[i]]]<>"[1] = params[\""<>ToString[paramList[[2]][[i]]]<>"2\"];"];
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[2]][[i]]]<>"3\") != params.end()) this->"<>ToString[paramList[[2]][[i]]]<>"[2] = params[\""<>ToString[paramList[[2]][[i]]]<>"3\"];"], 
		{i,1,Length[paramList[[2]]]}]	
	];
	
	If[Length[paramList[[3]]]!=0,
		Do[
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"11\") != params.end()) this->"<>ToString[paramList[[3]][[i]]]<>"[0][0] = params[\""<>ToString[paramList[[3]][[i]]]<>"11\"];"];
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"12\") != params.end()) this->"<>ToString[paramList[[3]][[i]]]<>"[0][1] = params[\""<>ToString[paramList[[3]][[i]]]<>"12\"];"];
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"13\") != params.end()) this->"<>ToString[paramList[[3]][[i]]]<>"[0][2] = params[\""<>ToString[paramList[[3]][[i]]]<>"13\"];"];
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"21\") != params.end()) this->"<>ToString[paramList[[3]][[i]]]<>"[1][0] = params[\""<>ToString[paramList[[3]][[i]]]<>"21\"];"];
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"22\") != params.end()) this->"<>ToString[paramList[[3]][[i]]]<>"[1][1] = params[\""<>ToString[paramList[[3]][[i]]]<>"22\"];"];
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"23\") != params.end()) this->"<>ToString[paramList[[3]][[i]]]<>"[1][2] = params[\""<>ToString[paramList[[3]][[i]]]<>"23\"];"];
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"31\") != params.end()) this->"<>ToString[paramList[[3]][[i]]]<>"[2][0] = params[\""<>ToString[paramList[[3]][[i]]]<>"31\"];"];
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"32\") != params.end()) this->"<>ToString[paramList[[3]][[i]]]<>"[2][1] = params[\""<>ToString[paramList[[3]][[i]]]<>"32\"];"];
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"33\") != params.end()) this->"<>ToString[paramList[[3]][[i]]]<>"[2][2] = params[\""<>ToString[paramList[[3]][[i]]]<>"33\"];"], 
		{i,1,Length[paramList[[3]]]}]	
	];
	
	WriteLine[line, StringJoin["}"]];
];

BuildUpdater[className_, paramList_, line_]:=Module[{},
	WriteLine[line, ""];
	WriteLine[line, StringJoin["void ", className, "::updateParams", "(std::unordered_map<std::string, double> params) {"]];
	
	If[Length[paramList[[1]]]!=0,
		Do[WriteLine[line, "    if (params.find(\""<>ToString[paramList[[1]][[i]]]<>"\") != params.end()) this->"<>ToString[paramList[[1]][[i]]]<>" = params[\""<>ToString[paramList[[1]][[i]]]<>"\"];"], 
		{i,1,Length[paramList[[1]]]}]	
	];
	
	If[Length[paramList[[2]]]!=0,
		Do[
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[2]][[i]]]<>"1\") != params.end()) this->"<>ToString[paramList[[2]][[i]]]<>"[0] = params[\""<>ToString[paramList[[2]][[i]]]<>"1\"];"];
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[2]][[i]]]<>"2\") != params.end()) this->"<>ToString[paramList[[2]][[i]]]<>"[1] = params[\""<>ToString[paramList[[2]][[i]]]<>"2\"];"];
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[2]][[i]]]<>"3\") != params.end()) this->"<>ToString[paramList[[2]][[i]]]<>"[2] = params[\""<>ToString[paramList[[2]][[i]]]<>"3\"];"], 
		{i,1,Length[paramList[[2]]]}]	
	];
	
	If[Length[paramList[[3]]]!=0,
		Do[
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"11\") != params.end()) this->"<>ToString[paramList[[3]][[i]]]<>"[0][0] = params[\""<>ToString[paramList[[3]][[i]]]<>"11\"];"];
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"12\") != params.end()) this->"<>ToString[paramList[[3]][[i]]]<>"[0][1] = params[\""<>ToString[paramList[[3]][[i]]]<>"12\"];"];
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"13\") != params.end()) this->"<>ToString[paramList[[3]][[i]]]<>"[0][2] = params[\""<>ToString[paramList[[3]][[i]]]<>"13\"];"];
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"21\") != params.end()) this->"<>ToString[paramList[[3]][[i]]]<>"[1][0] = params[\""<>ToString[paramList[[3]][[i]]]<>"21\"];"];
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"22\") != params.end()) this->"<>ToString[paramList[[3]][[i]]]<>"[1][1] = params[\""<>ToString[paramList[[3]][[i]]]<>"22\"];"];
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"23\") != params.end()) this->"<>ToString[paramList[[3]][[i]]]<>"[1][2] = params[\""<>ToString[paramList[[3]][[i]]]<>"23\"];"];
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"31\") != params.end()) this->"<>ToString[paramList[[3]][[i]]]<>"[2][0] = params[\""<>ToString[paramList[[3]][[i]]]<>"31\"];"];
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"32\") != params.end()) this->"<>ToString[paramList[[3]][[i]]]<>"[2][1] = params[\""<>ToString[paramList[[3]][[i]]]<>"32\"];"];
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"33\") != params.end()) this->"<>ToString[paramList[[3]][[i]]]<>"[2][2] = params[\""<>ToString[paramList[[3]][[i]]]<>"33\"];"], 
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


BuildFunctionWarsaw[modelName_,WCname_,expr_,ComplexPars_]:=Module[{returnExpr,fileName,path,line},
	returnExpr = ConvertFullExpression[expr,ComplexPars];
	path = FileNameJoin[{ParentDirectory[NotebookDirectory[]],"lib"}];
	fileName = First[StringSplit[WCname,"("]];
	line = Which[
		MemberQ[WCList[[1]],fileName],
		OpenWrite[path<>"/bosonic/"<>fileName<>".cpp"],
		MemberQ[WCList[[2]],fileName],
		OpenWrite[path<>"/two_fermions/"<>fileName<>".cpp"],
		MemberQ[WCList[[3]],fileName],
		OpenWrite[path<>"/four_fermions/"<>fileName<>".cpp"]
	];
	
	WriteLine[line, "#include \"OperatorImport.h\""];
	WriteLine[line, "#include \""<>modelName<>".h\""];
	WriteLine[line,""];
	WriteLine[line, "double "<>modelName<>"::"<>WCname<> " {"];		
	WriteLine[line, "    return ("<>returnExpr<>");"];
	WriteLine[line, "}"];
	Close[line];
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


SourceFileBuilder[modelName_, paramList_, ComplexPars_, matchingOutput_]:=Module[{keyList, exprList, path,line1},
	keyList=Keys[matchingOutput];
	exprList=Values[matchingOutput];
	
	path = FileNameJoin[{ParentDirectory[NotebookDirectory[]],"lib"}];
	line1 = OpenWrite[path<>"/"<>modelName<>".cpp"];
	
	BuildPreprocessorDirectives[modelName,line1];
	BuildConstructor[modelName,paramList,line1];
	BuildUpdater[modelName,paramList,line1];
	BuildPrinter[modelName,paramList,line1];
	
	Close[line1];
	
	Do[
		BuildFunctionWarsaw[modelName,WarsawAll[keyList[[k]]],exprList[[k]],ComplexPars],
	{k,1,Length[matchingOutput]}];
	
]



(* ::Subsection:: *)
(*Builder for .pyi file with python type information for all WC functions*)


GeneratePythonDeclarations[modelName_]:= Module[{path,line},
	path = FileNameJoin[{ParentDirectory[NotebookDirectory[]],"py"}];
	line = OpenWrite[path<>"/"<>"match_to_py.pyi"];
	WriteLine[line, "class "<>modelName<>":"];
	WriteLine[line, "    def __init__(self, param_dict: dict[str, float]) -> None: ..."];
	WriteLine[line, "    def updateParams(self, param_dict: dict[str, float]) -> None: ..."];
	WriteLine[line, "    def cllHH(self, i1: int, i2: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cG(self, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cW(self, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cGt(self, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cWt(self, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cH(self, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cHBox(self, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cHD(self, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cHG(self, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cHW(self, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cHB(self, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cHWB(self, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cHGt(self, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cHWt(self, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cHBt(self, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cHWtB(self, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def ceH(self, i1: int, i2: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cuH(self, i1: int, i2: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cdH(self, i1: int, i2: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def ceW(self, i1: int, i2: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def ceB(self, i1: int, i2: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cuG(self, i1: int, i2: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cuW(self, i1: int, i2: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cuB(self, i1: int, i2: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cdG(self, i1: int, i2: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cdW(self, i1: int, i2: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cdB(self, i1: int, i2: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cHl1(self, i1: int, i2: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cHl3(self, i1: int, i2: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cHe(self, i1: int, i2: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cHq1(self, i1: int, i2: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cHq3(self, i1: int, i2: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cHu(self, i1: int, i2: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cHd(self, i1: int, i2: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cHud(self, i1: int, i2: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cll(self, i1: int, i2: int, i3: int, i4: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cqq1(self, i1: int, i2: int, i3: int, i4: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cqq3(self, i1: int, i2: int, i3: int, i4: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def clq1(self, i1: int, i2: int, i3: int, i4: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def clq3(self, i1: int, i2: int, i3: int, i4: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cee(self, i1: int, i2: int, i3: int, i4: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cuu(self, i1: int, i2: int, i3: int, i4: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cdd(self, i1: int, i2: int, i3: int, i4: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def ceu(self, i1: int, i2: int, i3: int, i4: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def ced(self, i1: int, i2: int, i3: int, i4: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cud1(self, i1: int, i2: int, i3: int, i4: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cud8(self, i1: int, i2: int, i3: int, i4: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cle(self, i1: int, i2: int, i3: int, i4: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def clu(self, i1: int, i2: int, i3: int, i4: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cld(self, i1: int, i2: int, i3: int, i4: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cqe(self, i1: int, i2: int, i3: int, i4: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cqu1(self, i1: int, i2: int, i3: int, i4: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cqu8(self, i1: int, i2: int, i3: int, i4: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cqd1(self, i1: int, i2: int, i3: int, i4: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cqd8(self, i1: int, i2: int, i3: int, i4: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cledq(self, i1: int, i2: int, i3: int, i4: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cquqd1(self, i1: int, i2: int, i3: int, i4: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cquqd8(self, i1: int, i2: int, i3: int, i4: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def clequ1(self, i1: int, i2: int, i3: int, i4: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def clequ3(self, i1: int, i2: int, i3: int, i4: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cduq(self, i1: int, i2: int, i3: int, i4: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cqqu(self, i1: int, i2: int, i3: int, i4: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cqqq(self, i1: int, i2: int, i3: int, i4: int, mubarsq: float, hbar: float) -> float: ..."];
	WriteLine[line, "    def cduu(self, i1: int, i2: int, i3: int, i4: int, mubarsq: float, hbar: float) -> float: ..."];
	Close[line];
]


(* ::Section:: *)
(*Dictionaries (WCs and Miscellaneous) *)


(* ::Subsection:: *)
(*A dictionary containing all Warsaw basis operator names and the corresponding C++ function prototype*)


(* Here, we are treating all coefficients to be real*)

WarsawAll = Association[
	"cllHH" -> "cllHH(int i1, int i2, double mubarsq, double hbar)",
	
	"cG" -> "cG(double mubarsq, double hbar)",
	"cW" -> "cW(double mubarsq, double hbar)",
	"cGt" -> "cGt(double mubarsq, double hbar)",
	"cWt" -> "cWt(double mubarsq, double hbar)",
	
	"cH" -> "cH(double mubarsq, double hbar)",
	"cHBox" -> "cHBox(double mubarsq, double hbar)",
	"cHD" -> "cHD(double mubarsq, double hbar)",
	
	"cHG" -> "cHG(double mubarsq, double hbar)",
	"cHW" -> "cHW(double mubarsq, double hbar)",
	"cHB" -> "cHB(double mubarsq, double hbar)",
	"cHWB" -> "cHWB(double mubarsq, double hbar)",
	"cHGt" -> "cHGt(double mubarsq, double hbar)",
	"cHWt" -> "cHWt(double mubarsq, double hbar)",
	"cHBt" -> "cHBt(double mubarsq, double hbar)",
	"cHWtB" -> "cHWtB(double mubarsq, double hbar)",
	
	"ceH" -> "ceH(int i1, int i2, double mubarsq, double hbar)",
	"cuH" -> "cuH(int i1, int i2, double mubarsq, double hbar)",
	"cdH" -> "cdH(int i1, int i2, double mubarsq, double hbar)",
	
	"ceW" -> "ceW(int i1, int i2, double mubarsq, double hbar)",
	"ceB" -> "ceB(int i1, int i2, double mubarsq, double hbar)",
	"cuG" -> "cuG(int i1, int i2, double mubarsq, double hbar)",
	"cuW" -> "cuW(int i1, int i2, double mubarsq, double hbar)",
	"cuB" -> "cuB(int i1, int i2, double mubarsq, double hbar)",
	"cdG" -> "cdG(int i1, int i2, double mubarsq, double hbar)",
	"cdW" -> "cdW(int i1, int i2, double mubarsq, double hbar)",
	"cdB" -> "cdB(int i1, int i2, double mubarsq, double hbar)",
	
	"cHl1" -> "cHl1(int i1, int i2, double mubarsq, double hbar)",
	"cHl3" -> "cHl3(int i1, int i2, double mubarsq, double hbar)",
	"cHe" -> "cHe(int i1, int i2, double mubarsq, double hbar)",
	"cHq1" -> "cHq1(int i1, int i2, double mubarsq, double hbar)",
	"cHq3" -> "cHq3(int i1, int i2, double mubarsq, double hbar)",
	"cHu" -> "cHu(int i1, int i2, double mubarsq, double hbar)",
	"cHd" -> "cHd(int i1, int i2, double mubarsq, double hbar)",
	"cHud" -> "cHud(int i1, int i2, double mubarsq, double hbar)",
	
	"cll" -> "cll(int i1, int i2, int i3, int i4, double mubarsq, double hbar)",
	"cqq1" -> "cqq1(int i1, int i2, int i3, int i4, double mubarsq, double hbar)",
	"cqq3" -> "cqq3(int i1, int i2, int i3, int i4, double mubarsq, double hbar)",
	"clq1" -> "clq1(int i1, int i2, int i3, int i4, double mubarsq, double hbar)",
	"clq3" -> "clq3(int i1, int i2, int i3, int i4, double mubarsq, double hbar)",
	"cee" -> "cee(int i1, int i2, int i3, int i4, double mubarsq, double hbar)",
	"cuu" -> "cuu(int i1, int i2, int i3, int i4, double mubarsq, double hbar)",
	"cdd" -> "cdd(int i1, int i2, int i3, int i4, double mubarsq, double hbar)",
	"ceu" -> "ceu(int i1, int i2, int i3, int i4, double mubarsq, double hbar)",
	"ced" -> "ced(int i1, int i2, int i3, int i4, double mubarsq, double hbar)",
	"cud1" -> "cud1(int i1, int i2, int i3, int i4, double mubarsq, double hbar)",
	"cud8" -> "cud8(int i1, int i2, int i3, int i4, double mubarsq, double hbar)",
	"cle" -> "cle(int i1, int i2, int i3, int i4, double mubarsq, double hbar)",
	"clu" -> "clu(int i1, int i2, int i3, int i4, double mubarsq, double hbar)",
	"cld" -> "cld(int i1, int i2, int i3, int i4, double mubarsq, double hbar)",
	"cqe" -> "cqe(int i1, int i2, int i3, int i4, double mubarsq, double hbar)",
	"cqu1" -> "cqu1(int i1, int i2, int i3, int i4, double mubarsq, double hbar)",
	"cqu8" -> "cqu8(int i1, int i2, int i3, int i4, double mubarsq, double hbar)",
	"cqd1" -> "cqd1(int i1, int i2, int i3, int i4, double mubarsq, double hbar)",
	"cqd8" -> "cqd8(int i1, int i2, int i3, int i4, double mubarsq, double hbar)",
	
	"cledq" -> "cledq(int i1, int i2, int i3, int i4, double mubarsq, double hbar)",
	"cquqd1" -> "cquqd1(int i1, int i2, int i3, int i4, double mubarsq, double hbar)",
	"cquqd8" -> "cquqd8(int i1, int i2, int i3, int i4, double mubarsq, double hbar)",
	"clequ1" -> "clequ1(int i1, int i2, int i3, int i4, double mubarsq, double hbar)",
	"clequ3" -> "clequ3(int i1, int i2, int i3, int i4, double mubarsq, double hbar)" ,
	
	"cduq" -> "cduq(int i1, int i2, int i3, int i4, double mubarsq, double hbar)",
	"cqqu" -> "cqqu(int i1, int i2, int i3, int i4, double mubarsq, double hbar)",
	"cqqq" -> "cqqq(int i1, int i2, int i3, int i4, double mubarsq, double hbar)",
	"cduu" -> "cduu(int i1, int i2, int i3, int i4, double mubarsq, double hbar)"
];


(* ::Subsection:: *)
(*List of Wilson Coefficients segregated based on the number of fermions*)


WCList = { 
	{"cG","cW","cGt","cWt","cH","cHBox","cHD","cHG","cHW","cHB","cHWB","cHGt","cHWt","cHBt","cHWtB"},
	{"cllHH","ceH","cuH","cdH","ceW","ceB","cuG","cuW","cuB","cdG","cdW","cdB","cHl1","cHl3","cHe","cHq1","cHq3","cHu","cHd","cHud"},
	{"cll","cqq1","cqq3","clq1","clq3","cee","cuu","cdd","ceu","ced","cud1","cud8",
	 "cle","clu","cld","cqe","cqu1","cqu8","cqd1","cqd8","cledq","cquqd1","cquqd8",
	 "clequ1","clequ3","cduq","cqqu","cqqq","cduu"}
};


(* ::Chapter:: *)
(*-- Package Context --*)


End[];
EndPackage[];
