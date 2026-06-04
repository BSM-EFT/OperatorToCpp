(* ::Package:: *)

(* ::Chapter:: *)
(*-- Package Context --*)


BeginPackage["OperatorExport`", {"Global`"}];


(* ::Subsection:: *)
(*Description of public functions*)


SimplifyOutput::usage = "SimplifyOutput[matchedResult] creates a dictionary containing only Wilson coefficients and the corresponding matching conditions, while unwraping the Matchete API and converting the expressions to a lighter intermediate form.";
SegregateParams::usage = "SegregateParams[SimplifiedOutput, ComplexParams] collects all unique parameters within the simplified matching output and segregates them based on their dimensionality. Keeps track of the complex parameters to ensure that their conjugates are not counted as distinct parameters.";
HeaderFileBuilder::usage = "HeaderFileBuilder[ModelName,ModelParams] creates a ModelName.h file that defines a ModelName class with the model parameters as member variables and Warsaw basis Wilson coefficients as methods. Also, declares methods for initalizing, updating and printing the parameters of a ModelName object.";
SourceFileBuilder::usage = "SourceFileBuilder[ModelName,ModelParams,ComplexParams,SimplifiedOutput] creates a ModelName.cpp file with implementations for each Wilson coefficient method corresponding to the ModelName class. The method bodies are built by further manipulating the SimplifiedOutput.";
GeneratePythonDeclarations::usage = "GeneratePythonDeclarations[ModelName] creates a ModelName.pyi file with declarations for the model class and each Wilson coefficient method.";
BuildFiles::usage = "BuildFiles[ModelName,SimplifiedOutput,ComplexParams] offers a single command to build C++ header (.h) and source (.cpp) files, along with Python declarations of the Model class and its methods.";


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
      "OperatorExport  v0.2 \[LongDash] by ",
      "Suraj Prakash",
      " (", 
      Style["suraj.prakash@ific.uv.es", "Hyperlink", FontColor -> Blue],
      ")"
    }],
    Row[{
      "Affiliation: IFIC (Universitat de Valencia - CSIC)"
    }],
    Row[{
      "Repository: ", 
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
	dict = dict /. {Index[x_,__]:>ToExpression[StringJoin["r",StringPart[ToString[x],4;;-1]]]};
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
	Table[ToExpression["r"<>ToString[j]]->j,{j,1,12}], (* this upper limit of 12 has been added manually *)
	Table[ToExpression["i"<>ToString[j]]->100+j,{j,1,12}]
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
	(*newStr=StringReplace[newStr,ReplaceCCString[ComplexPars]];*)
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
	WriteLine[line, "#include <complex>"];
	WriteLine[line, "#include <unordered_map>"];
	WriteLine[line, "#include <map>"];
	WriteLine[line, "#include <functional>"];
	WriteLine[line, "#include <cmath>"];
	WriteLine[line, ""];
	WriteLine[line, "const double pi = 3.14159265;"];
];

HeaderTaskStruct[line_] := Module[{},
	WriteLine[line, "struct Task {"];
	WriteLine[line, "    std::string name;"];
	WriteLine[line, "    std::function<std::complex<double>()> work;"];
	WriteLine[line, "};"];
]

HeaderModelClass[className_,paramList_,ComplexPars_,line_] := Module[{},
	WriteLine[line, "class "<>className<>" {"];
	WriteLine[line, "    private:"];
	
	(* define and initialize the scale and loop-order *)
	WriteLine[line, "        double hbar = 1/(16 * pow(pi,2));"];
	WriteLine[line, "        double mubarsq = 1.0;"];
	(* define and initialize 0-dimensional parameters *)
	If[Length[paramList[[1]]]!=0,
		Do[
			WriteLine[line, "        "<>"std::complex<double> "<>ToString[paramList[[1]][[i]]]<>" = 0.0;"];
			If[
				MemberQ[ComplexPars,paramList[[1]][[i]]],
				WriteLine[line, "        "<>"std::complex<double> "<>ToString[paramList[[1]][[i]]]<>"c = 0.0;"];
			], {i,1,Length[paramList[[1]]]}
		];
	];
	WriteLine[line, ""];
	
	(* define and initialize 1-dimensional parameters (masses) *)
	If[Length[paramList[[2]]]!=0,
		Do[
			WriteLine[line, "        "<>"std::vector<std::complex<double> > "<>ToString[paramList[[2]][[i]]]<>" = {0.0, 0.0, 0.0};"];
			If[
				MemberQ[ComplexPars,paramList[[2]][[i]]],
				WriteLine[line, "        "<>"std::vector<std::complex<double> > "<>ToString[paramList[[2]][[i]]]<>"c = {0.0, 0.0, 0.0};"];
			], {i,1,Length[paramList[[2]]]}
		];
	];
	WriteLine[line, ""];
	
	(* define and initialize 2-dimensional parameters *)
	If[Length[paramList[[3]]]!=0,
		Do[
			WriteLine[line, "        "<>"std::vector<std::vector<std::complex<double> > > "<>ToString[paramList[[3]][[i]]]<>" = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};"];
			If[
				MemberQ[ComplexPars,paramList[[3]][[i]]],
				WriteLine[line, "        "<>"std::vector<std::vector<std::complex<double> > > "<>ToString[paramList[[3]][[i]]]<>"c = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};"];
			],{i,1,Length[paramList[[3]]]}
		];
	];
	WriteLine[line, ""];
		
	WriteLine[line, "    public:"];
	(* declaration for the default constructor*)
	WriteLine[line, "        "<>className<>"() = default;"];
	WriteLine[line, ""];
	
	(* declaration for the overloaded constructor*)
	WriteLine[line, StringJoin["        ", className, "(std::unordered_map<std::string, std::complex<double> > params, double scale, bool loop);"]];
	WriteLine[line, ""];

	(* declaration for the updater method *)
	WriteLine[line, StringJoin["        void updateParams", "(std::unordered_map<std::string, std::complex<double> > params);"]];
	WriteLine[line, ""];
	
	(* declaration for the scale getter method *)
	WriteLine[line, StringJoin["        double getScale", "();"]];
	WriteLine[line, ""];
	
	(* declaration for the scale updater method *)
	WriteLine[line, StringJoin["        void setScale", "(double scale);"]];
	WriteLine[line, ""];
	
	(* declaration for the loop-order updater method *)
	WriteLine[line, StringJoin["        void loopContributions", "(bool loop);"]];
	WriteLine[line, ""];
	
	(* declaration for a getter method for model parameters *)
	WriteLine[line, StringJoin["        std::unordered_map<std::string, std::complex<double> > getParams", "();"]];
	WriteLine[line, ""];
		
	(* declarations for WC methods (Warsaw) *)
	Do[WriteLine[line,"        std::complex<double> "<>Values[WarsawAll][[i]]<>";"],{i,1,Length[WarsawAll]}];
	
	WriteLine[line, ""];
	
	(* declaration for the batch evaluator *)
	WriteLine[line, "        std::map<std::string, std::complex<double> > batch_eval(const std::vector<Task>& tasks);"];
	WriteLine[line, "};"];
	WriteLine[line, ""];
];


(* ::Subsubsection:: *)
(*Master builder*)


HeaderFileBuilder[modelName_,paramList_,ComplexPars_]:=Module[{path,line1,line2},
	path = FileNameJoin[{ParentDirectory[NotebookDirectory[]],"include"}];
	
	line1 = OpenWrite[path<>"/"<>modelName<>".h"];
	HeaderPreprocessorDirectives[line1];
	WriteLine[line1,""];
	HeaderTaskStruct[line1];
	WriteLine[line1,""];
	HeaderModelClass[modelName,paramList,ComplexPars,line1];
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
	WriteLine[line, "#include \""<>modelName<>".h\""];
	WriteLine[line, "#include <omp.h>"];
];

BuildConstructor[className_, paramList_, ComplexPars_, line_]:=Module[{},
	WriteLine[line, ""];
	WriteLine[line, StringJoin[className, "::", className, "(std::unordered_map<std::string, std::complex<double>> params, double scale, bool loop) {"]];
	WriteLine[line, "    this->mubarsq = scale * scale;"];
	WriteLine[line, "    if (!loop) {"];
	WriteLine[line, "        this->hbar = 0.0;"];
	WriteLine[line, "    } else {"];
	WriteLine[line, "        this->hbar = 1/(16 * pow(pi,2));"];
	WriteLine[line, "    }"];
	WriteLine[line, ""];
	
	If[Length[paramList[[1]]]!=0,
		Do[
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[1]][[i]]]<>"\") != params.end()) {"];
			WriteLine[line, "        this->"<>ToString[paramList[[1]][[i]]]<>" = params[\""<>ToString[paramList[[1]][[i]]]<>"\"];"];
			If[
				MemberQ[ComplexPars,paramList[[1]][[i]]],
				WriteLine[line, "        this->"<>ToString[paramList[[1]][[i]]]<>"c = std::conj(params[\""<>ToString[paramList[[1]][[i]]]<>"\"]);"];];
			WriteLine[line, "    }"], {i,1,Length[paramList[[1]]]}]	
	];
	
	If[Length[paramList[[2]]]!=0,
		Do[
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[2]][[i]]]<>"1\") != params.end()) {"];
			WriteLine[line, "        this->"<>ToString[paramList[[2]][[i]]]<>"[0] = params[\""<>ToString[paramList[[2]][[i]]]<>"1\"];"];
			If[
				MemberQ[ComplexPars,paramList[[2]][[i]]],
				WriteLine[line, "        this->"<>ToString[paramList[[2]][[i]]]<>"c[0] = std::conj(params[\""<>ToString[paramList[[2]][[i]]]<>"1\"]);"];
			];
			WriteLine[line, "    }"];
			
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[2]][[i]]]<>"2\") != params.end()) {"];
			WriteLine[line, "        this->"<>ToString[paramList[[2]][[i]]]<>"[1] = params[\""<>ToString[paramList[[2]][[i]]]<>"2\"];"];
			If[
				MemberQ[ComplexPars,paramList[[2]][[i]]],
				WriteLine[line, "        this->"<>ToString[paramList[[2]][[i]]]<>"c[1] = std::conj(params[\""<>ToString[paramList[[2]][[i]]]<>"2\"]);"];
			];
			WriteLine[line, "    }"];
			
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[2]][[i]]]<>"3\") != params.end()) {"];
			WriteLine[line, "        this->"<>ToString[paramList[[2]][[i]]]<>"[2] = params[\""<>ToString[paramList[[2]][[i]]]<>"3\"];"];
			If[
				MemberQ[ComplexPars,paramList[[2]][[i]]],
				WriteLine[line, "        this->"<>ToString[paramList[[2]][[i]]]<>"c[2] = std::conj(params[\""<>ToString[paramList[[2]][[i]]]<>"3\"]);"];
			];
			WriteLine[line, "    }"], 
			{i,1,Length[paramList[[2]]]}]	
	];
	
	If[Length[paramList[[3]]]!=0,
		Do[
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"11\") != params.end()) {"];
			WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"[0][0] = params[\""<>ToString[paramList[[3]][[i]]]<>"11\"];"];
			If[
				MemberQ[ComplexPars,paramList[[3]][[i]]],
				WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"c[0][0] = std::conj(params[\""<>ToString[paramList[[3]][[i]]]<>"11\"]);"];
			];
			WriteLine[line, "    }"];
			
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"12\") != params.end()) {"];
			WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"[0][1] = params[\""<>ToString[paramList[[3]][[i]]]<>"12\"];"];
			If[
				MemberQ[ComplexPars,paramList[[3]][[i]]],
				WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"c[1][0] = std::conj(params[\""<>ToString[paramList[[3]][[i]]]<>"12\"]);"];
			];
			WriteLine[line, "    }"];
			
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"13\") != params.end()) {"];
			WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"[0][2] = params[\""<>ToString[paramList[[3]][[i]]]<>"13\"];"];
			If[
				MemberQ[ComplexPars,paramList[[3]][[i]]],
				WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"c[2][0] = std::conj(params[\""<>ToString[paramList[[3]][[i]]]<>"13\"]);"];
			];
			WriteLine[line, "    }"];
			
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"21\") != params.end()) {"];
			WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"[1][0] = params[\""<>ToString[paramList[[3]][[i]]]<>"21\"];"];
			If[
				MemberQ[ComplexPars,paramList[[3]][[i]]],
				WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"c[0][1] = std::conj(params[\""<>ToString[paramList[[3]][[i]]]<>"21\"]);"];
			];
			WriteLine[line, "    }"];
			
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"22\") != params.end()) {"];
			WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"[1][1] = params[\""<>ToString[paramList[[3]][[i]]]<>"22\"];"];
			If[
				MemberQ[ComplexPars,paramList[[3]][[i]]],
				WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"c[1][1] = std::conj(params[\""<>ToString[paramList[[3]][[i]]]<>"22\"]);"];
			];
			WriteLine[line, "    }"];
			
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"23\") != params.end()) {"];
			WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"[1][2] = params[\""<>ToString[paramList[[3]][[i]]]<>"23\"];"];
			If[
				MemberQ[ComplexPars,paramList[[3]][[i]]],
				WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"c[2][1] = std::conj(params[\""<>ToString[paramList[[3]][[i]]]<>"23\"]);"];
			];
			WriteLine[line, "    }"];
			
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"31\") != params.end()) {"];
			WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"[2][0] = params[\""<>ToString[paramList[[3]][[i]]]<>"31\"];"];
			If[
				MemberQ[ComplexPars,paramList[[3]][[i]]],
				WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"c[0][2] = std::conj(params[\""<>ToString[paramList[[3]][[i]]]<>"31\"]);"];
			];
			WriteLine[line, "    }"];
			
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"32\") != params.end()) {"];
			WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"[2][1] = params[\""<>ToString[paramList[[3]][[i]]]<>"32\"];"];
			If[
				MemberQ[ComplexPars,paramList[[3]][[i]]],
				WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"c[1][2] = std::conj(params[\""<>ToString[paramList[[3]][[i]]]<>"32\"]);"];
			];
			WriteLine[line, "    }"];
			
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"33\") != params.end()) {"];
			WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"[2][2] = params[\""<>ToString[paramList[[3]][[i]]]<>"33\"];"];
			If[
				MemberQ[ComplexPars,paramList[[3]][[i]]],
				WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"c[2][2] = std::conj(params[\""<>ToString[paramList[[3]][[i]]]<>"33\"]);"];
			];
			WriteLine[line, "    }"];
			
			
		,{i,1,Length[paramList[[3]]]}]	
	];
	
	WriteLine[line, StringJoin["}"]];
];

BuildUpdater[className_, paramList_, ComplexPars_, line_]:=Module[{},
	WriteLine[line, ""];
	WriteLine[line, StringJoin["void ", className, "::updateParams", "(std::unordered_map<std::string, std::complex<double>> params) {"]];
	
	If[Length[paramList[[1]]]!=0,
		Do[
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[1]][[i]]]<>"\") != params.end()) {"];
			WriteLine[line, "        this->"<>ToString[paramList[[1]][[i]]]<>" = params[\""<>ToString[paramList[[1]][[i]]]<>"\"];"];
			If[
				MemberQ[ComplexPars,paramList[[1]][[i]]],
				WriteLine[line, "        this->"<>ToString[paramList[[1]][[i]]]<>"c = std::conj(params[\""<>ToString[paramList[[1]][[i]]]<>"\"]);"];];
			WriteLine[line, "    }"], {i,1,Length[paramList[[1]]]}]	
	];
	
	If[Length[paramList[[2]]]!=0,
		Do[
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[2]][[i]]]<>"1\") != params.end()) {"];
			WriteLine[line, "        this->"<>ToString[paramList[[2]][[i]]]<>"[0] = params[\""<>ToString[paramList[[2]][[i]]]<>"1\"];"];
			If[
				MemberQ[ComplexPars,paramList[[2]][[i]]],
				WriteLine[line, "        this->"<>ToString[paramList[[2]][[i]]]<>"c[0] = std::conj(params[\""<>ToString[paramList[[2]][[i]]]<>"1\"]);"];
			];
			WriteLine[line, "    }"];
			
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[2]][[i]]]<>"2\") != params.end()) {"];
			WriteLine[line, "        this->"<>ToString[paramList[[2]][[i]]]<>"[1] = params[\""<>ToString[paramList[[2]][[i]]]<>"2\"];"];
			If[
				MemberQ[ComplexPars,paramList[[2]][[i]]],
				WriteLine[line, "        this->"<>ToString[paramList[[2]][[i]]]<>"c[1] = std::conj(params[\""<>ToString[paramList[[2]][[i]]]<>"2\"]);"];
			];
			WriteLine[line, "    }"];
			
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[2]][[i]]]<>"3\") != params.end()) {"];
			WriteLine[line, "        this->"<>ToString[paramList[[2]][[i]]]<>"[2] = params[\""<>ToString[paramList[[2]][[i]]]<>"3\"];"];
			If[
				MemberQ[ComplexPars,paramList[[2]][[i]]],
				WriteLine[line, "        this->"<>ToString[paramList[[2]][[i]]]<>"c[2] = std::conj(params[\""<>ToString[paramList[[2]][[i]]]<>"3\"]);"];
			];
			WriteLine[line, "    }"], 
			{i,1,Length[paramList[[2]]]}]	
	];
	
	If[Length[paramList[[3]]]!=0,
		Do[
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"11\") != params.end()) {"];
			WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"[0][0] = params[\""<>ToString[paramList[[3]][[i]]]<>"11\"];"];
			If[
				MemberQ[ComplexPars,paramList[[3]][[i]]],
				WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"c[0][0] = std::conj(params[\""<>ToString[paramList[[3]][[i]]]<>"11\"]);"];
			];
			WriteLine[line, "    }"];
			
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"12\") != params.end()) {"];
			WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"[0][1] = params[\""<>ToString[paramList[[3]][[i]]]<>"12\"];"];
			If[
				MemberQ[ComplexPars,paramList[[3]][[i]]],
				WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"c[1][0] = std::conj(params[\""<>ToString[paramList[[3]][[i]]]<>"12\"]);"];
			];
			WriteLine[line, "    }"];
			
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"13\") != params.end()) {"];
			WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"[0][2] = params[\""<>ToString[paramList[[3]][[i]]]<>"13\"];"];
			If[
				MemberQ[ComplexPars,paramList[[3]][[i]]],
				WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"c[2][0] = std::conj(params[\""<>ToString[paramList[[3]][[i]]]<>"13\"]);"];
			];
			WriteLine[line, "    }"];
			
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"21\") != params.end()) {"];
			WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"[1][0] = params[\""<>ToString[paramList[[3]][[i]]]<>"21\"];"];
			If[
				MemberQ[ComplexPars,paramList[[3]][[i]]],
				WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"c[0][1] = std::conj(params[\""<>ToString[paramList[[3]][[i]]]<>"21\"]);"];
			];
			WriteLine[line, "    }"];
			
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"22\") != params.end()) {"];
			WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"[1][1] = params[\""<>ToString[paramList[[3]][[i]]]<>"22\"];"];
			If[
				MemberQ[ComplexPars,paramList[[3]][[i]]],
				WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"c[1][1] = std::conj(params[\""<>ToString[paramList[[3]][[i]]]<>"22\"]);"];
			];
			WriteLine[line, "    }"];
			
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"23\") != params.end()) {"];
			WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"[1][2] = params[\""<>ToString[paramList[[3]][[i]]]<>"23\"];"];
			If[
				MemberQ[ComplexPars,paramList[[3]][[i]]],
				WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"c[2][1] = std::conj(params[\""<>ToString[paramList[[3]][[i]]]<>"23\"]);"];
			];
			WriteLine[line, "    }"];
			
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"31\") != params.end()) {"];
			WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"[2][0] = params[\""<>ToString[paramList[[3]][[i]]]<>"31\"];"];
			If[
				MemberQ[ComplexPars,paramList[[3]][[i]]],
				WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"c[0][2] = std::conj(params[\""<>ToString[paramList[[3]][[i]]]<>"31\"]);"];
			];
			WriteLine[line, "    }"];
			
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"32\") != params.end()) {"];
			WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"[2][1] = params[\""<>ToString[paramList[[3]][[i]]]<>"32\"];"];
			If[
				MemberQ[ComplexPars,paramList[[3]][[i]]],
				WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"c[1][2] = std::conj(params[\""<>ToString[paramList[[3]][[i]]]<>"32\"]);"];
			];
			WriteLine[line, "    }"];
			
			WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>"33\") != params.end()) {"];
			WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"[2][2] = params[\""<>ToString[paramList[[3]][[i]]]<>"33\"];"];
			If[
				MemberQ[ComplexPars,paramList[[3]][[i]]],
				WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"c[2][2] = std::conj(params[\""<>ToString[paramList[[3]][[i]]]<>"33\"]);"];
			];
			WriteLine[line, "    }"];
			
			
		,{i,1,Length[paramList[[3]]]}]	
	];
	
	WriteLine[line, "}"];
];

BuildScaleGetter[className_,line_]:=Module[{},
	WriteLine[line, ""];
	WriteLine[line, "double "<>className<>"::getScale() {"];
	WriteLine[line, "    return sqrt(this->mubarsq);"];
	WriteLine[line, "}"];
]

BuildScaleSetter[className_,line_]:=Module[{},
	WriteLine[line, ""];
	WriteLine[line, "void "<>className<>"::setScale(double scale) {"];
	WriteLine[line, "    this->mubarsq = scale * scale;"];
	WriteLine[line, "}"];
]

BuildLoopSwitch[className_,line_]:=Module[{},
	WriteLine[line, ""];
	WriteLine[line, "void "<>className<>"::loopContributions(bool loop) {"];
	WriteLine[line, "    if (!loop) {"];
	WriteLine[line, "        this->hbar = 0.0;"];
	WriteLine[line, "    } else {"];
	WriteLine[line, "        this->hbar = 1/(16 * pow(pi,2));"];
	WriteLine[line, "    }"];
	WriteLine[line, "}"];
]

BuildParamDictGetter[className_, paramList_, line_]:=Module[{},
	WriteLine[line, ""];
	WriteLine[line, "std::unordered_map<std::string, std::complex<double> > "<>className<>"::getParams() {"];
	WriteLine[line, "    std::unordered_map<std::string, std::complex<double> > param_dict = {"];
	If[Length[paramList[[1]]]!=0,
		Do[WriteLine[line, "        {\""<>ToString[paramList[[1]][[i]]]<>"\", this->"<>ToString[paramList[[1]][[i]]]<>"},"], {i,1,Length[paramList[[1]]]}]	
	];
	If[Length[paramList[[2]]]!=0,
		Do[
			WriteLine[line, "        {\""<>ToString[paramList[[2]][[i]]]<>"1\", this->"<>ToString[paramList[[2]][[i]]]<>"[0]},"];
			WriteLine[line, "        {\""<>ToString[paramList[[2]][[i]]]<>"2\", this->"<>ToString[paramList[[2]][[i]]]<>"[1]},"];
			WriteLine[line, "        {\""<>ToString[paramList[[2]][[i]]]<>"3\", this->"<>ToString[paramList[[2]][[i]]]<>"[2]},"], {i,1,Length[paramList[[2]]]}]	
	];
	If[Length[paramList[[3]]]!=0,
		Do[
			WriteLine[line, "        {\""<>ToString[paramList[[3]][[i]]]<>"11\", this->"<>ToString[paramList[[3]][[i]]]<>"[0][0]},"];
			WriteLine[line, "        {\""<>ToString[paramList[[3]][[i]]]<>"12\", this->"<>ToString[paramList[[3]][[i]]]<>"[0][1]},"];
			WriteLine[line, "        {\""<>ToString[paramList[[3]][[i]]]<>"13\", this->"<>ToString[paramList[[3]][[i]]]<>"[0][2]},"];
			WriteLine[line, "        {\""<>ToString[paramList[[3]][[i]]]<>"22\", this->"<>ToString[paramList[[3]][[i]]]<>"[1][1]},"];
			WriteLine[line, "        {\""<>ToString[paramList[[3]][[i]]]<>"23\", this->"<>ToString[paramList[[3]][[i]]]<>"[1][2]},"];
			WriteLine[line, "        {\""<>ToString[paramList[[3]][[i]]]<>"33\", this->"<>ToString[paramList[[3]][[i]]]<>"[2][2]},"], {i,1,Length[paramList[[3]]]}]	
	];
	WriteLine[line, "    };"];
	WriteLine[line, "    return param_dict;"];
	WriteLine[line, "}"];

]

BuildBatchEvaluator[className_, line_]:=Module[{},
	WriteLine[line, ""];
	WriteLine[line, "std::map<std::string, std::complex<double> > "<>className<>"::batch_eval(const std::vector<Task>& tasks) {"];
	WriteLine[line, "    int n = tasks.size();"];
	WriteLine[line, "    std::vector<std::complex<double> > results_temp(n);"];
	WriteLine[line, ""];
	WriteLine[line, "    #pragma omp parallel for schedule(dynamic)"];
	WriteLine[line, "    for (int i = 0; i < n; ++i) {"];
	WriteLine[line, "        results_temp[i] = tasks[i].work();"];
	WriteLine[line, "    }"];
	WriteLine[line, ""];
	WriteLine[line, "    std::map<std::string, std::complex<double> > results;"];
	WriteLine[line, "    for (int i = 0; i < n; ++i) {"];
	WriteLine[line, "        results[tasks[i].name] = results_temp[i];"];
	WriteLine[line, "    }"];
	WriteLine[line, ""];
	WriteLine[line, "    return results;"];
	WriteLine[line, "}"];
]


(* ::Subsubsection:: *)
(*Builder for a single WC function (Warsaw Basis)*)


BuildFunctionWarsaw[modelName_,WCname_,expr_,ComplexPars_]:=Module[{returnExpr,fileName,path,line},
	returnExpr = ConvertFullExpression[expr,ComplexPars];
	path = FileNameJoin[{ParentDirectory[NotebookDirectory[]],"lib"}];
	fileName = First[StringSplit[WCname,"("]];
	line = Which[
		MemberQ[WCList[[1]],fileName],
		OpenWrite[path<>"/0f/"<>fileName<>".cpp"],
		MemberQ[WCList[[2]],fileName],
		OpenWrite[path<>"/2f/"<>fileName<>".cpp"],
		MemberQ[WCList[[3]],fileName],
		OpenWrite[path<>"/4f/"<>fileName<>".cpp"]
	];
	
	WriteLine[line, "#include \"OperatorImport.h\""];
	WriteLine[line, "#include \"complex_math.h\""];
	WriteLine[line, "#include \""<>modelName<>".h\""];
	WriteLine[line,""];
	WriteLine[line, "std::complex<double> "<>modelName<>"::"<>WCname<> " {"];		
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
	BuildConstructor[modelName,paramList,ComplexPars,line1];
	BuildUpdater[modelName,paramList,ComplexPars,line1];
	BuildScaleGetter[modelName,line1];
	BuildScaleSetter[modelName,line1];
	BuildLoopSwitch[modelName,line1];
	BuildParamDictGetter[modelName,paramList,line1];
	BuildBatchEvaluator[modelName,line1];
	
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
	WriteLine[line, "class Task: ..."];
	WriteLine[line, ""];
	WriteLine[line, "class "<>modelName<>":"];
	WriteLine[line, "    def __init__(self, param_dict: dict[str, complex], scale: float, loop: bool) -> None: ..."];
	WriteLine[line, "    def updateParams(self, param_dict: dict[str, complex]) -> None: ..."];
	WriteLine[line, "    def getScale(self) -> float: ..."];
	WriteLine[line, "    def setScale(self, scale: float) -> None: ..."];
	WriteLine[line, "    def loopContributions(self, loop: bool) -> None: ..."];
	WriteLine[line, "    def getParams(self) -> dict[str, complex]: ..."];
	WriteLine[line, "    def batch_eval(self, tasks: list[Task]) -> dict[str, complex]: ..."];
	WriteLine[line, "    def wrap_0f(self, name: str, method_name: str) -> Task: ..."];
	WriteLine[line, "    def wrap_2f(self, name: str, method_name: str, i1: int, i2: int) -> Task: ..."];
	WriteLine[line, "    def wrap_4f(self, name: str, method_name: str, i1: int, i2: int, i3: int, i4: int) -> Task: ..."];
	WriteLine[line, "    def cllHH(self, i1: int, i2: int) -> complex: ..."];
	WriteLine[line, "    def cG(self) -> complex: ..."];
	WriteLine[line, "    def cW(self) -> complex: ..."];
	WriteLine[line, "    def cGt(self) -> complex: ..."];
	WriteLine[line, "    def cWt(self) -> complex: ..."];
	WriteLine[line, "    def cH(self) -> complex: ..."];
	WriteLine[line, "    def cHBox(self) -> complex: ..."];
	WriteLine[line, "    def cHD(self) -> complex: ..."];
	WriteLine[line, "    def cHG(self) -> complex: ..."];
	WriteLine[line, "    def cHW(self) -> complex: ..."];
	WriteLine[line, "    def cHB(self) -> complex: ..."];
	WriteLine[line, "    def cHWB(self) -> complex: ..."];
	WriteLine[line, "    def cHGt(self) -> complex: ..."];
	WriteLine[line, "    def cHWt(self) -> complex: ..."];
	WriteLine[line, "    def cHBt(self) -> complex: ..."];
	WriteLine[line, "    def cHWtB(self) -> complex: ..."];
	WriteLine[line, "    def ceH(self, i1: int, i2: int) -> complex: ..."];
	WriteLine[line, "    def cuH(self, i1: int, i2: int) -> complex: ..."];
	WriteLine[line, "    def cdH(self, i1: int, i2: int) -> complex: ..."];
	WriteLine[line, "    def ceW(self, i1: int, i2: int) -> complex: ..."];
	WriteLine[line, "    def ceB(self, i1: int, i2: int) -> complex: ..."];
	WriteLine[line, "    def cuG(self, i1: int, i2: int) -> complex: ..."];
	WriteLine[line, "    def cuW(self, i1: int, i2: int) -> complex: ..."];
	WriteLine[line, "    def cuB(self, i1: int, i2: int) -> complex: ..."];
	WriteLine[line, "    def cdG(self, i1: int, i2: int) -> complex: ..."];
	WriteLine[line, "    def cdW(self, i1: int, i2: int) -> complex: ..."];
	WriteLine[line, "    def cdB(self, i1: int, i2: int) -> complex: ..."];
	WriteLine[line, "    def cHl1(self, i1: int, i2: int) -> complex: ..."];
	WriteLine[line, "    def cHl3(self, i1: int, i2: int) -> complex: ..."];
	WriteLine[line, "    def cHe(self, i1: int, i2: int) -> complex: ..."];
	WriteLine[line, "    def cHq1(self, i1: int, i2: int) -> complex: ..."];
	WriteLine[line, "    def cHq3(self, i1: int, i2: int) -> complex: ..."];
	WriteLine[line, "    def cHu(self, i1: int, i2: int) -> complex: ..."];
	WriteLine[line, "    def cHd(self, i1: int, i2: int) -> complex: ..."];
	WriteLine[line, "    def cHud(self, i1: int, i2: int) -> complex: ..."];
	WriteLine[line, "    def cll(self, i1: int, i2: int, i3: int, i4: int) -> complex: ..."];
	WriteLine[line, "    def cqq1(self, i1: int, i2: int, i3: int, i4: int) -> complex: ..."];
	WriteLine[line, "    def cqq3(self, i1: int, i2: int, i3: int, i4: int) -> complex: ..."];
	WriteLine[line, "    def clq1(self, i1: int, i2: int, i3: int, i4: int) -> complex: ..."];
	WriteLine[line, "    def clq3(self, i1: int, i2: int, i3: int, i4: int) -> complex: ..."];
	WriteLine[line, "    def cee(self, i1: int, i2: int, i3: int, i4: int) -> complex: ..."];
	WriteLine[line, "    def cuu(self, i1: int, i2: int, i3: int, i4: int) -> complex: ..."];
	WriteLine[line, "    def cdd(self, i1: int, i2: int, i3: int, i4: int) -> complex: ..."];
	WriteLine[line, "    def ceu(self, i1: int, i2: int, i3: int, i4: int) -> complex: ..."];
	WriteLine[line, "    def ced(self, i1: int, i2: int, i3: int, i4: int) -> complex: ..."];
	WriteLine[line, "    def cud1(self, i1: int, i2: int, i3: int, i4: int) -> complex: ..."];
	WriteLine[line, "    def cud8(self, i1: int, i2: int, i3: int, i4: int) -> complex: ..."];
	WriteLine[line, "    def cle(self, i1: int, i2: int, i3: int, i4: int) -> complex: ..."];
	WriteLine[line, "    def clu(self, i1: int, i2: int, i3: int, i4: int) -> complex: ..."];
	WriteLine[line, "    def cld(self, i1: int, i2: int, i3: int, i4: int) -> complex: ..."];
	WriteLine[line, "    def cqe(self, i1: int, i2: int, i3: int, i4: int) -> complex: ..."];
	WriteLine[line, "    def cqu1(self, i1: int, i2: int, i3: int, i4: int) -> complex: ..."];
	WriteLine[line, "    def cqu8(self, i1: int, i2: int, i3: int, i4: int) -> complex: ..."];
	WriteLine[line, "    def cqd1(self, i1: int, i2: int, i3: int, i4: int) -> complex: ..."];
	WriteLine[line, "    def cqd8(self, i1: int, i2: int, i3: int, i4: int) -> complex: ..."];
	WriteLine[line, "    def cledq(self, i1: int, i2: int, i3: int, i4: int) -> complex: ..."];
	WriteLine[line, "    def cquqd1(self, i1: int, i2: int, i3: int, i4: int) -> complex: ..."];
	WriteLine[line, "    def cquqd8(self, i1: int, i2: int, i3: int, i4: int) -> complex: ..."];
	WriteLine[line, "    def clequ1(self, i1: int, i2: int, i3: int, i4: int) -> complex: ..."];
	WriteLine[line, "    def clequ3(self, i1: int, i2: int, i3: int, i4: int) -> complex: ..."];
	WriteLine[line, "    def cduq(self, i1: int, i2: int, i3: int, i4: int) -> complex: ..."];
	WriteLine[line, "    def cqqu(self, i1: int, i2: int, i3: int, i4: int) -> complex: ..."];
	WriteLine[line, "    def cqqq(self, i1: int, i2: int, i3: int, i4: int) -> complex: ..."];
	WriteLine[line, "    def cduu(self, i1: int, i2: int, i3: int, i4: int) -> complex: ..."];
	Close[line];
]


(* ::Subsection:: *)
(*Single function to call all file builders*)


BuildFiles[modelName_,simplifiedOutput_,complexPars_]:=Module[{modelPars},
	modelPars=SegregateParams[simplifiedOutput,complexPars];
	HeaderFileBuilder[modelName,modelPars,complexPars];
	Print["C++ header (.h) files successfully created."];
	SourceFileBuilder[modelName,modelPars,complexPars,simplifiedOutput];
	Print["C++ source (.cpp) files successfully created."];
	GeneratePythonDeclarations[modelName];
	Print["Python class and method declarations generated."]
]


(* ::Section:: *)
(*Dictionaries (WCs and Miscellaneous) *)


(* ::Subsection:: *)
(*A dictionary containing all Warsaw basis operator names and the corresponding C++ function prototype*)


(* Here, we are treating all coefficients to be real*)

WarsawAll = Association[
	"cllHH" -> "cllHH(int i1, int i2)",
	
	"cG" -> "cG()",
	"cW" -> "cW()",
	"cGt" -> "cGt()",
	"cWt" -> "cWt()",
	
	"cH" -> "cH()",
	"cHBox" -> "cHBox()",
	"cHD" -> "cHD()",
	
	"cHG" -> "cHG()",
	"cHW" -> "cHW()",
	"cHB" -> "cHB()",
	"cHWB" -> "cHWB()",
	"cHGt" -> "cHGt()",
	"cHWt" -> "cHWt()",
	"cHBt" -> "cHBt()",
	"cHWtB" -> "cHWtB()",
	
	"ceH" -> "ceH(int i1, int i2)",
	"cuH" -> "cuH(int i1, int i2)",
	"cdH" -> "cdH(int i1, int i2)",
	
	"ceW" -> "ceW(int i1, int i2)",
	"ceB" -> "ceB(int i1, int i2)",
	"cuG" -> "cuG(int i1, int i2)",
	"cuW" -> "cuW(int i1, int i2)",
	"cuB" -> "cuB(int i1, int i2)",
	"cdG" -> "cdG(int i1, int i2)",
	"cdW" -> "cdW(int i1, int i2)",
	"cdB" -> "cdB(int i1, int i2)",
	
	"cHl1" -> "cHl1(int i1, int i2)",
	"cHl3" -> "cHl3(int i1, int i2)",
	"cHe" -> "cHe(int i1, int i2)",
	"cHq1" -> "cHq1(int i1, int i2)",
	"cHq3" -> "cHq3(int i1, int i2)",
	"cHu" -> "cHu(int i1, int i2)",
	"cHd" -> "cHd(int i1, int i2)",
	"cHud" -> "cHud(int i1, int i2)",
	
	"cll" -> "cll(int i1, int i2, int i3, int i4)",
	"cqq1" -> "cqq1(int i1, int i2, int i3, int i4)",
	"cqq3" -> "cqq3(int i1, int i2, int i3, int i4)",
	"clq1" -> "clq1(int i1, int i2, int i3, int i4)",
	"clq3" -> "clq3(int i1, int i2, int i3, int i4)",
	"cee" -> "cee(int i1, int i2, int i3, int i4)",
	"cuu" -> "cuu(int i1, int i2, int i3, int i4)",
	"cdd" -> "cdd(int i1, int i2, int i3, int i4)",
	"ceu" -> "ceu(int i1, int i2, int i3, int i4)",
	"ced" -> "ced(int i1, int i2, int i3, int i4)",
	"cud1" -> "cud1(int i1, int i2, int i3, int i4)",
	"cud8" -> "cud8(int i1, int i2, int i3, int i4)",
	"cle" -> "cle(int i1, int i2, int i3, int i4)",
	"clu" -> "clu(int i1, int i2, int i3, int i4)",
	"cld" -> "cld(int i1, int i2, int i3, int i4)",
	"cqe" -> "cqe(int i1, int i2, int i3, int i4)",
	"cqu1" -> "cqu1(int i1, int i2, int i3, int i4)",
	"cqu8" -> "cqu8(int i1, int i2, int i3, int i4)",
	"cqd1" -> "cqd1(int i1, int i2, int i3, int i4)",
	"cqd8" -> "cqd8(int i1, int i2, int i3, int i4)",
	
	"cledq" -> "cledq(int i1, int i2, int i3, int i4)",
	"cquqd1" -> "cquqd1(int i1, int i2, int i3, int i4)",
	"cquqd8" -> "cquqd8(int i1, int i2, int i3, int i4)",
	"clequ1" -> "clequ1(int i1, int i2, int i3, int i4)",
	"clequ3" -> "clequ3(int i1, int i2, int i3, int i4)" ,
	
	"cduq" -> "cduq(int i1, int i2, int i3, int i4)",
	"cqqu" -> "cqqu(int i1, int i2, int i3, int i4)",
	"cqqq" -> "cqqq(int i1, int i2, int i3, int i4)",
	"cduu" -> "cduu(int i1, int i2, int i3, int i4)"
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
