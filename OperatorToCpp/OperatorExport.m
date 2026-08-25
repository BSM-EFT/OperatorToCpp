(* ::Package:: *)

(* ::Chapter:: *)
(*-- Package Context --*)


BeginPackage["OperatorExport`", {"Global`"}];


(* ::Subsection:: *)
(*Description of public functions*)


SimplifyOutput::usage = "SimplifyOutput[matchedResult,WCInfo] creates a dictionary containing only Wilson coefficients and the corresponding matching conditions, while unwraping the Matchete API and converting the expressions to a lighter intermediate form.";
BuildFiles::usage = "BuildFiles[ModelName,SimplifiedOutput,ComplexParams,WCInfo] offers a single command to build C++ header (.h) and source (.cpp) files, and a `pyBindings.cpp` file to generate Python bindings, along with Python declarations of the Model class and its methods.";


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
      "OperatorExport  v1.0.0\[Beta] \[LongDash] by ",
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


createMatchingDict[Output_,WCInfo_] := Module[{d6Dict, keyList, key1, key2},
	d6Dict= Association[];
	keyList = Keys[Output];
	
	Do[
		key1=StringCases[ToString[keyList[[i]]], RegularExpression["(?<=Coupling\\[c)"]~~x:WordCharacter..:>"c"<>x];
		(* this assumption is crucial - Matchete produced Wilson Coefficients have names starting with "c" *)
		If[Length[key1]==0,Continue[], AppendTo[d6Dict,key1[[1]]->Values[Output][[i]]]],
	{i,Length[keyList]}];
	
	Do[
		key2=Keys[WCInfo][[i]];
		If[KeyExistsQ[d6Dict,key2],
			Continue[],
			AppendTo[d6Dict,key2->0.0]],
		{i,Length[WCInfo]}];
		
	Return[d6Dict];
]


(* ::Subsection:: *)
(*Simplify the output through suitable replacement of symbols*)


SimplifyOutput[output_,WCInfo_]:=Module[{dict,newDict},
	dict = createMatchingDict[output,WCInfo];
	dict = dict /. {Index[i1,Flavor]->i1,Index[i2,Flavor]->i2,Index[i3,Flavor]->i3,Index[i4,Flavor]->i4 };
	dict = dict /. {Index[x_,__]:>ToExpression[StringJoin["r",StringPart[ToString[x],4;;-1]]]};
	dict = dict /. Bar[Coupling[x_,{y___},z_]] :> Coupling[ToString[x]<>"_c",{y},z];
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
		Do[AppendTo[rules,ToString[list[[i]]]<>"_c"->list[[i]]],{i,1,Length[list]}];
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


(*RulesFile=FileNameJoin[{NotebookDirectory[],"./LFRules.m"}];
ReplLoopFuncExpns=Import[RulesFile];*)


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
		Return[rest*Times@@rel]
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
	ReplMassVecinTF[expr/.$ReplLoopFuncExpns/.WrapLFinTF]/.ExtractFIdxfromTF
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
	If[Head[expr]===EinsSum,Return[StringReplace[ToString[expr],{"["->"(","]"->")"}]]];
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
(*Helper function to generate WC method declarations*)


CppFnDeclaration[wc_,WCInfo_]:=Module[{nflav, decl},
	nflav=WCInfo[wc]["Nf"];
	decl=wc<>"(";
	If[nflav>0,
		Do[decl=StringJoin[decl,"int i",ToString[i],", "],{i,1,nflav}];
		decl=StringDrop[decl,-2];
	];
	decl=decl<>")";
	Return[decl];
]


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

HeaderModelClass[className_,paramList_,ComplexPars_,WCInfo_,line_] := Module[{},
	WriteLine[line, "class "<>className<>" {"];
	WriteLine[line, "    private:"];
	
	(* define and initialize the scale and loop-order *)
	WriteLine[line, "        double hbar = 1/(16 * pow(pi,2));"];
	WriteLine[line, "        double mubarsq = 1000.0;"];
	(* define and initialize 0-dimensional parameters *)
	If[Length[paramList[[1]]]!=0,
		Do[
			WriteLine[line, "        "<>"std::complex<double> "<>ToString[paramList[[1]][[i]]]<>" = 0.0;"];
			If[
				MemberQ[ComplexPars,paramList[[1]][[i]]],
				WriteLine[line, "        "<>"std::complex<double> "<>ToString[paramList[[1]][[i]]]<>"_c = 0.0;"];
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
				WriteLine[line, "        "<>"std::vector<std::complex<double> > "<>ToString[paramList[[2]][[i]]]<>"_c = {0.0, 0.0, 0.0};"];
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
				WriteLine[line, "        "<>"std::vector<std::vector<std::complex<double> > > "<>ToString[paramList[[3]][[i]]]<>"_c = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};"];
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
	Do[WriteLine[line,"        std::complex<double> "<>CppFnDeclaration[Keys[WCInfo][[i]],WCInfo]<>";"],{i,1,Length[Keys[WCInfo]]}];	
	WriteLine[line, ""];
	
	(* declaration for the batch evaluator *)
	WriteLine[line, "        std::map<std::string, std::complex<double> > batch_eval(const std::vector<Task>& tasks);"];
	WriteLine[line, "};"];
	WriteLine[line, ""];
];


(* ::Subsubsection:: *)
(*Master builder*)


HeaderFileBuilder[modelName_,paramList_,ComplexPars_,WCInfo_]:=Module[{path,fullPath,line1},
	path = FileNameJoin[{NotebookDirectory[],"include"}];
	If[!DirectoryQ[path], CreateDirectory[path]];
	fullPath = FileNameJoin[{path,modelName<>".h"}];
	line1 = OpenWrite[fullPath];
	HeaderPreprocessorDirectives[line1];
	WriteLine[line1,""];
	HeaderTaskStruct[line1];
	WriteLine[line1,""];
	HeaderModelClass[modelName,paramList,ComplexPars,WCInfo,line1];
	Close[line1];
];



(* ::Subsection:: *)
(*Build the source (.cpp) file*)


(* ::Subsubsection:: *)
(*Builder for preprocessor directives, constructor and updater methods*)


BuildPreprocessorDirectives[modelName_,line_] := Module[{},
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
				WriteLine[line, "        this->"<>ToString[paramList[[1]][[i]]]<>"_c = std::conj(params[\""<>ToString[paramList[[1]][[i]]]<>"\"]);"];];
			WriteLine[line, "    }"], {i,1,Length[paramList[[1]]]}]	
	];
	
	If[Length[paramList[[2]]]!=0,
		Do[
			Do[
				WriteLine[line, "    if (params.find(\""<>ToString[paramList[[2]][[i]]]<>ToString[j]<>"\") != params.end()) {"];
				WriteLine[line, "        this->"<>ToString[paramList[[2]][[i]]]<>"["<>ToString[j-1]<>"] = params[\""<>ToString[paramList[[2]][[i]]]<>ToString[j]<>"\"];"];
				If[MemberQ[ComplexPars,paramList[[2]][[i]]],
					WriteLine[line, "        this->"<>ToString[paramList[[2]][[i]]]<>"_c["<>ToString[j-1]<>"] = std::conj(params[\""<>ToString[paramList[[2]][[i]]]<>ToString[j]<>"\"]);"];
				];
				WriteLine[line, "    }"];,
			{j,1,3}], 
		{i,1,Length[paramList[[2]]]}]
	];
	
	If[Length[paramList[[3]]]!=0,
		Do[
			Do[
				WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>ToString[j]<>ToString[k]<>"\") != params.end()) {"];
				WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"["<>ToString[j-1]<>"]["<>ToString[k-1]<>"] = params[\""<>ToString[paramList[[3]][[i]]]<>ToString[j]<>ToString[k]<>"\"];"];
				If[MemberQ[ComplexPars,paramList[[3]][[i]]],
					WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"_c["<>ToString[k-1]<>"]["<>ToString[j-1]<>"] = std::conj(params[\""<>ToString[paramList[[3]][[i]]]<>ToString[j]<>ToString[k]<>"\"]);"];
				];
				WriteLine[line, "    }"];,
			{j,1,3},{k,1,3}],
		{i,1,Length[paramList[[3]]]}]	
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
				WriteLine[line, "        this->"<>ToString[paramList[[1]][[i]]]<>"_c = std::conj(params[\""<>ToString[paramList[[1]][[i]]]<>"\"]);"];];
			WriteLine[line, "    }"], {i,1,Length[paramList[[1]]]}]	
	];
	
	If[Length[paramList[[2]]]!=0,
		Do[
			Do[
				WriteLine[line, "    if (params.find(\""<>ToString[paramList[[2]][[i]]]<>ToString[j]<>"\") != params.end()) {"];
				WriteLine[line, "        this->"<>ToString[paramList[[2]][[i]]]<>"["<>ToString[j-1]<>"] = params[\""<>ToString[paramList[[2]][[i]]]<>ToString[j]<>"\"];"];
				If[MemberQ[ComplexPars,paramList[[2]][[i]]],
					WriteLine[line, "        this->"<>ToString[paramList[[2]][[i]]]<>"_c["<>ToString[j-1]<>"] = std::conj(params[\""<>ToString[paramList[[2]][[i]]]<>ToString[j]<>"\"]);"];
				];
				WriteLine[line, "    }"];,
			{j,1,3}], 
		{i,1,Length[paramList[[2]]]}]
	];
	
	If[Length[paramList[[3]]]!=0,
		Do[
			Do[
				WriteLine[line, "    if (params.find(\""<>ToString[paramList[[3]][[i]]]<>ToString[j]<>ToString[k]<>"\") != params.end()) {"];
				WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"["<>ToString[j-1]<>"]["<>ToString[k-1]<>"] = params[\""<>ToString[paramList[[3]][[i]]]<>ToString[j]<>ToString[k]<>"\"];"];
				If[MemberQ[ComplexPars,paramList[[3]][[i]]],
					WriteLine[line, "        this->"<>ToString[paramList[[3]][[i]]]<>"_c["<>ToString[k-1]<>"]["<>ToString[j-1]<>"] = std::conj(params[\""<>ToString[paramList[[3]][[i]]]<>ToString[j]<>ToString[k]<>"\"]);"];
				];
				WriteLine[line, "    }"];,
			{j,1,3},{k,1,3}],
		{i,1,Length[paramList[[3]]]}]	
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


BuildFunction[modelName_,WCname_,expr_,ComplexPars_,WCInfo_]:=Module[{returnExpr,path,fullPath,line},
	returnExpr = ConvertFullExpression[expr,ComplexPars];
	path = FileNameJoin[{NotebookDirectory[],"lib",ToString[WCInfo[WCname]["Nf"]]<>"f"}];
	If[!DirectoryQ[path], CreateDirectory[path]];
	fullPath = FileNameJoin[{path,WCname<>".cpp"}];
	line = OpenWrite[fullPath];
		
	WriteLine[line, "#include \"OperatorImport.h\""];
	WriteLine[line, "#include \"complex_math.h\""];
	WriteLine[line, "#include \""<>modelName<>".h\""];
	WriteLine[line,""];
	WriteLine[line, "std::complex<double> "<>modelName<>"::"<>CppFnDeclaration[WCname,WCInfo]<> " {"];		
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


SourceFileBuilder[modelName_, paramList_, ComplexPars_, matchingOutput_,WCInfo_]:=Module[{keyList,exprList,path,fullPath,line1},
	keyList=Keys[matchingOutput];
	exprList=Values[matchingOutput];
	
	path = FileNameJoin[{NotebookDirectory[],"lib"}];
	If[!DirectoryQ[path], CreateDirectory[path]];
	fullPath = FileNameJoin[{path,modelName<>".cpp"}];
	line1 = OpenWrite[fullPath];
	
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
		BuildFunction[modelName,keyList[[k]],exprList[[k]],ComplexPars,WCInfo],
	{k,1,Length[matchingOutput]}];
	
]



(* ::Subsection:: *)
(*Builder for pyBindings.cpp*)


(* ::Input::Initialization:: *)
(* to generate py::arg("i1") py::arg("i1") *)
PyArgs[nflav_]:=Module[{args},
args ="";
Do[args=args<>"py::arg(\"i"<>ToString[i]<>"\"), ",{i,1,nflav}];
args=StringDrop[args,-2];	
Return[args];
]


(* ::Input::Initialization:: *)
(*to generate i1, i2, i3, i4 *)
SimpleArgs[nflav_]:=Module[{args},
args ="";
Do[args=args<>"i"<>ToString[i]<>", ",{i,1,nflav}];
args=StringDrop[args,-2];	
Return[args];
]


(* ::Input::Initialization:: *)
(*to generate int i1, int i2, int i3, int i4 *)
TypedArgs[nflav_]:=Module[{args},
args ="";
Do[args=args<>"int i"<>ToString[i]<>", ",{i,1,nflav}];
args=StringDrop[args,-2];	
Return[args];
]


(* ::Input::Initialization:: *)
AddWrapperLambdas[modelName_,WCInfo_,line_]:=Module[{NfUniq,wcs,wc,nf,lambdaArg},  
NfUniq=DeleteDuplicates[Table[WCInfo[[i]]["Nf"],{i,1,Length[WCInfo]}]];
wcs=Table[Select[WCInfo,#["Nf"]===NfUniq[[j]]&],{j,1,Length[NfUniq]}];
WriteLine[line, "        // wrapper lambdas"];
Do[
nf=NfUniq[[k]];
lambdaArg="[]("<>modelName<>"& self, std::string name, std::string method_name"<>If[nf>0,", "<>TypedArgs[nf],""]<>")";
WriteLine[line, "        .def(\"wrap_"<>ToString[nf]<>"f\", " <>lambdaArg<>" {" ];
WriteLine[line, "            static const std::unordered_map<std::string, std::complex<double> ("<>modelName<>"::*)("<>If[nf>0,StringDrop[StringRepeat["int, ",nf],-2],""]<>")> map_"<>ToString[nf]<>"f = {"];
Do[
wc=Keys[wcs[[k]]][[l]];
WriteLine[line, "                {\""<>wc<>"\", &"<>modelName<>"::"<>wc<>"},"],{l,1,Length[wcs[[k]]]}];
WriteLine[line, "            };"];
WriteLine[line, "            auto it = map_"<>ToString[nf]<>"f.find(method_name);"];
WriteLine[line, "            if (it == map_"<>ToString[nf]<>"f.end()) throw std::runtime_error(\"Method not found: \" + method_name);"];
WriteLine[line, ""];
WriteLine[line, "            auto fn_ptr = it->second;"];
WriteLine[line, "            return Task{name, [&self, fn_ptr"<>If[nf>0,", "<>SimpleArgs[nf],""]<>"](){ return (self.*fn_ptr)("<>If[nf>0,SimpleArgs[nf],""]<>"); }};"];
WriteLine[line, "        }, py::arg(\"name\"), py::arg(\"method_name\")"<>If[nf>0,", "<>PyArgs[nf],""]<>")"];
WriteLine[line, ""];,{k,1,Length[NfUniq]}];
]


(* ::Input::Initialization:: *)
AddWCMethods[modelName_,WCInfo_,line_]:=Module[{wc,wcLast},
WriteLine[line, "        // Wilson coefficient methods"];
Do[
wc=Keys[WCInfo][[i]];
WriteLine[line, "        .def(\""<>wc<>"\", &"<>modelName<>"::"<>wc<>If[WCInfo[wc]["Nf"]>0, ", "<>PyArgs[WCInfo[wc]["Nf"]], ""]<>", py::call_guard<py::gil_scoped_release>())"];,{i,1,Length[Keys[WCInfo]]-1}];
wcLast=Keys[WCInfo][[-1]];
WriteLine[line, "        .def(\""<>wcLast<>"\", &"<>modelName<>"::"<>wcLast<>If[WCInfo[wcLast]["Nf"]>0, ", "<>PyArgs[WCInfo[wcLast]["Nf"]], ""]<>", py::call_guard<py::gil_scoped_release>());"];
]


(* ::Input::Initialization:: *)
GeneratePyBindings[modelName_,WCInfo_]:=Module[{path,fullPath,line1},
path = FileNameJoin[{NotebookDirectory[],"src"}];
If[!DirectoryQ[path], CreateDirectory[path]];
fullPath = FileNameJoin[{path,"pyBindings.cpp"}];
line1 = OpenWrite[fullPath];
WriteLine[line1, "#include <pybind11/pybind11.h>"];
WriteLine[line1, "#include <pybind11/cast.h>"];
WriteLine[line1, "#include <pybind11/detail/common.h>"];
WriteLine[line1, "#include <pybind11/stl.h>"];
WriteLine[line1, "#include <pybind11/complex.h>"];
WriteLine[line1, "#include <pybind11/functional.h>"];
WriteLine[line1, "#include <pybind11/pytypes.h>"];
WriteLine[line1, "#include <stdexcept>"];
WriteLine[line1, "#include <unordered_map>"];
WriteLine[line1, "#include <string>"];
WriteLine[line1, "#include <complex>"];
WriteLine[line1, "#include \""<>modelName<>".h\""];
WriteLine[line1, ""];
WriteLine[line1, "namespace py = pybind11;"];
WriteLine[line1, ""];
WriteLine[line1, "PYBIND11_MODULE(match_to_py, m, py::mod_gil_not_used()) {"];
WriteLine[line1, "    auto task = py::class_<Task>(m, \"Task\");"];
WriteLine[line1, ""];
WriteLine[line1, "    py::class_<"<>modelName<>">(m, \""<>modelName<>"\")"];
WriteLine[line1, "        .def(py::init<std::unordered_map<std::string, std::complex<double> >, double, bool>())"];
WriteLine[line1, "        .def(\"updateParams\", &"<>modelName<>"::updateParams, py::arg(\"param_dict\"))"];
WriteLine[line1, "        .def(\"getScale\", &"<>modelName<>"::getScale)"];
WriteLine[line1, "        .def(\"setScale\", &"<>modelName<>"::setScale, py::arg(\"scale\"))"];
WriteLine[line1, "        .def(\"loopContributions\", &"<>modelName<>"::loopContributions, py::arg(\"loop\"))"];
WriteLine[line1, "        .def(\"getParams\", &"<>modelName<>"::getParams)"];
WriteLine[line1, ""];
WriteLine[line1, "        // batch evaluator"];
WriteLine[line1, "        .def(\"batch_eval\", &"<>modelName<>"::batch_eval, py::arg(\"tasks\"), py::call_guard<py::gil_scoped_release>())"];
WriteLine[line1, ""];
AddWrapperLambdas[modelName,WCInfo,line1];
AddWCMethods[modelName,WCInfo,line1];
WriteLine[line1, "}"];
Close[line1];
]


(* ::Subsection:: *)
(*Builder for .pyi file with python type information for all WC functions*)


PyFnDeclaration[wc_,WCInfo_]:=Module[{nflav, decl},
	nflav=WCInfo[wc]["Nf"];
	decl="def "<>wc<>"(self, ";
	If[nflav>0,
		Do[decl=StringJoin[decl,"i",ToString[i],": int, "],{i,1,nflav}];
	];
	decl=StringDrop[decl,-2];
	decl=decl<>") -> complex: ...";
	Return[decl];
]

GeneratePythonDeclarations[modelName_,WCInfo_]:= Module[{path,fullPath,line},
	path = FileNameJoin[{NotebookDirectory[],"py"}];
	If[!DirectoryQ[path], CreateDirectory[path]];
	fullPath = FileNameJoin[{path,"match_to_py.pyi"}];
	line = OpenWrite[fullPath];
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
	Do[WriteLine[line, "    "<>PyFnDeclaration[Keys[WCInfo][[i]],WCInfo]],{i,1,Length[WCInfo]}];
	Close[line];
]


(* ::Subsection:: *)
(*Single function to call all file builders*)


BuildFiles[modelName_,simplifiedOutput_,complexPars_,WCInfo_]:=Module[{modelPars},
	modelPars=SegregateParams[simplifiedOutput,complexPars];
	HeaderFileBuilder[modelName,modelPars,complexPars,WCInfo];
	Print["C++ header (.h) files successfully created and placed in the ./include/ subdirectory."];
	SourceFileBuilder[modelName,modelPars,complexPars,simplifiedOutput,WCInfo];
	Print["C++ source (.cpp) files successfully created and placed in the ./lib/ subdirectory."];
	GeneratePyBindings[modelName,WCInfo];
	Print["pyBindings.cpp successfully created and placed in the ./src/ subdirectory"];
	GeneratePythonDeclarations[modelName,WCInfo];
	Print["match_to_py.pyi successfully created and placed in the ./py/ subdirectory"]
]


(* ::Chapter:: *)
(*-- Package Context --*)


End[];
EndPackage[];
