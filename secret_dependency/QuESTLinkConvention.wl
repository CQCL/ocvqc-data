(* ::Package:: *)

Stacking::usage = "Pick row or colum stacking. By default, it's column stacking";
Vectorize::usage = "Vectorization of a densitry matrix, Option: Stacking -> column or row"
Options[Tensorize] = {Stacking -> "column"};
Tensorize::usage = "Turns a vector into 2x2 matrix. By default, it's column stacking";


Options[Vectorize] = {Stacking -> "column"};
Vectorize[m_, OptionsPattern[]] := ArrayReshape[
	If[OptionValue[Stacking] == "row",
		TensorTranspose[m , {2, 1}],
	m]
,
 {Times @@ Dimensions @ m}]
 
 Tensorize[v_, OptionsPattern[]] := With[ {d = Round @ Sqrt @ First @ Times @ Dimensions @ v},
	 If[ OptionValue @ Stacking === "row", 
		 TensorTranspose[ArrayReshape[v, {d, d}], {2, 1}]
		 ,
		 ArrayReshape[v, {d, d}]
	 ]
 ]
 
(* Superopeators-related*)
bpSWAP::usage="bpSWAP[matrix]. Bipartite swap for 4x4 matrix: ";
rowShuffle::usage="Row-shuffling in LSB";
columnShuffle::usage="Column-shuffling in LSB";
swap::usage="swap[vector]. SWAP the vector to convert column \[DoubleLeftRightArrow] row stacking convention.";
bitReorder::usage="bitReorder[matrix]. Reorder the LSB \[DoubleLeftRightArrow] MSB of a square matrix.";

bpSWAP[m_] /; (Dimensions[m] == {4,4}):=Module[{newm, mt},
	mt = ArrayReshape[m, {2, 2, 2, 2}]; (* i, j -- u -- k, l*)
	newm = TensorTranspose[m, {2, 1, 4, 3}]; (* j, i -- u -- l, j  *)
	(* shuffle back *)
	ArrayReshape[newm, {4, 4}]
]

rowShuffle[m_] := With[{mt = ArrayReshape[m, {2, 2, 2, 2}]}, 
	ArrayReshape[TensorTranspose[mt, {4, 2, 3, 1}], {4, 4}]
]

columnShuffle[m_] := With[{mt = ArrayReshape[m, {2, 2, 2, 2}]},
	ArrayReshape[TensorTranspose[mt, {1, 3, 2, 4}], {4, 4}]
]

swap[vec_] := With[{d = Round @ Sqrt @ First @ Times @ Dimensions @ vec}, 
	ArrayReshape[TensorTranspose[ArrayReshape[vec, {d, d}], {2, 1}], d*d]
]

bitReorder[m_] := Module[{ dims, newmat, newdims }, 
	dims = Dimensions @ m;
	newdims = Join @@ (ConstantArray[2, Log2 @ #]& /@ dims);
	newmat = ArrayReshape[m, newdims];
	ArrayReshape[TensorTranspose[newmat, Reverse @ Range[Length @ newdims]], dims]
] 
 
(*Choi-matrix representations*)
choiMatrix::usage = "choiMatrix[kraus_list] Returns choi matrix with LSB: (\!\(\*SubscriptBox[\(I\), \(0\)]\) \[TensorProduct] \!\(\*SubscriptBox[\(E\), \(1\)]\))(\[CapitalPhi]) for column and (\!\(\*SubscriptBox[\(E\), \(0\)]\) \[TensorProduct] \!\(\*SubscriptBox[\(I\), \(1\)]\))(\[CapitalPhi]) for row -stacking";
superoperatorMatrix::usage = "superoperatorMatrix[kraus_list]. Returns superoperator \[Sum] \!\(\*SuperscriptBox[SubscriptBox[\(K\), \(0\)], \(*\)]\)\[TensorProduct] \!\(\*SubscriptBox[\(K\), \(1\)]\) for column-stacking, \[Sum] \!\(\*SuperscriptBox[SubscriptBox[\(K\), \(1\)], \(*\)]\)\[TensorProduct] \!\(\*SubscriptBox[\(K\), \(0\)]\) for row-stacking, with LSB convention.";
choiOperate::usage = "choiOperate[choi_matrix, density_matrix]. Perform choi operation. ";


Options[choiMatrix] = {Stacking -> "column"};
choiMatrix[kraus_List, OptionsPattern[]] := With[{\[Phi]={{1,0,0,1}, {0,0,0,0}, {0,0,0,0}, {1,0,0,1}} (* |00> + |11> *)},
	If[ OptionValue @ Stacking == "row",
		Total @ Table[KroneckerProduct[IdentityMatrix@2, k] . \[Phi] . KroneckerProduct[IdentityMatrix@2, ConjugateTranspose@k], {k, kraus }]
		,
		Total @ Table[KroneckerProduct[k, IdentityMatrix@2] . \[Phi] . KroneckerProduct[ConjugateTranspose@k, IdentityMatrix@2], {k, kraus }]
	]
]

Options[superoperatorMatrix] = {Stacking -> "column"};
superoperatorMatrix[kraus_List, OptionsPattern[]] :=
If[ OptionValue @ Stacking == "row",
	Total @ Table[ KroneckerProduct[Conjugate[k], k] , {k, kraus}]
,
	Total @ Table[ KroneckerProduct[k, Conjugate[k]] , {k, kraus}]
]

choiOperate[\[CapitalLambda]_, \[Rho]_] := TensorContract[
		TensorProduct[ArrayReshape[\[CapitalLambda], {2, 2, 2, 2}], \[Rho]]
	, {{2, 5}, {4, 6}}]
