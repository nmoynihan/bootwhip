(* ::Package:: *)

BeginPackage["CHYTools`"]
PTFactor::usage = "PTFactor[n] outputs the Parke-Taylor factor up to point n: (z[1] - z[2])(z[2]-z[3])...(z[n]-z[1])";
RemoveRowsColumns::usage = "RemoveRowsColumns[M, r, c] removes rows r and columns c from matrix M. For example, to remove rows 1 and 2 and columns 3,5: RemoveRowsColumns[M, {1,2},{3,5}]";
GluonMatrix::usage = "BuildPsi[n] generates the (singular) 2n by 2n matrix of kinematic and polarisation factors for n gluons";
BuildPhi::usage = "BuildPhi[n] builds the matrix \[CapitalPhi]";
BuildPT::usage = "BuildPT[particleLabels] builds the Parke-Taylor factor specifically out of the specified labels. For example: BuildPT[{1,3,5}] = (z[1]-z[3])(z[3]-z[5])(z[5]-z[1])";
BuildPhiPrime::usage = "BuildPhiPrime[n,r,c] builds the Jacobian factor det`\[CapitalPhi]";
Pf::usage = "Pf[A] gives the pfaffian of the skew symmetric A.";
Pf2::usage = "Pf2[A] gives the (alternate) pfaffian of the skew symmetric A.";
FactorByVariable::usage = "FactorByVariable[expr,f] pulls the factor f out of expr to give f*\!\(\*FractionBox[\(expr\), \(f\)]\)";


PTFactor[m_] := Product[(z[i]-z[i+1]),{i,1,m-1}]*(z[m]-z[1]);
(* Build a Pure Gluon (or graviton) Matrix *)
GluonMatrix[num_] := Module[{M,MA,MB,MC,i,j},
(* Construct Matrices *)
\[CapitalDelta][i_,j_] := If[i == 1 && j == num, m^2,0] + If[i == num && j == 1, -m^2,0];
a[i_,j_] := If[i == j,0,If[i < j,2 Dot[k[i],k[j]]/(z[i]-z[j]), 2 Dot[k[j],k[i]]/(z[i]-z[j])] + \[CapitalDelta][i,j]];
b[i_,j_] := If[i == j,0,If[i < j,2 Dot[\[Epsilon][i],\[Epsilon][j]]/(z[i]-z[j]), 2 Dot[\[Epsilon][j],\[Epsilon][i]]/(z[i]-z[j])]];
c[i_,j_] := If[i == j,-Sum[2 Dot[\[Epsilon][i],k[l]]/(z[i]-z[l]),{l,Delete[Table[k,{k,1,num}],i]}],2 Dot[\[Epsilon][i],k[l]]/(z[i]-z[j])];
MA = Table[a[i,j] , {i,1,num},{j,1,num}];
MB= Table[b[i,j] , {i,1,num},{j,1,num}];
MC = Table[c[i,j] , {i,1,num},{j,1,num}];
(* Build Skew-Symmetric Matrix *)
M = ArrayFlatten[{{MA,-1*Transpose[MC]},{MC,MB}}];
If[AntisymmetricMatrixQ[M] == True, M, Message["Operation Failed, the resulting matrix is not skew-symmetric."]; M]
];
(* Build the Matrix *)
BuildPsi[num_] := Module[{M,MA,MB,MC,i,j},
(* Construct Matrices *)
\[CapitalDelta][i_,j_] := If[i == 1 && j == num, m^2,0] + If[i == num && j == 1, -m^2,0];
a[i_,j_] := If[i == j,0,If[i < j,2 Dot[Subscript[k, i],Subscript[k, j]]/(Subscript[z, i]-Subscript[z, j]), 2 Dot[Subscript[k, j],Subscript[k, i]]/(Subscript[z, i]-Subscript[z, j])] + \[CapitalDelta][i,j]];
b[i_,j_] := If[i == j,0,If[i < j,2 Dot[Subscript[\[Epsilon], i],Subscript[\[Epsilon], j]]/(Subscript[z, i]-Subscript[z, j]), 2 Dot[Subscript[\[Epsilon], j],Subscript[\[Epsilon], i]]/(Subscript[z, i]-Subscript[z, j])]];
c[i_,j_] := If[i == j,-Sum[2 Dot[Subscript[\[Epsilon], i],Subscript[k, l]]/(Subscript[z, i]-Subscript[z, l]),{l,Delete[Table[k,{k,1,num}],i]}],2 Dot[Subscript[\[Epsilon], i],Subscript[k, l]]/(Subscript[z, i]-Subscript[z, j])];
MA = Table[a[i,j] , {i,1,num},{j,1,num}];
MB= Table[b[i,j] , {i,1,num},{j,1,num}];
MC = Table[c[i,j] , {i,1,num},{j,1,num}];
(* Build Skew-Symmetric Matrix *)
M = ArrayFlatten[{{MA,-1*Transpose[MC]},{MC,MB}}];
If[AntisymmetricMatrixQ[M] == True, M, Print["Operation Failed."]; M]
];
(* GluonScalarMatrix does not yet function. Should be a 2(n-m) matrix (because the scalar rows should not be present). Build later *)
GluonScalarMatrix[n_,m_] := Module[{M,MA,MB,MC,i,j},
(* Construct Matrices *)
\[CapitalDelta][i_,j_] := If[i == 1 && j == num, m^2,0] + If[i == num && j == 1, -m^2,0];
a[i_,j_] := If[i == j,0,If[i < j,Dot[k[i],k[j]]/(z[i]-z[j]), Dot[k[j],k[i]]/(z[i]-z[j])] + \[CapitalDelta][i,j]];
b[i_,j_] := If[i == j,0,If[i < j,Dot[\[Epsilon][i],\[Epsilon][j]]/(z[i]-z[j]), Dot[\[Epsilon][j],\[Epsilon][i]]/(z[i]-z[j])]];
c[i_,j_] := If[i == j,-Sum[Dot[\[Epsilon][i],k[l]]/(z[i]-z[l]),{l,Delete[Table[k,{k,1,num}],i]}],Dot[\[Epsilon][i],k[l]]/(z[i]-z[j])];
MA = Table[a[i,j] , {i,1,num},{j,1,num}];
MB= Table[b[i,j] , {i,1,num},{j,1,num}];
MC = Table[c[i,j] , {i,1,num},{j,1,num}];
(* Build Skew-Symmetric Matrix *)
M = ArrayFlatten[{{MA,-1*Transpose[MC]},{MC,MB}}];
If[AntisymmetricMatrixQ[M] == True, M, Print["Operation Failed."]; M]
];

BuildPhi[num_] := Module[{M,MA,MB,MC,i,j},
(* Construct Matrices *)
\[CapitalDelta][i_,j_] := If[i == 1 && j == num, m^2,0] + If[i == num && j == 1, -m^2,0];
p[i_,j_] := If[i == j,-Sum[2 (Dot[k[i],k[l]]+\[CapitalDelta][i,l])/(z[i]-z[l])^2,{l,Delete[Table[k,{k,1,num}],i]}],If[i < j,2 (Dot[k[i],k[j]]+\[CapitalDelta][i,j])/(z[i]-z[j])^2, 2 (Dot[k[j],k[i]]+\[CapitalDelta][j,i])/(z[i]-z[j])^2]];
\[CapitalPhi] = Table[p[i,j] , {i,1,num},{j,1,num}];
\[CapitalPhi]
];
BuildPT[array_] := Product[(z[array[[i]]]-z[array[[i+1]]]),{i,1,Length[array] - 1}]*(z[array[[Length[array]]]]-z[array[[1]]]);

BuildPhiPrime[num_, rows_, cols_] := Det[RemoveRowsColumns[BuildPhi[num], rows, cols]] /(BuildPT[rows]BuildPT[cols]);
Begin["`Private`"]
RemoveRowsColumns[m_,rows_,cols_]:=Module[{r,c},
r=Complement[Range@Length@m,rows];
c=Complement[Range@Length@First@m,cols];
m[[r]][[All,c]]
];

FactorByVariable[p_,c_] := c Expand[p/c];

Pf[A_] := Switch[Length[A], 0, 1, _?OddQ, 0, _?EvenQ, xPf[A, 1]];
xPf[A_, p0_] := Module[{A0, n, pivot, sign = 1, A1, p1},
   n = Length[A]/2;
   If[n != 1, A0 = A;
    pivot = First[Ordering[Normal[Abs[A0[[2 n - 1, All]]]], -1]];
    If[pivot != 2 n, A0[[{pivot, 2 n}, All]] = A0[[{2 n, pivot}, All]];
     A0[[All, {pivot, 2 n}]] = A0[[All, {2 n, pivot}]];
     sign = -1;];
    p1 = A0[[2 n - 1, 2 n]];
    A1 = p1 A0[[1 ;; 2 n - 2, 1 ;; 2 n - 2]];
    A1 += (# - Transpose[#]) &@
     Outer[Times, A0[[1 ;; 2 n - 2, 2 n]], 
      A0[[1 ;; 2 n - 2, 2 n - 1]]];
    A1 /= p0;
    sign xPf[A1, p1], A[[1, 2]]]
];

Pf2[Mat_]:=Module[{A, N, ip, pfaff},
A=Mat;
N=Dimensions[A][[1]];
If[OddQ[N], Return[0]];
pfaff=1;
For[i=1, i<N-1, i+=2,
(*find out the maximum entry in the column i, starting from row i+1*)
ip=i+Position[Abs[A[[i+1;;, i]]], Max[Abs[A[[i+1;;,i]]]]][[1,1]];
(*if the maximum entry is not at i+1, permute the matrix so that it is*)
If[i+1 != ip, 
(*Interchange rows and columns in A*)
A[[{i+1, ip},;;]]=A[[{ip,i+1},;;]]; 
A[[;;,{i+1, ip}]]=A[[;;,{ip,i+1}]];
(*interchange contributes det(P)=-1*)
pfaff=-pfaff;
];
(*Multiply with every other entry on the diagonal*)
pfaff = pfaff*A[[i,i+1]];
(*Build the Gauss vector*)
A[[i+2;;,i]]=A[[i+2;;,i]]/A[[i+1,i]];
(*Update the remainder of the matrix using an outer product skew-symmetric update. Note that
column and row i+1 are not affected by the update*)
A[[i+2;; , i+2;; ]]+=(#-Transpose[#])&@Outer[Times,A[[i+2;;,i]],A[[i+2;;,i+1]]]
(* The above is much faster than this construct for me: 
Transpose[{A[[i+2;;,i]]}].{A[[i+2;;,i+1]]}-Transpose[{A[[i+2;;,i+1]]}].{ A[[i+2;;,i]]};*)
];
Return[pfaff*A[[N-1,N]]]
];
End[ ]
EndPackage[ ]









