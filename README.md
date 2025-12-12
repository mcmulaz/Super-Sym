# Super-Sym
Contains mathematica code for L and R matrices and option to test the closure of the GRDN algebra. 

The file "Garden Algebra" contains the definitions of the L and R matrices. Below are commands that allow one to see the explicit entries in each L and R matrix.  Commands are given to test the anticommutator relations, and to see the entries of the holoraumy tensor associated with the commutator equations. 

File "Garden Algebra" contains the code that runs the verification of the garden algebra. It also has the definitions of the L and R matrices. 
The other algebra RL + RL = delta + (non-closure) is also given. Explicit representation of any left (right) matrix is given by L[1] (R[1]) // Matrix Form (replace 1 with whichever index). 

Numerical Representation of Matrices (Replace 1 with desired matrix): 
L[1] // Matrix Form
R[1] // Matrix Form 

Explicit Holoraumy Representation (Replace 1 and 2 with desired Holoraumy matrix): 
MatrixForm[ (R[1].L[2] - R[2].L[1] ]
MatrixForm[ (L[1].R[2] - L[2].R[1] ]

Explicit form of non-closure term:
MatrixForm[ (R[1].L[2] + R[2].L[1] - 2IdentityMatrix[176])/2 ]

Closure Relation Check for LR + LR = 2I: 
(* =========Relation checker:L_i R_j+L_j R_i=2 \[Delta]_{ij} I=========*)
Clear[Lmat, Rmat];
Lmat[a_Integer] := Lmat[a] = L[a];
Rmat[a_Integer] := Rmat[a] = R[a];

(*Single-pair check;returns {i,j,exactTrueQ,frobError}*)
checkPair[i_Integer, j_Integer] := 
  Module[{lhs, rhs, diff}, lhs = Lmat[i] . Rmat[j] + Lmat[j] . Rmat[i];
   rhs = 2 KroneckerDelta[i, j] IdentityMatrix[colPhi];
   diff = Simplify[lhs - rhs];
   {i, j, diff === ConstantArray[0, {colPhi, colPhi}], 
    Norm[Chop@N@diff, "Frobenius"]}];

(*Sweep all i,j\[Element]{1..16}*)
results = 
  Table[checkPair[i, j], {i, 1, 16}, {j, 1, 16}] // Flatten[#, 1] &;

(*Quick summaries*)
numOK = Count[results[[All, 3]], True];
numTot = Length[results];
maxErr = Max[results[[All, 4]]];

Print["Pairs passing exactly: ", numOK, " / ", numTot];
Print["Max Frobenius error (numeric, after Chop): ", maxErr];

(*List any failures with their numeric error*)
failures = Select[results, Not@#[[3]] &];
If[failures === {}, 
 Print["All pairs satisfy L_i R_j + L_j R_i = 2 \[Delta]_{ij} I."], 
 Print["Violations (i,j, exact?, FrobeniusError):"];
 failures]


 
