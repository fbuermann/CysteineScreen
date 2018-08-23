(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)

(* :Title: CysteineScreen *)
(* :Context: CysteineScreen` *)
(* :Author: Frank Buermann <fburmann@mrc-lmb.cam.ac.uk> *)

BeginPackage["CysteineScreen`"];
ClearAll["`*"];

(***********************************************************************)
(* :Documentation: *)
(***********************************************************************)
makePrimerPair::usage =
    "makePrimerPair[\"cds\", res, codon] gives PCR primers for\
 mutagenesis of a residue using a specified codon.";

(***********************************************************************)
Begin["`Private`"];

(***********************************************************************)
(**** Configuration ****)

(* Minimum melting temperature *)
melt = Quantity[50, "DegreesCelsius"];

(*
 Codon replacement rules for the residue in front the of the mutation.
 These codons end on T or G, limiting BsaI overhangs to
 TXYZ and GXYZ where XYZ is the specified substitution codon.
*)
precedingCodons =
    {
      "R" -> "CGT",
      "L" -> "CTT",
      "S" -> "TCT",
      "A" -> "GCT",
      "G" -> "GGT",
      "P" -> "CCT",
      "T" -> "ACT",
      "V" -> "GTT",
      "I" -> "ATT",
      "N" -> "AAT",
      "D" -> "GAT",
      "C" -> "TGT",
      "Q" -> "CAG",
      "E" -> "GAG",
      "H" -> "CAT",
      "K" -> "AAG",
      "F" -> "TTT",
      "Y" -> "TAT",
      "M" -> "ATG",
      "W" -> "TGG"
    };

(***********************************************************************)
(**** Utility functions ****)

(* Make reverse complementary of a DNA sequence *)
dnaReverse[seq_String] :=
    StringReplace[
      StringReverse[seq],
      {
        "A" -> "T",
        "C" -> "G",
        "G" -> "C",
        "T" -> "A",
        "a" -> "t",
        "c" -> "g",
        "g" -> "c",
        "t" -> "a"
      }
    ];

translationRules =
    {
      "UCA" -> "S",
      "UCC" -> "S",
      "UCG" -> "S",
      "UCU" -> "S",
      "UUC" -> "F",
      "UUU" -> "F",
      "UUA" -> "L",
      "UUG" -> "L",
      "UAC" -> "Y",
      "UAU" -> "Y",
      "UAA" -> "Z",
      "UAG" -> "Z",
      "UGC" -> "C",
      "UGU" -> "C",
      "UGA" -> "Z",
      "UGG" -> "W",
      "CUA" -> "L",
      "CUC" -> "L",
      "CUG" -> "L",
      "CUU" -> "L",
      "CCA" -> "P",
      "CCC" -> "P",
      "CCG" -> "P",
      "CCU" -> "P",
      "CAC" -> "H",
      "CAU" -> "H",
      "CAA" -> "Q",
      "CAG" -> "Q",
      "CGA" -> "R",
      "CGC" -> "R",
      "CGG" -> "R",
      "CGU" -> "R",
      "AUA" -> "I",
      "AUU" -> "I",
      "AUC" -> "I",
      "AUG" -> "M",
      "ACA" -> "T",
      "ACC" -> "T",
      "ACG" -> "T",
      "ACU" -> "T",
      "AAC" -> "N",
      "AAU" -> "N",
      "AAA" -> "K",
      "AAG" -> "K",
      "AGC" -> "S",
      "AGU" -> "S",
      "AGA" -> "R",
      "AGG" -> "R",
      "GUA" -> "V",
      "GUC" -> "V",
      "GUG" -> "V",
      "GUU" -> "V",
      "GCA" -> "A",
      "GCC" -> "A",
      "GCG" -> "A",
      "GCU" -> "A",
      "GAC" -> "D",
      "GAU" -> "D",
      "GAA" -> "E",
      "GAG" -> "E",
      "GGA" -> "G",
      "GGC" -> "G",
      "GGG" -> "G",
      "GGU" -> "G"
    };

(* Translate a coding sequence. *)
dnaTranslateCDS[seq_String] :=
    StringJoin[
      StringJoin /@
          Partition[Characters[ToUpperCase@seq] /. "T" -> "U", 3]
          /. translationRules
    ];

(* Compute the melting temperature for an oligonucleotide *)
(*
 The nearest-neighbor method with unified parameters from
 Santa Lucia, 1998 (PNAS) is employed. The correction for divalent
 cations from von Ahsen et al., 2001 (Clin Chem) is used.
*)
dnaMeltingTemperature[
  seq_String,
  conc_ : 500, (* Oligo concentration in nM *)
  monovalent_ : 50, (* Monovalent salt in mM *)
  divalent_ : 1.5, (* Divalent salt in mM *)
  dNTPs_ : 1 (* dNTP concentration in mM *)
] :=
    Module[
      {
        oligo, labels, dHmatrix, dSmatrix,
        dinucleotides1, dinucleotides2,
        dn, ends, nnMatrix, dH, dS, dG, tM
      },

      oligo = ToUpperCase@seq;
      labels =
          {
            "AA-TT", "AT-TA", "TA-AT",
            "CA-GT", "GT-CA", "CT-GA",
            "GA-CT", "CG-GC", "GC-CG",
            "GG-CC", "initGC", "initAT"
          };

      (* Enthalpy *)
      dHmatrix =
          {
            -7.9, -7.2, -7.2, -8.5, -8.4, -7.8,
            -8.2, -10.6, -9.8, -8.0, 0.1, 2.3
          };
      (* Entropy *)
      dSmatrix =
          {
            -22.2, -20.4, -21.3, -22.7, -22.4, -21.0,
            -22.2, -27.2, -24.4, -19.9, -2.8, 4.1
          };

      dinucleotides1 =
          {"AA", "AT", "TA", "CA", "GT", "CT", "GA", "CG", "GC", "GG"};
      dinucleotides2 =
          {"TT", "xx", "xx", "TG", "AC", "AG", "TC", "xx", "xx", "CC"};

      (* Count dinucleotide instances. *)
      dn =
          Map[
            StringCount[oligo, #, Overlaps -> True]&,
            {dinucleotides1, dinucleotides2},
            {2}
          ] // Total;

      ends =
          Map[
            StringCount[
              StringTake[oligo, 1] <> StringTake[oligo, -1],
              Alternatives @@ #
            ]&,
            {{"G", "C"}, {"A", "T"}}
          ];

      (* Construct nearest-neighbour matrix. *)
      nnMatrix = Join[dn, ends];

      (* Compute enthalpy (dH) and entropy (dS) contributions,
      binding energy (dG) and melting temperature (tM). *)
      dH = nnMatrix.dHmatrix;
      dS =
          nnMatrix.dSmatrix
              + 0.368 * (StringLength[seq] - 1)
              *
              Log[(monovalent + 120 * Sqrt[divalent - dNTPs]) *
                  10^-3];
      dG = dH - dS * (37 + 273.15) * 10^-3;
      tM = dH / ((dS + 1.987 * Log[conc * 10^-9 / 4]) * 10^-3) - 273.15;

      Quantity[tM, "DegreesCelsius"]
    ];

(* Mutate a codon corresponding to an amino acid position *)
dnaMutateCDS[
  cds_String, (* The CDS *)
  res_Integer, (* Amino acid position to be mutated *)
  codon_String (* Codon for the substitution *)
] :=
    StringReplacePart[cds, codon, {res * 3 - 2, res * 3}];

(***********************************************************************)
(**** Primer construction ****)

(* Construct primers for mutagenesis *)
SyntaxInformation[makePrimerPair] =
    {"ArgumentsPattern" -> {_, _, _, _., _.}};
makePrimerPair[
  cds_String, (* Coding sequence *)
  res_Integer, (* Residue number to be mutated *)
  codon_String, (* Codon for substitution *)
  maxTrim_Integer : 3, (* Homology optimization parameter *)
  maxExtend_Integer : 6 (* Homology optimization parameter *)
] :=
    Module[
      {
        precedingCodon, (* Codon in front of the mutation *)
        precedingSubstitutionFlag, (* Flag for preceding codon *)
        mutatedSeq,
        fwdHomology,
        fwdExtension,
        revHomology,
        revExtension
      }
      ,
      precedingCodon = StringTake[cds, {0, 2} + 3 * res - 5];

      mutatedSeq =
          If[
            MemberQ[
              ToUpperCase@StringTake[precedingCodon, {-1}],
              {"G", "T"}
            ],
          (* Use sequence as is *)
            precedingSubstitutionFlag = False;
            ToLowerCase@cds
            ,
          (* Or make base in front of the mutation a G or T *)
            precedingSubstitutionFlag = True;
            dnaMutateCDS[
              ToLowerCase@cds,
              res - 1,
              dnaTranslateCDS@precedingCodon /. precedingCodons
            ]
          ]
          (* Do substitution *)
              // dnaMutateCDS[#, res, ToUpperCase@codon]&;
      fwdHomology =
          StringTake[mutatedSeq, {1, 21} + res * 3];

      revHomology =
          If[
            precedingSubstitutionFlag,
            dnaReverse@StringTake[mutatedSeq, {-26, -6} + res * 3],
            dnaReverse@StringTake[mutatedSeq, {-23, -3} + res * 3]
          ];

      fwdExtension =
          StringTake[mutatedSeq, {-3, 0} + res * 3];

      revExtension =
          If[
            precedingSubstitutionFlag,
            dnaReverse@StringTake[mutatedSeq, {-5, 0} + res * 3],
            dnaReverse@StringTake[mutatedSeq, {-3, 0} + res * 3]
          ];

      (* Results *)
      <|
        "Fwd" ->
            fwdExtension <>
                optimizeHomology[cds, fwdHomology, maxTrim, maxExtend],
        "Rev" ->
            revExtension <>
                optimizeHomology[dnaReverse@cds, revHomology, maxTrim, maxExtend]
      |>
    ];

(* Optimize primer homology by making it longer or shorter *)
optimizeHomology[
  template_String, (* PCR template sequence *)
  homology_String, (* Primer homology sequence *)
  maxTrim_Integer, (* Trim primer no more then this *)
  maxExtend_Integer (* Extend homology no more then this *)
] :=
    Catch@Module[
      {
        variants =
            homologyVariants[template, homology, maxTrim, maxExtend],
        testResults
      },

    (* Scan the test suite over the variants, throwing the first
    good primer *)

      testResults =
          With[{tests = primerTests@#},
            If[
              And @@ (Values@tests), (* All tests passed *)
              Throw[#],
              <|"Sequence" -> #, "Tests" -> tests|>
            ]
          ]& /@ variants;

      (* If no variant matches all requirements,
      then take the least bad one *)

      SortBy[
        testResults,
        {
          #["Tests", "Self-Complement"]&,
          #["Tests", "Melt"]&,
          #["Tests", "GC-Clamp"]&,
          #["Tests", "Terminal-T"]&,
          #["Tests", "GC"]&,
          dnaMeltingTemperature@#["Sequence"]&
        }
      ][[-1, "Sequence"]]
    ];

(* Construct a list of trimmed and extended homology regions *)
homologyVariants[template_, homology_, maxTrim_, maxExtend_] :=
    Module[
      {
        pos =
            StringPosition[
              template,
              homology
              , IgnoreCase -> True
            ][[1, 2]],
        extended,
        trimmed
      },
      extended =
          Table[
            homology <> StringTake[template, {pos + 1, pos + i}],
            {i, 0, maxExtend}
          ];

      trimmed =
          Table[
            StringDrop[homology, -i],
            {i, 0, maxTrim}
          ];

      DeleteDuplicates@Riffle[extended, trimmed]
    ];

(* Primer test suite *)
primerTests[homology_] :=
    <|
      "Melt" -> primerMeltTest@homology,
      "GC-Clamp" -> primerClampTest@homology,
      "Self-Complement" -> primerSelfComplementTest@homology,
      "Terminal-T" -> primerTerminalTTest@homology,
      "GC" -> primerGCTest@homology
    |>;

(* These checks return True if the sequence passes the test *)
(* Test for melting temperature *)
primerMeltTest[seq_] :=
    dnaMeltingTemperature@seq > melt;

(* Test for GC clamps *)
primerClampTest[seq_] :=
    Not@MemberQ[{"GC", "CG"}, ToUpperCase@StringTake[seq, -2]];

(* Test for self complementarity *)
primerSelfComplementTest[seq_] :=
    Nor @@ Table[
      With[{s = StringTake[seq, -i]},
        ToUpperCase@s === ToUpperCase@dnaReverse@s
      ],
      {i, 4, 17}
    ];

(* Test for 3' T *)
primerTerminalTTest[seq_] :=
    ToUpperCase@StringTake[seq, {-1}] =!= "T";

(* Test for presence of G/C in the 3' region *)
primerGCTest[seq_] :=
    MemberQ[Characters@ToUpperCase@StringTake[seq, -4], "G" | "C"];

(***********************************************************************)
End[]; (* `Private` *)

EndPackage[];