Package["CysteineScreen`"]

(* Mathematica Package *)
(* Created by Mathematica Plugin for IntelliJ IDEA *)

(* :Title: CysteineScreen *)
(* :Context: CysteineScreen` *)
(* :Author: Frank Buermann <fburmann@mrc-lmb.cam.ac.uk> *)

(***********************************************************************)
(* :Documentation: *)
(***********************************************************************)
PackageExport["makeSingleSitePrimers"]
Unprotect[makeSingleSitePrimers];
makeSingleSitePrimers::usage =
    "makeSingleSitePrimers[\"cds\", res, \"codon\"] gives PCR primers\
 for the mutagenesis of a residue using a specified substitution codon.";

PackageExport["makeMultiSitePrimers"]
Unprotect[makeMultiSitePrimers];
makeMultiSitePrimers::usage =
    "makeMultiSitePrimers[\"cds\", {res1, \"codon\"}, {res2, \"codon\"},\
 ...]\ generates a primer set for multiple site-directed mutagenesis.
makeMultiSitePrimers[\"cds\", {{res, \"codon\"}, ...}] generates a set\
 of mutagenesis primers sharing a common reverse primer.
makeMultiSitePrimers[\"cds\", {{res1, \"codon\"}, ...}, {{resn,\
 \"codon\"}, ...}, ...] generates primer sets for multiple sites with\
 compatible overhangs.";

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

(* Check if DNA is palindromic *)
palindromicQ[seq_] := ToLowerCase@seq === ToLowerCase@dnaReverse@seq;

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

      (* Count di-nucleotide instances. *)
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
SyntaxInformation[makeSingleSitePrimers] =
    {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};

makeSingleSitePrimers::palin = "Palindromic overhang `1`. Consider changing\
 the substitution codon or setting the \"UpstreamSubstitutionRules\"\
 option.";

Options[makeSingleSitePrimers] =
    SortBy[First]@{
      "AdjustUpstreamCodon" -> True,
      "UpstreamSubstitutionRules" -> precedingCodons,
      "MaxTrim" -> 3, (* Homology optimization parameter *)
      "MaxExtend" -> 6(* Homology optimization parameter *)
    };

makeSingleSitePrimers[
  cds_String, (* Coding sequence *)
  res_Integer, (* Residue number to be mutated *)
  codon_String, (* Codon for substitution *)
  op : (_Integer | Automatic) : Automatic, (* Overhang position *)
  opts : OptionsPattern[]
] :=
    Module[
      {
        precedingCodon, (* Codon in front of the mutation *)
        precedingSubstitutionFlag, (* Preceding codon substituted? *)
        mutatedCds,
        fwdHomology,
        fwdExtension,
        overhang,
        overhangPos =
            If[op === Automatic, res * 3 - 3, op],
        extensionRegion,
        revHomology,
        revExtension,
        maxTrim = OptionValue["MaxTrim"],
        maxExtend = OptionValue["MaxExtend"],
        upstreamSubstRules = OptionValue["UpstreamSubstitutionRules"]
      }
      ,
      precedingCodon = StringTake[cds, {0, 2} + 3 * res - 5];

      mutatedCds =
          If[
            IntegerQ@op
                ||
                Not@OptionValue["AdjustUpstreamCodon"]
                ||
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
              dnaTranslateCDS@precedingCodon /. upstreamSubstRules
            ]
          ](* Do substitution *)
              // dnaMutateCDS[#, res, ToUpperCase@codon]&;

      overhang =
          StringTake[mutatedCds, {0, 3} + overhangPos];

      (* Check if overhang is palindromic *)
      If[
        palindromicQ@overhang,
        Message[makeSingleSitePrimers::palin, overhang]
      ];

      (* Find the minimal region containing overhang and mutation *)
      extensionRegion =
          MinMax[
            {
              (* First base of the overhang and and upstream substitution *)
              overhangPos - 2 * Boole@precedingSubstitutionFlag,
              overhangPos + 3, (* Last base of the overhang *)
              res * 3 - 2, (* First base of the mutant codon *)
              res * 3 (* Last base of the mutant codon *)
            }
          ];

      fwdHomology =
          StringTake[mutatedCds, {1, 21} + extensionRegion[[2]]];

      revHomology =
          dnaReverse@
              StringTake[mutatedCds, {-21, -1} + extensionRegion[[1]]];

      fwdExtension =
          StringTake[mutatedCds, {overhangPos, extensionRegion[[2]]}];

      revExtension =
          dnaReverse@
              StringTake[
                mutatedCds,
                {extensionRegion[[1]], overhangPos + 3}
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

Protect[makeSingleSitePrimers];

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
(* The functions return True if the sequence passes the test *)
primerTests[homology_] :=
    <|
      "Melt" -> primerMeltTest@homology,
      "GC-Clamp" -> primerClampTest@homology,
      "Self-Complement" -> primerSelfComplementTest@homology,
      "Terminal-T" -> primerTerminalTTest@homology,
      "GC" -> primerGCTest@homology
    |>;

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
        palindromicQ@s
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
(**** Multiple site-directed mutagenesis ****)

(* End ligation fidelity table *)
fidelityTable = Get["CysteineScreen`FidelityTable_37deg`"];

(* Function for making a lookup table for ligation efficiency
and failure fraction *)
makeFidelityLookupTable[fidelityTable_] :=
    Module[
      {
        m = fidelityTable["LigationMatrix"],
        scores, failureFraction
      },
      scores = Rescale[Diagonal@m, MinMax@Diagonal@m];
      failureFraction =
          Total[m * UnitStep@DiagonalMatrix@Table[-1, Length@m]]
              / Total[m];

      (* Make <|"AAAA" -> <|"EfficiencyScore" -> 0, "FailureFraction" -> 0.4|>, ... |> association *)
      AssociationThread[
        fidelityTable["ColumnLabels"] ->
            (
              AssociationThread[
                {"EfficiencyScore", "FailureFraction"} -> #
              ]&
                  /@ N@Transpose[{scores, failureFraction}]
            )
      ]
    ];

(* Function for making a table of cross-reactivity
 (ratio between off-target reactivity and on-target reactivity) *)
makeCrossReactivityLookupTable[fidelityTable_] :=
    AssociationThread[
      fidelityTable["RowLabels"] ->
          (
            AssociationThread[
              fidelityTable["ColumnLabels"] -> Rescale@#
            ]& /@ N@fidelityTable["LigationMatrix"]
          )
    ];

(* Filter function for defining the ends of a overhang search interval *)
designFilter["Forward", minmax_, dist_] :=
    First@minmax - 4 - dist <= # < First@minmax &;

designFilter["Reverse", minmax_, dist_] :=
    Last@minmax < # <= Last@minmax + 4 + dist &;

designFilter["Mixed", minmax_, dist_] :=
    First@minmax - 4 - dist <= # <= Last@minmax + 4 + dist &;

(* Generate a sorted list of all possible overhangs *)
Options[proposeOverhangs] =
    SortBy[First]@{
      "OverhangSearchDistance" -> 2,
      "Design" -> "Forward"
    };

proposeOverhangs[
  cds_,
  sites_, (* Positions of first base of targeted codons *)
  opts : OptionsPattern[]
] :=
    Module[
      {
        cdsUpperCase = ToUpperCase@cds,
        forbidden, (* Sites that cannot be part of an overhang *)
        allowed (* Sites that can be part of an overhang *)
      },

      forbidden =
          Union[
            Flatten[Range[#, # + 2]& /@ sites]
          ];

      allowed =
          Select[
            Complement[ (* Discard sites that will be targeted by mutagenesis *)
              Range[StringLength@cds],
              forbidden
            ],
            designFilter[ (* Trim allowed region as specified by the "Design" option *)
              OptionValue["Design"],
              MinMax@forbidden,
              OptionValue["OverhangSearchDistance"]
            ]
          ];

      (* Generate a list of {site, overhang} pairs *)
      If[
        Abs@*Subtract @@ MinMax@# == 3, (* Check for disruption by forbidden regions *)
        {First@#, StringTake[cdsUpperCase, MinMax@#]},
        Nothing
      ]&
          /@ Partition[allowed, 4, 1]
    ];

(* Check if a set of overhangs is ok *)
overhangValidatorF[
  crossReactivityTable_,
  maxCrossReactivity_
][overhangs__] :=
    With[
      {
        (* If there is only one site, only check for palindromes *)
        (* If there are more sites, check cross-reactivity as well *)
        validQ =
            If[
              Length@{overhangs} <= 1,
              Not@palindromicQ@#2&, (* Only one site *)
              ( (* More than one site *)
                Not@palindromicQ@Last@#1 && Not@palindromicQ@Last@#2
                    &&
                    crossReactivityTable[[Last@#1, dnaReverse@Last@#2]]
                        <= maxCrossReactivity
              )&
            ]
      },
      (* Do the checking *)
      Apply[
        And,
        validQ @@@
            If[
              Length@{overhangs} > 1,
              Subsets[{overhangs}, {2}], (* Check all pairwise combinations *)
              {overhangs} (* Check only one site *)
            ]
      ]
    ];

(* From sets of candidate overhangs, choose a good one *)
Options[chooseOverhangSet] =
    SortBy[First]@{
      "CrossReactivityThreshold" -> 0.05,
      "MinimumEfficiency" -> 0
    };

chooseOverhangSet[
  candidateGroups_,
  fidelityTable_,
  opts : OptionsPattern[]
] :=
    Module[
      {
        fidelityLookupTable =
            makeFidelityLookupTable@fidelityTable,
        crossReactivityLookupTable =
            makeCrossReactivityLookupTable@fidelityTable,
        scoringFunction,
        selectorFunction,
        sortedOverhangs,
        found = False,
        result
      },

      (* Function for scoring by normalized ligation efficiency *)
      scoringFunction =
          With[
            {
              (* For normalization *)
              maxEfficiency =
                  Max@fidelityLookupTable[[All, "EfficiencyScore"]]
            },
            (* Normalized efficiency *)
            #EfficiencyScore / maxEfficiency &
          ];

      (* Function for selecting by minimum ligation efficiency *)
      selectorFunction =
          With[
            {
              maxEfficiency =
                  Max@fidelityLookupTable[[All, "EfficiencyScore"]]
            },
            #EfficiencyScore / maxEfficiency
                >= OptionValue["MinimumEfficiency"]&
          ];

      (* Sort overhangs and drop inefficient ones *)
      sortedOverhangs =
          Reverse@*SortBy[
            scoringFunction@fidelityLookupTable@ToUpperCase@Last@#&
          ]
              @*Select[
            selectorFunction@fidelityLookupTable@ToUpperCase@Last@#&
          ]
              /@ candidateGroups;

      (* Look at all combinations of overhangs between the groups and
      return the first good set *)
      result =
          Catch@
              Outer[
                If[
                  overhangValidatorF[
                    crossReactivityLookupTable,
                    OptionValue["CrossReactivityThreshold"]
                  ]@##,
                  found = True; Throw[List@##]
                ]&,
                Sequence @@ sortedOverhangs,
                1 (* Level *)
              ];

      If[found, result, {}]
    ];

(* Top level function for multiple site-directed mutagenesis *)
SyntaxInformation[makeMultiSitePrimers] =
    {"ArgumentsPattern" -> {_, __, OptionsPattern[]}};

makeMultiSitePrimers::noovh =
    "No suitable overhangs found. Try adjusting the\
 \"OverhangSearchDistance\" option.";

Options[makeMultiSitePrimers] =
    DeleteDuplicates@
        Join[
          {
            "OverhangSearchDistance" -> 2, (* Maximum distance between overhang and outermost mutation *)
            "Design" -> "Forward",
            "OverhangFidelityTable" -> Automatic,
            "CrossReactivityThreshold" -> 0.05,
            "MinimumEfficiency" -> 0
          },
          Options[makeSingleSitePrimers]
        ];

makeMultiSitePrimers[
  cds_String,
  mutations : {_Integer, _String}..,
  opts : OptionsPattern[]
] :=
    Flatten@makeMultiSitePrimers[cds, Sequence @@ (List /@ {mutations}), opts];

makeMultiSitePrimers[
  cds_String,
  groups : {{_Integer, _String}..}..,
  opts : OptionsPattern[]
] :=
    Module[
      {
        targetSiteGroups = {groups}[[All, All, 1]] * 3 - 2, (* First base positions of target codons *)
        overhangCandidateGroups,
        fidelity =
            If[
              OptionValue["OverhangFidelityTable"] === Automatic,
              fidelityTable,
              OptionValue["OverhangFidelityTable"]
            ],
        overhangs,
        tasks
      },

      (* Propose possible overhangs for each group *)
      overhangCandidateGroups =
          proposeOverhangs[
            cds, #,
            FilterRules[
              {opts, Options[makeMultiSitePrimers]},
              Options@proposeOverhangs
            ]
          ]& /@ targetSiteGroups;

      (* Check if there are groups with no candidates *)
      If[
        MemberQ[overhangCandidateGroups, {}],
        Message[makeMultiSitePrimers::noovh];
        Return[$Failed]
      ];

      (* Choose one overhang per group *)
      overhangs =
          chooseOverhangSet[
            overhangCandidateGroups,
            fidelity,
            FilterRules[
              {opts, Options["makePrimerSets"]},
              Options[chooseOverhangSet]
            ]
          ];

      (* Check if the search failed *)
      If[
        overhangs === {},
        Message[makeMultiSitePrimers::noovh];
        Return[$Failed]
      ];

      (* Pair overhangs with mutations *)
      tasks =
          <|"CDS" -> cds, "Overhang" -> #1, "Mutations" -> #2|>&
              @@@ Transpose[{overhangs, {groups}}];

      (* Construct primer pairs for each mutation *)
      Map[
        With[
          {seq = #CDS, op = #Overhang[[1]]},
          makeSingleSitePrimers[
            seq, #[[1]], #[[2]], op,
            FilterRules[
              {opts, Options[makeMultiSitePrimers]},
              Options[makeSingleSitePrimers]
            ]
          ]& /@ #Mutations
        ]&,
        tasks
      ]
    ];

Protect[makeMultiSitePrimers];

(***********************************************************************)
