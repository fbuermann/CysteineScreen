# CysteineScreen Package

## Tutorial 1: Primer design for high-throughput cysteine mutagenesis

High-throughput mutagenesis is a very useful tool for setting up cysteine cross-linking screens. This tutorial shows how to automatically design primers for mutating hundreds of residues in the *smc* gene of *B. subtilis* to cysteine. Most of the code resides in the accompanying **CysteineScreen** Wolfram Language package which exposes the function **makeSingleSitePrimers**.

### Package installation

Make sure that the package folder has been placed in the directory opened by



```mathematica
SystemOpen@FileNameJoin[{$UserBaseDirectory, "Applications"}]
```

Alternatively, you can download the .paclet file from the **CysteineScreen** release tab on GitHub and install it like this:



```mathematica
PacletInstall["please/enter/path/to/.paclet"]
```

Once the package is installed, this tutorial notebook can be accessed from the “Add-Ons and Packages” section in the Documentation Center or by entering “CysteineScreen” into the documentation search box.

Now load the package:



```mathematica
<< CysteineScreen`
```

The package exposes the function **makeSingleSitePrimers**:



```mathematica
? makeSingleSitePrimers
```

[//]: # (No rules defined for Output)

### Mutagenesis strategy

We will generate primers for mutating more than 400 individual coiled-coil residues of the Smc protein to cysteine. The targeted region of the *smc* gene will be split into two fragments. The first fragment is generated by PCR with a “constant” forward primer and a “variable” reverse primer. The second fragment is generated with a variable forward primer and a constant reverse primer. Matching PCR fragments are joined by Golden Gate assembly using *Bsa*I sites in the variable primers. Their *Bsa*I overhangs contain the cysteine mutation and are constrained to GTGC or TTGC to reduce variations in ligation efficiency.  For some positions, this requires recoding of the codon upstream of the mutation. The **makeSingleSitePrimers** function from the **CysteineScreen** package will take care of this. The two constant primers are used to clone the recombinant DNA into a gene targeting construct and have to be designed manually (not covered here).

### Choice of target residues

First, let’s construct a list of residue numbers that we want to mutate to cysteine. We will target residues in the Smc coiled-coil arm. The arm is an antiparallel intramolecular coiled coil, which means that it is comprised of an N-terminal and a C-terminal helix that follow the canonical heptad repeat pattern of coiled coils. We will mutagenize the residues that do not mediate a contact between the helices, i.e. which are not at the **a** and **d** positions of the heptad repeat.

Let’s define a function that takes a residue range and a starting heptad position and writes out each residue together with its heptad position:



```mathematica
heptadTable[start_, stop_, heptad_] :=
  With[
    {
      residues = Range[start, stop],
      hpos = FirstPosition[CharacterRange["a", "g"], heptad]
    }
    ,
    Riffle[
      residues,
      RotateLeft[CharacterRange["a", "g"], hpos - 1], 
      {2, 2 * Length@residues, 2}
     ]
       // Partition[#, 2]&
  ];
```

We use this function and knowledge about the coiled-coil register (see Waldman *et al.*, 2015) to assign heptad positions to Smc residues:



```mathematica
heptadPositions = 
  Join[
    heptadTable[174, 209, "a"],
    heptadTable[227, 381, "a"],
    heptadTable[395, 493, "d"],
    heptadTable[679, 777, "a"],
    heptadTable[794, 948, "d"],
    heptadTable[990, 1025, "d"]
  ];
```

Now we use pattern matching to select residues that are not at positions **a** or **d**:



```mathematica
targets =
  heptadPositions /. {{_, "a" | "d"} -> Nothing, {res_Integer, _} :> res};
```

### Automated design of the variable primers

We make use of the **makeSingleSitePrimers** function from the **CysteineScreen** package to design the variable forward and reverse primers. The function takes the sequence template, the target residue and a substitution codon (“TGC” for cysteine in our case). It returns the parts of the forward and reverse primers that contain the homology regions plus a 5’ extension which introduces the mutation. The function does not add *Bsa*I restriction sites, so we need to take care of this.

The **makeSingleSitePrimers** function adheres to the following rules:

   + If required, use an alternative codon upstream of the mutation to give a GTGC or TTGC overhang upon *Bsa*I cleavage.

   + If possible, have a homolgy region with a melting temperature of more than 50 °C.

   + If possible, avoid self-complementarity of the last 4-17 bases.

   + If possible, do not end in GC clamps.

   + If possible, do not end with T.

   + If possible, have at least one G/C within the last four bases.

First, we load the coding sequence for a cysteine free variant of the Smc protein: 



```mathematica
smc = Import[FileNameJoin[{NotebookDirectory[], "ExampleData", "BsSmc_cysless.txt"}], "Text"];
```

Now we construct primers for the target residues:



```mathematica
primers = makeSingleSitePrimers[smc, #, "TGC"]& /@ targets;
```

The primers are still missing the *Bsa*I recognition site. Let’s fix that:



```mathematica
primersFwd = StringJoin["gttacaggtctca", #]& /@ primers[[All, "Fwd"]];
primersRev = StringJoin["tcattgggtctct", #]& /@ primers[[All, "Rev"]];
```

### Primer export

Construct primer names composed of the direction and the target residue number:



```mathematica
namesFwd = StringJoin["Fwd", ToString@#]& /@ targets;
namesRev = StringJoin["Rev", ToString@#]& /@ targets;
```

Construct a table of name/sequence pairs and export as a spreadsheet:



```mathematica
Export[
  FileNameJoin[{$HomeDirectory, "CysteineScreenPrimers.xlsx"}],
  "Sheets" ->
    {
      "Fwd" -> Transpose[{namesFwd, primersFwd}],
      "Rev" -> Transpose[{namesRev, primersRev}]
    }
    , "Rules"
]
```

Open the file and have a look:



```mathematica
SystemOpen@%
```

## Tutorial 2: Multiple site-directed mutagenesis

If you want to introduce multiple point mutations at different sites, you need to make sure that restriction overhangs are compatible. The **makeMultiSitePrimers** function takes care of this. It chooses sets of overhangs based on ligation efficiency and cross-reactivity data from Potapov *et. al.*, 2019, *ACS Synth Biol*.

Install and load the package as described above.



```mathematica
<< CysteineScreen`
```



```mathematica
?makeMultiSitePrimers
```

[//]: # (No rules defined for Output)

### Mutagenesis strategy

We will generate primers for screening a set of 5 cysteines in one region of the Smc protein for cross-linking to a set of 5 cysteines in another region. The mutations will be introduced by the forward primers, and each region will have a shared reverse primer.

### Choice of target residues

We will mutate residues located at two different sites, close to R558 and N634, respectively.
Let’s prepare lists of {residue number, codon} pairs:



```mathematica
site1 = Transpose[{Range[556,560], ConstantArray["TGC", 5]}];
site2 = Transpose[{Range[632,636], ConstantArray["TGC", 5]}];
```

### Automated primer design

As before, we load the coding sequence for a cysteine free variant of the Smc protein: 



```mathematica
smc = Import[FileNameJoin[{NotebookDirectory[], "ExampleData", "BsSmc_cysless.txt"}], "Text"];
```

Let’s call the **makeMultiSitePrimers** function on the residue lists. This will make sure that overhangs for the two sites can be used together in an assembly reaction. Forward primers will contain the mutations, reverse primers will be identical for each site (this can be changed by giving the “Design” -> “Reverse” option).



```mathematica
primers = makeMultiSitePrimers[smc, site1, site2];
```

The primers are still missing the *Bsa*I recognition site. Let’s fix that:



```mathematica
primersFwd = Map[StringJoin["gttacaggtctca", #]&, primers[[All, All, "Fwd"]], {2}] // Flatten;
primersRev = Map[StringJoin["tcattgggtctct", #]&, primers[[All, All, "Rev"]], {2}] // Flatten;
```

### Primer export

Construct primer names composed of the direction and the target residue number:



```mathematica
targets = Flatten[{site1, site2}\[LeftDoubleBracket]All, All, 1\[RightDoubleBracket]];
namesFwd = StringJoin["Fwd", ToString@#]& /@ targets;
namesRev = StringJoin["Rev", ToString@#]& /@ targets;
```

Construct a table of name/sequence pairs and export as a spreadsheet (this also removes duplicate reverse primers):



```mathematica
Export[
  FileNameJoin[{$HomeDirectory, "MultiSitePrimers.xlsx"}],
  "Sheets" ->
    {
      "Fwd" -> Transpose[{namesFwd, primersFwd}],
      "Rev" -> DeleteDuplicatesBy[Transpose[{namesRev, primersRev}], Last]
    }
    , "Rules"
]
```

Open the file and have a look:



```mathematica
SystemOpen@%
```

