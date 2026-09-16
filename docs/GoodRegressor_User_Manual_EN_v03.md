# GoodRegressor
## User manual and worked output guides

**GoodParser | GoodDesigner | GoodCurator | GoodRegressor**

User manual edition 0.3 - 16 September 2026

This manual describes GoodParser, GoodDesigner, GoodCurator and GoodRegressor, together with the worked examples `Example001_BCCDuctile` and `Example002_SocialPsychology`. The algorithm is introduced in Seong-Hoon Jang, *GoodRegressor: A Hierarchical Inductive Bias for Navigating High-Dimensional Compositional Space*, [doi:10.48550/arXiv.2510.18325](https://doi.org/10.48550/arXiv.2510.18325), including its Supplementary Information (SI).

The emphasis is practical: what each input controls, how files connect, what the example actually produced, and which conclusions the outputs do and do not support. Values labelled **Example** below are values in the supplied files, not universal defaults or recommended settings for every dataset.

**Essential terminology.** In `GoodRegressor.cpp` and its model reports, the legacy label `test` means **validation**: these data participate in model selection. It does not mean an independent held-out test set. Literal filenames, input keys and output headers are retained in this manual so that they remain searchable. The preprint's separately reserved outer benchmark test folds retain their genuine test-set meaning.

**Scope.** Chapter 1 explains the input files and parameters. Chapters 2 and 3 explain the saved outputs of the materials and social-psychology examples. The appendix provides an output-file dictionary, a source map and an explanation of the acceleration mechanisms. Example results are reported from the supplied calculations.

# Reader's map

Chapter 1 is the input reference. Read 1.1 before editing any configuration, then use the section for the module being run. Section 1.6 explains the small model-chain files, stacking coefficients and other auxiliary inputs that are easy to overlook.

Chapter 2 follows Example001 from raw formulae to the Yield and Fracture regression ensembles. It explains record counts, stage selection, fitted coefficients, model performance, interaction summaries and symbolic exports.

Chapter 3 follows Example002 from its numeric social-psychology records through cohort-specific model construction, interpretation, genuine outer-fold evaluation and the separate machine-learning benchmark. It explains why legacy validation results, full-cohort fit and coverage-filtered outer-test metrics must remain distinct.

Appendix A is an output-file dictionary. Appendix B maps the main operations to their source locations. Appendix C explains Cholesky updates, on-the-fly learning and active MPI, and is linked from the `TargetCPUTimeBeforeSwap` parameter. `CodeOnly/Utility` contains the selected-model and benchmarking scripts; `input_examples` accompanies this manual with the complete configuration fragments introduced below.

## Reading sources and examples

Source references use the following short keys. Line numbers refer to the extracted source snapshot; they will move when the program is edited. Function names are provided where helpful.

| Key | Source |
| --- | --- |
| P | [GoodRegressor preprint, doi:10.48550/arXiv.2510.18325](https://doi.org/10.48550/arXiv.2510.18325); main-text page numbers and printed SI page numbers |
| GP | `CodeOnly/GoodParser.cpp` |
| GD | `CodeOnly/GoodDesigner.cpp` |
| GC | `CodeOnly/GoodCurator.cpp` |
| GR | `CodeOnly/GoodRegressor.cpp` |
| E | A path relative to `Example001_BCCDuctile/` |
| E2 | A path relative to `Example002_SocialPsychology/` |
| BM | `CodeOnly/Utility/benchmark_v14.py` and `benchmark_pp_v14.py`; the saved example calculations also include copies under E2, `005_MLBM`. |
| ED | [Wolfram Language, ElementData](https://reference.wolfram.com/language/ref/ElementData.html), used in preparing the general atomic-property tables. |
| SP | Related social-psychology preprint, [doi:10.31234/osf.io/waqjn_v2](https://doi.org/10.31234/osf.io/waqjn_v2). |

The manual's workflow is **Parser -> Designer (descriptors) -> Curator -> Regressor -> Designer (analysis/prediction/design)**. A user with an already prepared numeric TSV can enter directly at Curator or Regressor. The parser and materials-specific descriptor machinery are not prerequisites for a general numerical regression problem. [P, SI S3-S9; GR, input reader.]

# 1. Input parameter reference

## 1.1 Shared conventions, compilation and execution

### Configuration files are ordered records

Each executable reads its own configuration from the current working directory:

| Executable | Required configuration filename | Parallel interface in this snapshot |
| --- | --- | --- |
| `GoodParser.x` | `basic_properties_GoodParser.txt` | Serial C++ |
| `GoodDesigner.x` | `basic_properties_GoodDesigner.txt` | MPI C++ |
| `GoodCurator.x` | `basic_properties_GoodCurator.txt` | Serial C++ |
| `GoodRegressor.x` | `basic_properties_GoodRegressor.txt` | MPI C++ |

The long text before the value on a configuration line is generally a descriptive label consumed by an ordered reader. It is **not** a free-form key/value configuration language. Keep the original order and required lines. Do not insert a comment or blank line among active fixed-position records. Whitespace separates tokens; TSV data tables use tabs. Sentinels such as `STOP`, `Functypes` and `BetaRegression` have control-flow significance. [GP 217-265; GC 103-234; GR 1290-1563; GD 3430-4050.]

Use a period as the decimal separator, whitespace-free file paths and numeric labels, and a final newline. Use LF line endings for Linux execution. Keep each fixed-position configuration record on its own physical line.

**One working directory per independent job.** The sources create temporary files and contain wildcard cleanup operations. Examples include Curator's `*.temporary.txt`, Regressor's `*.reg`, `*.pr.txt` and `smt`, and Designer's design-worker files. Never run multiple jobs in one directory or use a directory containing unrelated files matching these patterns. Work from a copy, not the only copy of archived results. [GC/GR/GD, cleanup branches.]

### Build commands and dependencies

The supplied `compile_example.txt` gives these commands. Run them inside `CodeOnly`, preserving the bundled Eigen directory structure:

```sh
g++ GoodParser.cpp -o GoodParser.x -std=c++11
mpicxx GoodDesigner.cpp -o GoodDesigner.x -std=c++11
g++ GoodCurator.cpp -o GoodCurator.x -std=c++11
mpicxx GoodRegressor.cpp -o GoodRegressor.x -std=c++11
```

Eigen is the external linear-algebra library identified by the package. The MPI modules additionally require a working MPI compiler wrapper and runtime.

Executables read configuration/data relative to their working directory. For example, in a prepared, isolated Parser folder:

```sh
/path/to/CodeOnly/GoodParser.x
```

In an allocated MPI environment, a generic launch has the following shape; replace the process count and executable path with values appropriate to the allocation:

```sh
mpirun -np 8 /path/to/CodeOnly/GoodRegressor.x
```

This is a launch pattern, not a scheduler submission script. The preprint's large-core benchmarks are not minimum hardware requirements and do not establish the runtime of Example001 on a workstation. A configured 1000-second sampling stage is also not a 1000-second wall-time guarantee for the entire pipeline. [P, pp. 24-27; E, Regressor configurations.]

### Keep numerical meaning and display names separate

`AveXensity` denotes the average gravimetric density (mass per unit volume), as distinct from the `MolDensity` descriptor. Retain the literal spelling `AveXensity` in numeric headers and `setX`. A LaTeX display label does not rename a numeric column. Record physical units separately: changing units changes the arguments of nonlinear transforms and is not automatically compensated by the program.

Treat zero, missing data and the finite `VInfinite` sentinel as different concepts. The normal Designer statistics omit zero-valued atomic attributes; an absent tabular value must not be silently converted to a meaningful zero. Check all transformed inputs and predictions for finiteness. Domain restrictions of powers, logarithms and divisions still apply. [GD, `get_mamssk`, 2019-2093; GR, `get_v`, 900-940.]

## 1.2 GoodParser: formulae and aligned target values

**Role.** Convert chemical-composition strings into an explicit block representation and append accompanying target values. Parser prepares the material records consumed by Designer; it does not calculate the full descriptor matrix. [P, SI S6; GP 217-265, 328-454, 616-673.]

### Fixed input records

| Position / parameter | Example | Meaning and operational consequence |
| --- | --- | --- |
| 1. `InputFilename` | `BCCYield.txt` | Raw composition file. Each data record must be parseable by this program's grammar. |
| 2. `PostscriptFilename` | `BCCYield_postscript.txt` | Text/target values appended to parsed compositions by **row position**, not by a join on identifiers. |
| 3. `LogFilename` | `BCCYield_log.txt` | Parsing diagnostics. Read this before accepting the final record count. |
| 4. `OutputFilename` | `BCCYield_output.txt` | Intermediate formula-parsing output and associated composition diagnostics. |
| 5. `FinalOutputFilename` | `BCCYield_dat.txt` | Final material records, including the appended postscript/target data; the handoff to Designer. |
| Following block list | `None`, then `Og`, then `STOP` | Assign species to composition blocks; these blocks carry through to subsequent material processing. |

The working BCC configuration is:

```text
InputFilename BCCYield.txt
PostscriptFilename BCCYield_postscript.txt
LogFilename BCCYield_log.txt
OutputFilename BCCYield_output.txt
FinalOutputFilename BCCYield_dat.txt
======AtomBlocks(IfNone,"None")======
None
Og
STOP
```

### Atom blocks and formula grammar

`None` provides the fallback block for species not assigned to a named block. The example places ordinary alloy species in the first block and `Og` in its own block. A block can contain a list of species. The reader also supports a `*` no-count marker for block-level aggregate diagnostics; that option is not used for `Og` in the archived example. It must not be confused with the multiplication symbol in Designer's user-defined descriptors. [GP, block reader and final-output writer.]

The source contains handling for numeric decimal stoichiometry, nested parentheses/brackets, composition mixtures and selected molar/atomic/weight-percentage notations. This is not a general symbolic chemistry engine: resolve algebraic variables such as `x` or `1-x` before submitting a formula, and verify any less common syntax in the intermediate output. Never infer successful parsing solely from an executable exit code. [GP 328-454.]

Examples of mixture notation accepted by the parsing workflow include:

```text
Mg85Al15 + 5 wt% Ni
0.2(CeO2) + 0.8(Fe2O3)
```

The first combines an Mg-Al composition with a weight-percentage Ni addition. The second uses numerical mixture coefficients for two parenthesized oxide compositions. Keep the stated weight or mixture notation in the input; examine the explicit species-content pairs in the Parser output to check the resulting interpretation.

The final representation is divided by literal `||` separators. It contains identifying/formula text, explicit species-content pairs grouped by the configured blocks, aggregate diagnostic fields and the appended postscript. Preserve the separators and block order. Downstream code does not expect this representation to be a simple two-column formula/target CSV.

### What `Og` means in Example001

The raw string `Al0.214Nb0.714Ta0.571Ti1V0.143Zr0.929Og298.15` encodes both alloy content and a condition value. Designer's descriptor declaration is `Og U Og * 1 L= T M= T`, so the `Og` amount is extracted as the feature displayed as `T`. It is a temperature carrier in this example, **not a proposed oganesson-containing alloy**. The numeric headers do not explicitly declare the units of `T`, `Yield` or `Fracture`; this manual does not infer units from their magnitudes. [E, 001_Parser/BCCYield raw input; 002_Designer/BCCYield configuration, line 44.]

Because the archived Parser block is `Og`, not a no-count `Og *` block, aggregate mass/content diagnostics include the pseudo-species. Do not read those totals as the physical mass or atomic fraction of the alloy. Descriptor construction later handles the carrier differently; neither representation licenses chemical substitution of that carrier during materials design.

### Checks before continuing

Compare raw-formula, postscript, successful-output and downstream-input record counts. Inspect at least one simple alloy, one fractional composition and one mixed/block case. Check that each target remains associated with the intended composition after any failures or reordering. The appending procedure is positional; manually deleting a failed formula but not its matching postscript can shift all later targets.

## 1.3 GoodDesigner: descriptors, ensemble analysis and design controls

**Role.** Designer has two substantially different uses. With no model files, it calculates descriptors from materials, atomic attributes and optional structural data. With model chains, it evaluates symbolic models, constructs or loads stacking ensembles, reports interaction statistics, generates feature-scan tables and can optimize compositions. The same configuration has records for both uses, so unused design values still occupy their prescribed positions. [P, SI S6, S8-S9; GD 3430-4050.]

### Files and data interpretation: records 1-6

| Parameter | Value type / example | Explanation |
| --- | --- | --- |
| `AtomDatListFilename` | Path; `BCC_general_atomdat.txt` | Atomic-property lookup table, prepared using Wolfram `ElementData` [ED]. Descriptor names must match its column names. **Always retain the `Valence` column**, including when it is not selected as a regression descriptor, because nominal-valence bookkeeping uses it. Keep the corresponding table with the trained models. |
| `FaceFilename` | Path or `-` | Optional face-dependent data; the labelled format requires `AtomName` and `Face` as its first two columns. Unused in Example001. |
| `StructureListFilename / UseStrName` | `BCC_strdat.txt 0` | Structural descriptor table and structure-name selection. The BCC example uses the `None` structure with `0`. For multiple named structures, the supplied reader uses `1`; see the format example below. |
| `MaterialsFilename / HeadExists / PostTagNames` | `BCCYield_dat_shuffled.txt 0 Yield` | Parsed materials, whether a header is present, and names of trailing reference/target columns. These names also associate observations with model groups. |
| `RegModelFilename` | List of chain files, `-`, or list ending in `!` | `-`: descriptors without models. A list: load model chains. A standalone `!` token on the list disables further composition design, allowing analysis/screening. |
| `MapOutSpaceFilename` | Path or `-` | Map-space input is parsed, but the associated main processing call is commented out in this source snapshot. Do not advertise this as a verified turnkey mode. |

Designer uses MPI in the supplied code, even though the preprint's high-level workflow particularly emphasizes MPI for the Regressor. Use the build command for the actual source, not an inferred serial Designer build.

`PostTagNames` are literal identifiers. The group name `Yield` in a model-chain file must match the `Yield` post-tag when observations are needed for fitting/evaluation or offset correction. A file prefix such as `Y005001.txt` does not itself define the target group. Grouping is based on metadata inside the chain file. [GD, `read_reg_chain` and material/model matching.]

### Structural data: `StructureListFilename`, `S` and `B`

The materials record carries its structure label separately from the parsed composition:

```text
(DataEntryName) || (StrName) || (Parsed composition) || (Target metrics)
```

The parsed-composition field may itself contain `||` separators for its site/group blocks. Match the `StrName` in each materials record to the corresponding row of `StructureListFilename` when using multiple named structures. In the supplied reader, `UseStrName = 1` enables that lookup; `0` uses the `None` structure. [GD, `read_strdat`; main setup near 9458-9462.]

The structure table begins with the descriptor list and the ordered site list. Its header then gives **descriptor-major, site-minor** columns: all sites for the first descriptor, then all sites for the next descriptor, and so on. Terminate the table with `STOP`.

```text
DescriptorNames Coor OmegaA OmegaS PV
SiteNames A B C
```

For the BCC/FCC example, the complete header order is:

```text
StrName Coor_A Coor_B Coor_C OmegaA_A OmegaA_B OmegaA_C OmegaS_A OmegaS_B OmegaS_C PV_A PV_B PV_C
```

The following table displays the same two structure rows vertically so that every value remains readable. The ready-to-use row-oriented file is `input_examples/structure_BCC_FCC.txt`.

| Column | BCC | FCC |
| --- | --- | --- |
| `Coor_A` | 14 | 12 |
| `Coor_B` | X | X |
| `Coor_C` | X | X |
| `OmegaA_A` | 0.866646249 | 0.897597901 |
| `OmegaA_B` | X | X |
| `OmegaA_C` | X | X |
| `OmegaS_A` | 0.248872795 | 0.399887002 |
| `OmegaS_B` | X | X |
| `OmegaS_C` | X | X |
| `PV_A` | 53.05512343 | 54.90571602 |
| `PV_B` | X | X |
| `PV_C` | X | X |

`X` means that the specified site does not exist in that structure; it is not a numerical zero. `Coor` is the coordination number. `OmegaA` is the mean site-to-face solid angle of the polyhedron formed by surrounding atoms, `OmegaS` is the standard deviation of those solid angles, and `PV` is the polyhedron volume. Keep the units and geometrical definitions used to prepare the supplied values consistent throughout the input.

Declare these structural quantities in the `AnalyzeDescriptors` block as follows:

```text
Coor S L= n_{c} M= Subscript["n", "c"]
OmegaA S L= \Omega_{m} M= Subscript["\[Omega]", "m"]
OmegaS S L= \Omega_{\sigma} M= Subscript["\[Omega]", "\[Sigma]"]
PV S L= V_{p} M= Subscript["V", "p"]
```

`S` means **StructureOnly**: the value is selected from the structural descriptor/site table.

For a descriptor that depends on both the species and its structural environment, use **`B` (Both)**. The following coordination table illustrates the site values used to select coordination-dependent ionic radii:

```text
DescriptorNames Coor
SiteNames A B C
```

| StrName | Coor_A | Coor_B | Coor_C |
| --- | --- | --- | --- |
| Perovskites | 12 | 06 | X |
| Perovskite-related,brownmillerites | 07 | 05 | X |
| Perovskite-related,Ruddlesden-Popper | 09 | 06 | X |
| Perovskite-related,Dion-Jacobson | 08 | 05 | 06 |
| Perovskite-related,hexagonal | 12 | 06 | X |
| Perovskite-related,others(V-layered) | 05 | 06 | X |
| Perovskite-related,others(Ti-layered) | 06 | 06 | X |
| Perovskite-related,others(Nd-layered) | 06 | 07 | 06 |
| Perovskite-related,others(A-site-vacancies) | 08 | 06 | X |
| Fluorites | 08 | X | X |
| d-Bi2O3 | 08 | X | X |
| BIMEVOX | 05 | 06 | X |
| Scheelites | 08 | 04 | X |
| Apatites | 08 | 04 | X |
| Cuspidines | 07 | 04 | X |
| Melilites | 08 | 04 | X |
| Garnets | 08 | 05 | X |
| Pyrochlores | 08 | 06 | X |
| Columbites | 06 | 06 | X |
| a-SrSiO3 | 08 | 04 | X |
| Bi2UO6 | 05 | 10 | X |

The corresponding row-oriented input, including its header and final `STOP`, is `input_examples/structure_coordination.txt`. These are the supplied descriptor assignments for the example, rather than a general crystallographic definition of every structure in a class.

Then add the following declaration to `AnalyzeDescriptors`:

```text
Radius B Coor
```

Designer uses the site's `Coor` value to select a property column from `AtomDatListFilename`: for example, `04` selects `Radius04`, `05` selects `Radius05`, and `06` selects `Radius06`. It then obtains that radius for the corresponding species. This is the mechanism for using coordination-dependent Shannon radii without assigning one radius to every environment. Preserve the zero-padded coordination strings and provide the matching `Radius04`, `Radius05`, etc. columns in the atomic table.

### Objective, reliability and scale: records 7-11

**`TargetVReg`** supplies the desired value for each model group, in the target's working units. It is a target-location objective, not an implicit instruction to maximize every property. For the revised single-target Fracture setup, the desired value is `10000`.

**`TargetStdDvReg`** provides a maximum permitted inter-model prediction dispersion for each group. Its units follow the corresponding modeled target, including any logarithmic transformation of that target. This is an ensemble-consistency screen, not a calibrated probability, confidence interval or error guarantee.

**`ZDistWeight`** supplies the per-target weights in the standardized distance to the desired target. Higher weight gives that target more influence in the distance. Weights belong to target groups, not to individual models in the ensemble. Use a nonnegative, deliberately chosen set for a distance objective.

**`UseAve`** and **`UseStdDv`** supply optional target-centering and scaling values; `-` leaves a value to be determined by the program. A user-supplied standard deviation must be a meaningful positive scale. Do not derive an evaluation scale from all-zero placeholders and then interpret it as a measured-data standardization.

Removing the standalone `!` from the model list enables further materials design. For that use, provide `TargetVReg`, `TargetStdDvReg`, `ZDistWeight`, `UseAve` and `UseStdDv` with one entry for each target metric, in the same model-group order. For two target metrics, each record therefore needs two entries. Further multi-target composition design with this setup is not demonstrated in these examples.

### Composition search: records 12-20

**`ModulateContentFactor`** is a list controlling allowed amount changes and doping proposals. Positive values are content multipliers, such as `0.9 1.0 1.1`. Nonpositive entries are read through the doping-factor branch using their absolute values; the archived post-processing example uses `1.0 -0.1`, corresponding to no positive scaling change plus a 0.1 doping fraction. Do not interpret `-0.1` as a negative atom count. Avoid adding unexplained zero or negative factors merely to populate the list. [GD, parameter reader and material modulation branches.]

**`UseOffset`** is `0` for direct model predictions or `1` for an observed-value-anchored change. With `1`, matching reference post-tags must be present. The latter prediction convention is:

```text
y_modified = y_observed_original
             + f(modified composition) - f(original composition)
```

This is distinct from `f(modified composition)`. State which convention is used when reporting candidates. It does not validate the extrapolated difference. [P, SI S42, Eq.S65; GD, offset branches.]

**`MaxIteration`** limits the composition-improvement iteration loop; Example001 uses `100`. It is not the number of Regressor bagging runs, MPI ranks or depth stages.

**`NeighborSelectSpace`** lists atomic-attribute names used for replacement-neighbor selection, such as `Electronegativity MetalRadius` or the same list with `MValEleN`. These are atomic attributes, not already aggregated `Ave...` descriptor names.

**`NeighborZCrit`** has two operational forms. A positive first value selects a scalar standardized-distance criterion. A nonpositive first value selects componentwise bounds, using the magnitudes of the listed entries in the corresponding attribute dimensions. The example `-0.3 -0.5` therefore means separate difference limits; it is not a negative Euclidean radius. Ensure the number and ordering match the selected neighbor attributes. [GD, neighbor-selection setup.]

**`ExcludeNeighbors`** removes specified species from the proposed replacement-neighbor list. It accepts explicit species or the predefined groups below. These are the program's literal lists, not general definitions of toxicity, rarity or electronic configuration.

| Input group | Species excluded |
| --- | --- |
| `Toxic` / `toxic` | Hg Pb Cd As Tl Be Cr Cs Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og |
| `Rare` / `rare` | Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu |
| `Non-quenched` / `NQ` | Mn Fe Co Ni Cu Ru Rh Pd Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Os Ir Pt Au |

Use `NQ` as the no-hyphen spelling for the last group in the supplied reader. A hyphen disables the `ExcludeNeighbors` record before that reader reaches its group-name matching. Explicit species lists can be used instead. [GD, exclusion-list reader, 3616-3715.]

**`DoNotFindNeighbors`** distinguishes protected species from partial-doping-only species. Unprefixed names are not substituted with other species. A name prefixed by `-` can remain unchanged or be partially substituted, but is not fully replaced. The migrating-ion example below illustrates both forms. Also inspect any `<species>.nb.txt` files supplying explicit neighbor lists.

**`AnionName`** selects the composition/charge bookkeeping. A single anion name uses that anion together with `AnionValence`; `-` disables that single-anion mechanism. The extended `*` form instead supplies species/nominal-valence pairs whose contents can adjust for charge compensation, as shown below.

**`AnionValence`** gives the charge in single-anion bookkeeping. Its numerical value is not used for the `AnionName * ...` example below, which supplies the valences directly. Retain the fixed-position record in the configuration even when its value is unused.

### Example: protect migrating ions while permitting charge compensation

```text
DoNotFindNeighbors(IfNone,"-") Li Na Mg Ag K Ca Zn -H -N -F -S -Cl -Se -Br -I -O
AnionName(IfNone,"-"//YouCanAddCationNamesAlso...Try) * Li 1 Na 1 Mg 2 Ag 1 K 1 Ca 2 Zn 2 limit 0.1
```

Keep each displayed record on one physical line in the input file; line wrapping on this page is for reading only. The same records are supplied in `input_examples/migrating_ion_controls.txt`.

In `DoNotFindNeighbors`, `Li Na Mg Ag K Ca Zn` are protected from substitution. The prefixed entries `-H -N -F -S -Cl -Se -Br -I -O` allow only partial substitution (doping) or no change.

The leading `*` in `AnionName` activates the species/nominal-valence list. When other atoms are doped or replaced with different nominal valences, the amounts of the listed species can be modulated for charge compensation. Their identities remain protected even though their contents can change. In this branch, `limit 0.1` bounds the candidate content relative to its starting content within factors 0.9 and 1.1. `AnionValence` is unnecessary for this calculation because the valences are given in the list. **Always leave `Valence` in the atomic-property table** for the remaining nominal-valence bookkeeping. [GD, extended anion reader and content check near 7465.]

The oxide constraints described in the preprint are settings for those demonstrations, not universal guarantees of phase stability or experimental feasibility. Neither neighbor similarity nor small ensemble dispersion replaces structural or experimental validation. [P, SI S42-S43.]

### Outputs, statistics and interaction analysis: records 21-27

| Parameter | Example | Meaning and caution |
| --- | --- | --- |
| `VInfinite` | `100000000000000000` | Finite sentinel used in exceptional/undefined numerical situations. Do not treat it as a physical measurement or as an accuracy setting. |
| `Verbose` | `0` | Diagnostic verbosity; `1` requests more detail. |
| `OutputFilename` | `BCCYield_designeroutput.txt` | Main report, including model/stacking information and design progress when enabled. |
| `AnalysisOutputFilename` | `BCCYield_anaoutput.txt` | Descriptor/analysis table. An optional integer on the same record supplies the analysis mesh. |
| `DoNotAnalyze` | `-` | Species excluded from descriptor statistics. This is not a list of descriptor columns to remove. |
| `AnalyzeInteractionLevel` | `2 3` | Maximum interaction-set level and optional number of chains. Use **`1 1`** when detailed descriptor/interaction analysis is not needed: that analysis can become a calculation bottleneck. |
| `SkipMinMaxKurtoInAnalyzeInteraction` | `1` | Exclude Min/Max/Kurto from interaction analysis. `0` retains them. A distinct `2` branch is mean-only analysis and must not be treated as an ordinary Boolean. It retains zero-valued attributes in `get_mamssk`; Example002 uses it for singleton numeric observations. |

The `1` setting does not remove Min/Max/Kurto columns from the descriptor table. Example001's initial analysis still contains all six statistics for each selected atomic attribute; Curator later chooses the intended modeling subset. [GD, analysis/output branches; E, 002_Designer tables.]

### Descriptor declaration block

After the fixed records, preserve the descriptor-block heading, add declaration rows, then terminate with `STOP`.

```text
Mass A L= M M= M
ShearModulus A L= G M= G
Og U Og * 1 L= T M= T
STOP
```

`A` denotes an atomic attribute. `S` denotes a structural attribute. `B` denotes a value depending on both atomic and structural information and requires the corresponding postfix/mapping information. These modes must agree with the atom/structure tables. A specialized `P` branch reads an index and converts to atomic-style handling internally; it is not demonstrated here and is not documented as a general portable descriptor recipe. [GD, descriptor reader and post-processing of descriptor definitions.]

`U` defines a quantity from species contents. The numerator and denominator/product side may select species, use `All`, exclude species using `!`, or contain a fixed number. Exactly one `*` or `/` operator separates the two sides. For example, `n U Na * 1` retains sodium content, while a selected-species sum divided by `All` defines a content ratio. The `Og` declaration extracts the stored condition without averaging it with alloy atom attributes. `L=` and `M=` are display/export strings for LaTeX and Mathematica; they do not change numerical evaluation.

### More `U` declarations: contents, fractions and ratios

```text
n U Na * 1
rTR U Ni Sc Ti V Cr Mn Fe Co Cu Zn Y Zr Nb Mo Tc Ru Rh Pd Ag Cd Hf Ta W Re Os Ir Pt Au Hg Rf Db Sg Bh Hs / All
rTRe U Ni Sc Ti V Cr Mn Fe Co Cu Zn Y Zr Nb Mo Tc Ru Rh Pd Ag Cd Hf Ta W Re Os Ir Pt Au Hg Rf Db Sg Bh Hs / ! Ni Sc Ti V Cr Mn Fe Co Cu Zn Y Zr Nb Mo Tc Ru Rh Pd Ag Cd Hf Ta W Re Os Ir Pt Au Hg Rf Db Sg Bh Hs
```

`n` is the Na content multiplied by 1. `rTR` is the sum of the explicitly listed species contents divided by the total content (`All`). `rTRe` uses the same numerator, but the denominator excludes that same species list using `!`, so it is the ratio of the listed content to the remaining content.

For **either side** of a `U` declaration, choose an explicit species list, `All`, an exclusion list beginning with `!`, or a fixed number. Use exactly one operator, either `*` or `/`, between the sides. A zero denominator is not a usable ratio. The long declarations above remain single-line records in `input_examples/AnalyzeDescriptors_U.txt`; page wrapping is only typographical.

### Atomic statistics: what the numbers actually mean

For an attribute `a` and retained species contents `c_j`, normal Designer processing uses stoichiometric weighting. In compact notation:

```text
w_j = c_j / sum(c_j over species retained for this attribute)
Ave(a)   = sum(w_j * a_j)
StdDv(a) = sqrt(sum(w_j * (a_j - Ave(a))^2))
Skew(a)  = sum(w_j * (a_j - Ave(a))^3) / StdDv(a)^3
Kurto(a) = sum(w_j * (a_j - Ave(a))^4) / StdDv(a)^4
```

The normal branch **excludes zero-valued atomic attributes** from the retained set, as well as excluded species. Thus the normalization may differ between properties. A zero can behave as missing/unavailable atomic information rather than as a true zero contribution. The separate mean-only branch changes this behavior. For zero spread the implementation handles skewness and kurtosis specially; a sentinel kurtosis is not a physical heavy-tailed distribution. [GD, `get_mamssk`, 2019-2093.]

Example001 selects 15 attributes: `Electronegativity`, `Mass`, `Valence`, `MValEleN`, `MetalRadius`, `Xensity`, `MolDensity`, `MFilling`, `BulkModulus`, `ShearModulus`, `PoissonRatio`, `ThermalConductivity`, `ThermalExpansion`, `DebyeTemp` and `MMSusc`. Normal processing yields 90 atomic-statistic columns, plus the user-defined `Og` feature. `Valence` and `MValEleN` are distinct supplied attributes; do not merge them merely because both relate to electrons.

## 1.4 GoodCurator: filtering and first-level interactions

**Role.** Read a numerical feature table, select target/feature spaces, optionally handle duplicated names or composition proximity, and construct an initial pool of products/ratios. Its output is a report containing a tabular block, not a pure ready-to-read TSV file from first line to last. [P, SI S7; GC 103-234 and dataset/output routines.]

### Parameter reference

| Position / parameter | Example | Explanation |
| --- | --- | --- |
| 1. `DataFilenamePrefix` | `BCCYield_anaoutput` | Input prefix; the program appends `.txt`. Supplying the extension twice points at the wrong path. |
| 2. `LogOutputFilename` | `BCCYield.txt` | Main Curator report, including the curated `DATASET` block. Despite the label, it is not merely a short log. |
| 3. `OutputPrecision` | `20` | Decimal output precision, not source-data or floating-point accuracy. |
| 4. `NameDuplicateHandling` | `-1` | `-1`: leave duplicates; `0`: retain minimum target; `1`: maximum target; `2`: modal/histogram-based handling with supplied minimum, maximum and increment. |
| 5. `XSpace` | `Ave StdDv Skew Og` | Header-substring selectors for the base features. Matching is not an exact-name-only list. |
| 6. `InteractXSpace` | `Ave Og` | Header-substring selectors for features eligible for first-level products/ratios; `-` disables this additional space. |
| 7. `ParseData` | `-` | Optional parsed-material file for composition-based handling; disabled in Example001. |
| 8. `ParseGroup` | `2` | Number of groups inside the parsed-composition field: one plus its internal `||` count. Not the number of targets or folds. |
| 9. `setY` | `Yield * *` | Target-column name and lower/upper bounds; `*` means no corresponding bound. |
| 10. `SubFilterX` | `-` | Optional feature-name/lower-bound/upper-bound triples applied as record filters. |
| Following sampling block | Heading, then `STOP` | Optional postfix/action/factor rows for composition-proximity selection flags. |

The Curator needs unambiguous identifiers. Its name-column detection uses a `Name` substring. Feeding a Designer table that still contains both `Name` and `CompositionName` is unsafe; prepare one intended `Name` column and numeric columns. Example001 removes the three other metadata columns before Curator. [GC, input header selection; E, 002/003 tables.]

`NameDuplicateHandling` operates on repeated data names, not on independent measurements automatically grouped by chemical similarity. Confirm that the identity policy suits repeated experiments at different conditions. A maximum-target policy is a scientific curation choice, not a generic remedy for every duplicate.

In this source, a hyphen can disable the entire `SubFilterX` line. Therefore a negative numerical bound is not automatically safe to enter under the displayed generic syntax. Check the reader before using such a bound. A working-copy filtering step outside the program may be clearer than relying on this ambiguous branch; record that preprocessing explicitly. [GC 103-234.]

### Products and ratios

Curator iterates over selected feature pairs in their existing order. The demonstrated pool contains a product and a forward ratio for a pair, not both possible ratio directions. Ratio eligibility requires fixed nonzero sign over the data for both participating variables in the implemented screening. A column with zeros or a sign change can therefore remove an expected ratio even when it was named in `InteractXSpace`. [P, SI S7; GC, interaction construction.]

In Example001, `XSpace` chooses 46 base features: 15 means, 15 standard deviations, 15 skewnesses and `Og`. `InteractXSpace` chooses 16 of these: 15 means plus `Og`. All pairs represented in this example produce the demonstrated product/ratio pool:

```text
46 base features + 2 * C(16,2) interactions = 286 predictors
286 predictors + Name + target = 288 table columns
```

Those are measured properties of this example's feature pool, not an unconditional feature-count formula for a new dataset.

### Optional composition-proximity sampling

`ParseData`, `ParseGroup` and the sampling block let Curator identify a preferred sample among compositionally close records, using either the minimum or maximum target value.

```text
ParseData(*.txt,IfNone,"-") materials_dat.txt
ParseGroup 3
```

Here `materials_dat.txt` is the parsed `MaterialsFilename` used by GoodDesigner. Leave the intervening fixed records in their original positions. In the sampling block, use:

```text
Postfix_SamplingAction(TowardsMin0;TowardsMax1)_SamplingFactor
U0.01 1 0.01
STOP
```

A materials record has the following structure:

```text
(DataEntryName) || (StrName or user-defined name) || (Parsed composition) || (Target metrics)
```

**`ParseGroup` counts the groups inside the parsed composition.** It equals one plus the number of `||` separators within that composition field, excluding the separators for the leading metadata and trailing targets. Thus three parsed composition groups require `ParseGroup = 3`.

Curator normalizes the composition contents so that their sum is 1. It compares records in that normalized-composition space. When their distance is below `SamplingFactor`, the action chooses which target value to favor: **`0` (TowardsMin) retains the lower value; `1` (TowardsMax) retains the higher value.**

In the example, `U0.01 1 0.01` favors the maximum target for compositions closer than 0.01 in that normalized space. The output named by `LogOutputFilename` gains a column `U0.01`: `1` marks a retained record and `0` marks a record not retained by this selection. Additional sampling rows can create additional postfix columns. The output flags can then be applied as data filters, for example `U0.01 1 1` in Regressor's `SubFilterX`. Producing the flags and filtering on them are separate operations. [GC, `pick_data` and sampling-indicator output.]

### Handoff to Regressor

Extract the rectangular block beginning after the `DATASET` marker, including its header and only its material rows. Exclude the report preamble and trailing statistics. Retain literal feature names, including parentheses and operators. Save the result as `input.txt`, then set Regressor's `setX` to the chosen predictor names and `setY` to the target.

The companion `tools/extract_curator_dataset.py` performs this extraction conservatively and checks that each selected row has the expected number of tab-separated columns. It refuses to overwrite an existing output. It does not train models, change descriptors or repair malformed numeric entries.

## 1.5 GoodRegressor: search, depth and statistical selection

**Role.** Search explicit symbolic forms under a controlled hierarchy. The preprint describes run-through, swap, transit, pick and repeated ensemble construction. The code searches a sampled lexicographic space and refines selected models; a reported best model is not proof of exhaustive enumeration or a globally optimal equation. [P, pp. 24-27; GR, search/refinement routines.]

### Data and output records

**`DataFilename`** identifies the numeric TSV input; the reader can accept a list of compatible input files. All included tables must have consistent predictor/target headers and numerical meaning. `Name` is the optional observation identifier. Avoid whitespace inside feature names and preserve exact expression-style headers.

**`OutputFilename`** is the report prefix, `output` in the example. The generated filenames additionally identify the target, depth stage and requested term-count mode. `output_Yield_2L4I.txt` is not a separate model-target called `2L4I`; it is the second stage of the Yield search with the indicated four-term limit and interaction mode.

**`OutputPrecision`**, here `20`, controls printed precision. Twenty printed digits do not imply twenty reliable digits in a measured target, an optimized coefficient or a cross-platform rerun.

**`tFTestMesh`**, here `1000`, controls the numerical mesh used by the implemented t/F significance calculations. It is not the number of folds, permutations or candidate models. Changing it changes the numerical significance calculation rather than the train/validation split. [GR 1290-1563 and statistical routines.]

### `NVariableLimit`: the depth schedule

The example uses:

```text
NVariableLimit(IfNegative,Limit) -5 -4I -3I -2I -1I
```

The magnitude is the number of active additive terms at the stage. A negative sign requests a strict limit, rather than permitting the corresponding additional-variable behavior; it does not request negative variables. The suffix selects the pool used at subsequent stages:

| Token form | Meaning |
| --- | --- |
| `-5` | Strict five-term stage using the original selected input pool. |
| `-4i` | Strict four-term stage using the preceding stage's interaction pool. |
| `-4I` | Strict four-term stage using that interaction pool together with the original input pool. |

The first stage has no preceding-stage interaction pool to reuse. Preserve capitalization: `i` and `I` are not synonyms. Output `L` marks the limited mode; the suffix identifies the interaction-pool mode. [GR, schedule parsing and inter-stage dataset construction.]

As the additive term count decreases, retained transformed expressions and their products/ratios can make individual terms more deeply composed. Thus fewer additive terms do **not** necessarily mean a simpler expression tree. The preprint's interaction-depth axis refers to this construction, not merely a generic feature-count sweep. [P, pp. 25-26; SI S8.]

The schedule supplies stages, not independent ensemble repeats. The example's 10 model members live in 10 run directories outside this one configuration. Do not describe `-5 -4I -3I -2I -1I` as five-fold cross-validation or five-member bagging.

### Significance, time and local-refinement controls

**`PickUpPLevel`** is the significance threshold, `0.05` in the example. The full-Fisher/full-frequentist criterion described in the preprint requires the overall F-test and coefficient/intercept t-tests to pass the selected threshold. This is a model-admission rule inside a search; it is not a causal test of the resulting descriptors or an independent post-selection discovery probability. Look for the report's `full_Frequentist` flag rather than assuming every provisional printed model satisfies it.

**`Verbose`** selects log detail (`0` or `1`). It does not change the conceptual distinction between training, validation and external testing.

**`TargetCPUTimeBeforeSwap`** controls the initial sampling budget and the refinement options. For most cases, start with one of the following two presets:

```text
TargetCPUTimeBeforeSwap -100 -1 !!
TargetCPUTimeBeforeSwap -1000 -1 !!
```

These are alternatives: keep **one** `TargetCPUTimeBeforeSwap` record in a configuration. Use `-100 -1 !!` for a shorter initial exploration, or `-1000 -1 !!` for a larger initial sampling budget. Both use the accelerated swap option and on-the-fly selection. The three fields have the following meanings:

| Field | Meaning |
| --- | --- |
| `-100` or `-1000` | Magnitude: run-through time target in seconds. Negative sign: enable the `useSuperFastSwap` branch. |
| `-1` | No positive time cap is imposed on the after-swap refinement by this field. |
| `!!` | Enable `useSomeShotMode`, the on-the-fly selection mode; do not enable `useExitMode`. |

For the punctuation alternatives, `!` enables only `useExitMode`, `!!` only `useSomeShotMode`, and `!!!` both. These flags are Regressor options, not shell syntax, and have a different meaning from Designer's standalone `!`.

[Appendix C: Calculation-speed features](#appendix-c-calculation-speed-features) explains Cholesky factorization updates, on-the-fly learning and the LUHMANN-to-HEGEL active-MPI workflow.

The run-through budget is not a total job timeout: later swaps, transformations, depth stages, I/O and MPI coordination add work. More CPU resources or altered time settings can change which models are visited. Use the two presets above as starting choices and retain the full configuration with the results. [GR, time-option reader; Appendix C.]

**`R2Threshold`**, `0.05`, is used in the additional-variable improvement logic, including a training-R2 improvement check and associated validation conditions. It is not a universal minimum accepted R2, and it should not be presented as the convergence threshold for every stage of swap/transit. Strict negative term limits may prevent the addition that this setting controls. [GR, approximately 4707-4728.]

**`MaxBFGSIterations`**, `10`, limits the relevant beta-regression numerical optimization. It does not limit the lexicographic search to 10 models or generate 10 bagged members. Baseline OLS is the demonstrated Example001 mode.

### Splitting, duplicates and robustness

**`TrainRatio`**, `0.8`, controls the training fraction. The remaining rows form the validation subset, although report fields call them `test`. For Example001 the resulting sizes are 215/54 for Yield and 122/31 for Fracture. Split membership can change during the staged computation; inspect identifiers rather than inferring membership from the ratio alone.

**`DuplicateX`**, `1`, is the multiplicity of predictor entries in the search machinery; `1` keeps one copy, whereas values above one create repeated entries available to the search. It is not a request to duplicate observations and is not the bagging count. Repeated or highly collinear transformed predictors require interpretation of the actual fitted design matrix and diagnostics.

**`NameDuplicateHandling`** uses the same displayed policy family as Curator: `-1` retains repeated names; `0`/`1` retain smaller/larger targets; `2` uses the provided histogram settings. Check name identity and condition information before discarding repeated names. Example001 uses `-1`.

**`NTestRobustness`**, `100`, controls the number of robustness resamplings. The current routine selects subsets from the combined training and validation predictions and evaluates the already fitted model. It does **not** fit 100 new models. It is neither 100-fold independent validation nor an external test set. [GR 7041-7147; P, SI S8.]

The `ROBUSTNESS` summary reports the minimum, mean, maximum and standard deviation of each metric over the resamples. The corrected implementation takes the square root of the corresponding mean squared deviations for MSE, RMSE, MAE and R2. These describe metric variability across resampled subsets, not the inter-model prediction dispersion reported by Designer.

### Target and predictor selection

**`setY`** supplies target names with lower/upper filtering bounds; `*` leaves a bound unrestricted. Example: `Yield * *`. Bounds filter records; they are not the optimizer's desired target. Designer's `TargetVReg` is a different input.

**`setX`** explicitly lists the predictor headers used by the Regressor. Unlike Curator's substring selectors, preserve the actual column/expression names here. The example lists all 286 curated predictors. A feature table can have extra unused columns, but numerical columns cannot be renamed without updating this list and any downstream model interpretation.

**`SubFilterX`** supplies additional feature-name/lower-bound/upper-bound triples, or `-` for no filter. For example:

```text
SubFilterX(IfNone,"-";{Name/Min/Max:IfNoLimit,Min/Max=*}*NSubFilterX) AveMass 0 50 U0.01 1 1
```

This retains only data entries satisfying **`0 <= AveMass <= 50` and `U0.01 = 1`**. The latter applies Curator's composition-proximity selection flag. The names must match columns in the supplied numeric table.

### `Functypes`: scalar transformations

The fixed records are followed by a `Functypes` heading and transformation rows. For example:

```text
Functypes
1 ^1 pow a= 0 b= 1
2 ^-1 pow a= 0 b= -1
3 log10 log a= 0 b= 10 c= 1
BetaRegression(IfNo,"none") none
```

This is a minimal **illustration**, not the Example001 catalogue. Each row contains an index/display label, a function family and ordered coefficients. Keep row order stable for archived model evaluation. The printed display label is not an independently interpreted mathematical expression: actual evaluation follows the function family and coefficients. [GR, transform reader; GD, reading original Regressor configuration.]

| Family | Implemented mathematical form |
| --- | --- |
| `pow` | `(x + a)^b` |
| `log` | `[ln(x + a) / ln(b)]^c` |
| `sin`, `cos`, `tan` | `[func((b * (x + a))^c)]^d` |
| `sinh`, `cosh`, `tanh` | `[func((b * (x + a))^c)]^d` |
| `exp` | `[d^(c * (x + a)^b)]^e` |
| `erf`, `erfc` | `[func(b * (x + a))]^c` |
| `abs` | `abs((x + a)^b)` |

The last form computes the power before the absolute value; it is not necessarily numerically equivalent to `abs(x+a)^b` for fractional powers on negative inputs. Some explanatory comments in the configuration differ from the implementation; the table above follows `get_v`. The `erfc` branch calls the complementary error function. [GR 900-940.]

**The `BetaRegression` line terminates this read.** Example001 has 67 active scalar-transform rows before that boundary at configuration line 86. Additional trigonometric transform listings later in the file are inactive text, not part of that run's search. The preprint's 109-transform setting must not be copied into the description of these archived runs. Read the rows before that boundary to identify the active families.

The example uses numerical constants such as `0.333333` and `2.71828`; they are the supplied numerical values, not exact symbolic `1/3` or the exact mathematical constant e. Use the original coefficients when reproducing predictions.

### `BetaRegression`

`none` leaves the demonstrated baseline regression mode. The source supports beta-link options including `logit`, `probit`, `cloglog`, `cauchit` and `nloglog` in the beta setting. They belong on the actual controlling record as required by the reader, not on a later standalone explanatory line. Beta regression has additional numerical and bounded-response considerations and is not exercised by these BCC archives. Do not describe it as interchangeable with the validated OLS route without a dedicated check. [P, SI S47; GR, beta setup and optimization.]

## 1.6 Auxiliary inputs that control post-processing

### Selected-model files for GoodDesigner

Files such as `005001.txt`, `Y005001.txt` and `F005001.txt` are small manifests, not datasets or single flat equations. They identify the target, beta mode, original Regressor configuration and chronological stage reports needed to evaluate the selected model.

These selected-model files are generated to be supplied to GoodDesigner as a post-process. Use `CodeOnly/Utility/find_maxR2.py` to read the search summaries and create the files for the recommended models. The worked examples show the resulting manifests and their use in `RegModelFilename`.

The archived manifests contain absolute paths beginning `/work/k0736/k073600/GoodRegressor/`. Before relocation, rebase both the configuration path and every listed stage-report path in a copy. Ensure each path exists, the target group name is preserved, and the original transform catalogue is still the one used to train the model. Blindly moving only the small manifest is insufficient.

Optional text after the main chain boundary supports additional post-processing conventions, but no such advanced configuration is needed for the demonstrated selected models. Preserve the known-good syntax instead of inventing an extra formula block.

### Utility scripts

The `CodeOnly/Utility` directory contains the companion scripts:

| Script | Purpose |
| --- | --- |
| `find_maxR2.py` | Read the Regressor recommendations and generate selected-model files for Designer post-processing. |
| `benchmark_v14.py` | Run the conventional machine-learning and symbolic-baseline benchmark workflow. |
| `benchmark_pp_v14.py` | Post-process the benchmark logs and prediction files into summary tables. |

Example002 demonstrates the benchmark inputs and outputs. Review the script's path and dataset settings for the intended working directory before execution.

### Stacking coefficient files: `.se.txt`

Designer forms a consensus prediction of the form:

```text
f_stack(x) = m0 + sum_i(m_i * f_i(x))
```

It is a least-squares stacking fit, not necessarily an arithmetic average. Coefficients can be negative. The progressive rows in an ensemble coefficient file correspond to accumulation of model members in their loaded order. Do not reorder the chain list while keeping the same coefficient file. [P, pp. 8-9; SI S8-S9; GD, ensemble construction.]

The source checks for target-named coefficient files such as `Yield.se.txt` and `Fracture.se.txt`. When these are present, they supply the archived stacking coefficients. Without them, the relevant path can fit coefficients to the current reference post-tags. For an unlabeled candidate grid with zero placeholders this would not be a meaningful training operation: carry the coefficients from the labeled fitting dataset instead. The dual-target folder includes these saved files.

The current Designer dispersion is calculated from individual-model predictions using absolute stacking weights around the unweighted member mean:

```text
mean_member = sum_i(f_i) / number_of_members
ModelStdDv = sqrt(sum_i(abs(m_i) * (f_i - mean_member)^2)
                  / sum_i(abs(m_i)))
```

This describes the supplied numerical routine, not a generic confidence interval. It is not the residual RMSE of the consensus predictor and it does not include every source of predictive uncertainty. [GD, `calc_bag_chain`, 3209-3233.]

### Optional `fold.txt`

The optional fold input can restrict which material rows are analyzed. After its heading the source reads `NFolds`, `TakeFold`, training/test-side selection and straight/rotating assignment choices, with an optional `TestModelStDevLimit`. `TakeFold` uses a one-based fold number. In straight assignment, membership is determined from row position in contiguous blocks; in rotating assignment it is row index modulo the fold count. [GD, optional fold reader near 9378 and material-selection logic.]

A useful workflow is to fit the stack with **`Train0Test1 = 0`**, then evaluate it with **`Train0Test1 = 1`**. The first run prints its stacking intercept and coefficients in the final block of the training analysis output. Copy that block into a target-named coefficient file, such as `satisfaction.se.txt`, in the test directory. GoodDesigner reads this file automatically and uses the fixed coefficients for the test rows instead of fitting a new stack to them. Preserve the model order. Example002 shows this operation explicitly in Section 3.10.

This file does not automatically train all outer-fold Regressor models or prevent data leakage. A genuine outer-test evaluation requires the fitted descriptors/choices, model search and stacking fit to exclude those test outcomes. Record the exact folds and model manifests used. Any dispersion-based filtering of test predictions must report its retained coverage, not only the accuracy of the accepted subset.

### Explicit neighbor files and structural metadata

A `<species>.nb.txt` file can override the automatically assembled replacement list. Search the working directory for these inputs when a prohibited or unexpected replacement appears. Configuration-only reasoning is incomplete when auxiliary files change the neighbor set.

`BCC_strdat.txt` in this example provides a minimal `None` structure and zero-valued coordination entries with named sites; it is not an atomistic BCC structure file and does not contain relaxed atomic coordinates. The alloy class comes from the example's data selection, not from a crystallographic reconstruction performed by the pipeline.

# 2. Example001: BCC Yield and Fracture

## 2.1 What is in the archive

`Example001_BCCDuctile` contains two labeled modeling tasks and a later candidate-screening use of their two ensembles. It is not the oxygen-ion-conductor benchmark shown in the preprint. Use the preprint for the algorithm and concepts, but use this archive for its actual inputs, model choices and results.

```text
Example001_BCCDuctile/
  001_Parser/
    BCCYield/                         raw formulae -> parsed materials
    BCCFracture/
  002_Designer/
    BCCYield/                         parsed materials -> descriptors
    BCCFracture/
  003_Curator/
    BCCYield/                         selected descriptors -> interactions
    BCCFracture/
  004_Regressor/
    BCCYield/005001 ... 005010/        ten staged searches
    BCCFracture/005001 ... 005010/
  005_Designer/
    BCCYield/                         Yield ensemble, interpretation, design
    BCCFracture/                      Fracture ensemble, interpretation, design
    BCCYieldFracture_DualDesign/       two-ensemble candidate screening
```

| Quantity | Yield | Fracture |
| --- | ---: | ---: |
| Labeled material records | 269 | 153 |
| First Designer descriptor-table columns | 96 | 96 |
| Prepared table columns supplied to Curator | 93 | 93 |
| Curated predictor columns | 286 | 286 |
| Regressor input columns including Name/target | 288 | 288 |
| Training / validation records per stage | 215 / 54 | 122 / 31 |
| Archived run directories | 10 | 10 |
| Depth stages per run | 5 | 5 |

These counts were read from the supplied files and their model reports. The table headers lack explicit physical units for the targets; all RMSE and MAE values below are in the **stored target units**. The manual deliberately avoids guessing that a raw number must be MPa, percent or kelvin.

## 2.2 Stage 001: read the parsed result before calculating features

The first Yield raw record and its target illustrate the separation of chemical composition, condition and response:

```text
Raw:     Al0.214Nb0.714Ta0.571Ti1V0.143Zr0.929Og298.15
Target:  1965
```

The corresponding parsed record separates the alloy species from the `Og` block and appends the target. Read the explicit atom-content pairs, not only the display name, because later output names may use shortened decimal strings. `Og298.15` remains a condition carrier rather than a physical species in the alloy interpretation.

`BCCYield_output.txt` is useful for checking parsing and block assignment. `BCCYield_dat.txt` is the material-record handoff. `BCCYield_log.txt` should be consulted for parsing failures; an empty failure section is not a substitute for checking row correspondence. The equivalent Fracture files contain 153 material records.

## 2.3 Stage 002: understand the descriptor table

The materials files used by Designer are named `BCCYield_dat_shuffled.txt` and `BCCFracture_dat_shuffled.txt`. Their ordering differs from the preceding unshuffled materials files. Preserve each record's identifier and response when shuffling; never shuffle formula and target columns separately.

The first Designer configuration has `RegModelFilename ... -`, so the primary purpose is descriptor generation. Its analysis table has four metadata columns, 90 atomic-statistic columns, the condition feature `Og`, and one target:

```text
No.  Name  CompositionName  Strcture
     [six statistics for each of 15 atomic attributes]
     Og  Yield
```

`Strcture` is the literal header spelling. The six statistics are mean, standard deviation, skewness, minimum, maximum and kurtosis. Their literal column names include the prefixes `Ave`, `StdDv`, `Skew`, `Min`, `Max` and `Kurto`. They should not all be assumed equally suitable for regression simply because they were generated.

### The handoff is prepared, not automatic

The tables in `003_Curator` are not byte-identical copies of the entire 96-column Designer table. The archived preparation removes `No.`, `CompositionName` and `Strcture`, leaving `Name` and 92 other columns: all 90 atomic statistics, `Og` and the target. It also reflects numerical formatting/rounding in the intermediate table. The feature selection to 46 base variables occurs later inside Curator.

This distinction matters when rebuilding the example: feeding all 96 columns directly to Curator can trigger ambiguous `Name` matching, while stripping every Min/Max/Kurto column at this handoff would be a different, additional preprocessing decision. A faithful reconstruction separates metadata removal, numerical formatting and Curator's own feature selection. [E, 002_Designer and 003_Curator inputs; GC, header matching.]

## 2.4 Stage 003: read the Curator report and export the dataset

The example selects `Ave`, `StdDv`, `Skew` and `Og` for the base pool and `Ave` plus `Og` for interactions. Consequently the regression input contains both an individual feature such as `AveMass` and expression-style columns such as `(AveMass*Og)` or `(AveMass/Og)`.

Open `003_Curator/BCCYield/BCCYield.txt`. The preamble reports processing information. The `DATASET` section contains the header and 269 data rows. The report then resumes statistical summaries. The full report has 335 lines; it is not a 334-observation TSV. The Fracture report has 219 lines, of which 153 are data rows in the dataset section.

The companion extractor can be run as follows from a location where the paths resolve:

```sh
python tools/extract_curator_dataset.py \
  Example001_BCCDuctile/003_Curator/BCCYield/BCCYield.txt \
  work_yield/input.txt --expected-rows 269 --expected-columns 288
```

The output directory must already exist; the tool does not overwrite a file. Check the reported row/column count before using the result. The archived Regressor `input.txt` remains the reference when comparing a historic run: a freshly extracted table can preserve more printed precision than another intermediate preparation.

Curator's ratio columns have directional names. `(AveMass/Og)` and `(Og/AveMass)` are different numeric predictors; the example does not include every reversed ratio automatically. This becomes important when tracing a fitted term back to its underlying attributes.

## 2.5 Stage 004: a run is a five-stage search, not five independent tests

Each target has directories `005001` through `005010`. Each contains an `input.txt`, a `basic_properties_GoodRegressor.txt`, and reports for the five requested stages. Within a target, the archived input/configuration files are identical across the ten run directories; the model results differ.

The nominal progression is five terms, then four, three, two and one, with the interaction pool expanded between stages according to the `I` schedule. It does not guarantee monotonically improving predictions. It also does not guarantee identical held-out memberships between stages.

### Yield run 005001: why stage 2 was retained

| Stage / report suffix | Training R2 | Validation R2 | Robustness mean R2 | Validation RMSE |
| --- | ---: | ---: | ---: | ---: |
| `1L5` | 0.76801 | 0.90224 | 0.79523 | 161.0996 |
| `2L4I` - selected | 0.84341 | 0.78071 | 0.83295 | 195.7809 |
| `3L3I` | 0.82845 | 0.79416 | 0.83162 | 176.7308 |
| `4L2I` | 0.77439 | 0.81967 | 0.78557 | 217.2629 |
| `5L1I` | 0.72558 | 0.68000 | 0.72554 | 274.5299 |

The small manifest retains stage 2 and its stage 1 dependency. Stage1 has the largest single validation R2 in this run, but stage 2 has the largest reported robustness mean R2. Therefore selecting the file with the largest `test` R2 would **not** reproduce the archived selected chain.

Do not read this table as a fixed-test-set learning curve. The stage reports can use different training/validation memberships, and the robustness statistic reuses fixed predictions over resampled subsets of the combined data. Its comparison is the archive's internal model-selection diagnostic, not an unbiased external accuracy estimate. [E, 004_Regressor/BCCYield/005001 reports and manifest; GR, robustness routine.]

### Fracture run 005001

| Stage / report suffix | Training R2 | Validation R2 | Robustness mean R2 | Validation RMSE |
| --- | ---: | ---: | ---: | ---: |
| `1L5` | 0.62201 | 0.66364 | 0.60838 | 9.4639 |
| `2L4I` - selected | 0.66305 | 0.43746 | 0.62938 | 10.6561 |
| `3L3I` | 0.55512 | 0.59469 | 0.54616 | 10.1247 |
| `4L2I` | 0.51509 | 0.47700 | 0.51081 | 9.8969 |
| `5L1I` | 0.35822 | 0.30333 | 0.35055 | 13.0969 |

Here, too, the archived selected stage is not the one with the highest single validation R2. The Fracture task remains more weakly reproduced by these individual selected models than the Yield task in the corresponding stored target scales; this observation does not by itself identify its physical cause.

### The ten selected chains

| Run | Yield selected stage / terms | Fracture selected stage / terms |
| --- | --- | --- |
| 005001 | 2 / 4 | 2 / 4 |
| 005002 | 3 / 3 | 2 / 4 |
| 005003 | 2 / 4 | 1 / 5 |
| 005004 | 3 / 3 | 1 / 5 |
| 005005 | 2 / 4 | 1 / 5 |
| 005006 | 2 / 4 | 2 / 4 |
| 005007 | 2 / 4 | 2 / 4 |
| 005008 | 3 / 3 | 2 / 4 |
| 005009 | 2 / 4 | 1 / 5 |
| 005010 | 3 / 3 | 2 / 4 |

The selection is read from the supplied manifests. Do not assume that all ensemble members have the same term count or the same number of predecessor reports.

### Split reproducibility: what was actually checked

For each target, the first-stage validation identifiers are identical across all ten archived runs. The later stages inspected have different memberships across the ten runs. The source uses its random-selection mechanism but exposes no seed field in the documented configuration. These facts support neither the statement that all splits are identical nor the statement that ten fully independent first-stage partitions were established.

For future reproducible ensembles, record random initialization, row order, process count, time/search settings and the actual split identifiers. The preprint's intended repeated train/validation procedure must be distinguished from the split evidence available in this particular archive. A time-limited parallel search can yield different models even when its initial data partition is the same.

## 2.6 How to read a Regressor model report

Use `004_Regressor/BCCYield/005001/output_Yield_2L4I.txt` as the reference example. The report includes an initial input/data section, search progress, a final answer, fitted-term statistics, training predictions, the legacy `TEST` block, robustness summaries and the next-stage interaction information.

### `FINAL_ANSWER` and training statistics

In this report the `FINAL_ANSWER` marker is at line 807. The ensuing training summary gives four model degrees of freedom, 210 error degrees of freedom, training R2 `0.8434054967531893`, RMSE `200.911693128815` and MAE `158.398188270382`. The flag `full_Frequentist` is `1`.

`ModelDF` describes the fitted model terms in this output; `ErrorDF` is the residual degree-of-freedom count used for inference. `SStot`, `SSres`, `MSE`, `RMSE`, `MAE`, `det_r2`, `det_radj2`, `Fvalue` and `p_Fvalue` are not interchangeable metrics. A printed p-value of zero should be read as a numerical/printing result at extreme significance, not as a mathematical proof of exactly zero probability.

### Coefficient table: separate prediction from interpretation

The fitted model is an intercept plus coefficients times symbolic terms. The report includes the raw coefficient, a standardized coefficient, statistical uncertainty/significance information and collinearity diagnostics. Use the raw coefficient and exact transformed predictor to reproduce predictions; standardized coefficients are for comparison within their specified normalization, not drop-in replacements in the raw-unit equation.

The selected Yield example has these raw coefficients, rounded here for reading:

| Term | Raw coefficient | Standardized coefficient |
| --- | ---: | ---: |
| 1 | 9175.32518184 | 0.87856535 |
| 2 | -2209.48771472 | -0.33665076 |
| 3 | -27221.24202102 | -0.33404909 |
| 4 | -2526.45783891 | -0.11786138 |
| Intercept | 6365.29441034 | Not a feature-effect ranking |

For example, term 3 is `(AveMValEleN/AveThermalConductivity)^2` in the actual fitted report. The sign of its coefficient applies to that **whole composite term**, with the remaining model terms held fixed. It is not by itself an unconditional causal effect of valence-electron number or thermal conductivity.

The leading term also combines temperature, thermal expansion, density and thermal conductivity through nested transformations. Expanding it by intuition from an abbreviated label is unsafe; preserve the stage-chain definitions and original `Functypes` parameters. Shared base variables can appear in several terms with opposing contributions.

### Training and validation prediction rows

`dataset.size()` gives 215 training rows. The `TEST` marker at line 1045 introduces 54 validation rows. The validation summary gives:

```text
R2    = 0.7807093932299739
RMSE  = 195.780927504701
MAE   = 164.667892647044
```

A validation record such as `204_MoNbTaTiVWOg1273` has stored target `752.8666667` and prediction about `546.0820898`. This illustrates why a model's overall R2 should not be substituted for a per-material reliability check. Keep the observation name, actual target and predicted value together when plotting residuals or parity diagrams.

**Terminology to use in captions:** "Validation predictions; the program's legacy output label is TEST." Do not label this panel "external test". The parameters and structure were selected using validation performance during the search. [P, pp. 24-26; user terminology clarification; GR selection logic.]

### `ROBUSTNESS`: sample size, mean and standard deviation

The example's block uses 100 resamples of size 54 from the combined fitted prediction pool. Its mean R2 is `0.832947975252808`. The standard deviation corresponding to this example's resampled R2 values is approximately `0.03565`.

| Quantity | Value / meaning |
| --- | --- |
| Number of resamples | 100 |
| Rows per resample | 54 |
| Mean R2 | 0.83294798 |
| `StDev` of R2 | Approximately 0.03565 |

This statistic samples already available prediction errors, including training-data predictions. It provides the internal robustness/selection summary. It is not a substitute for a completely excluded test fold or a newly acquired dataset. [GR, robustness routine; P, SI S8.]

### `INTERACTION` and predecessor reports

The interaction section supplies expressions that can become candidate terms in the next stage. Deleting it, trimming the report to a printed coefficient table, or replacing the original transform list can break reconstruction of a later chain. Archive reports and configurations as model dependencies, not merely human-readable logs.

## 2.7 Stage 005: build and inspect the stacking ensemble

The single-target post-processing folders load the ten selected manifests. Designer groups models by target and constructs the progressive stacking ensemble, or uses matching stored coefficients when supplied. Check the coefficient source before interpreting the output.

### Output layers

`BCCYield_designeroutput.txt` contains model/ensemble information and design progress. `BCCYield_anaoutput.txt` starts with material descriptors and ends with model-performance summaries. `__R__.BCCYield_anaoutput.txt` is the prediction-oriented table. Other files describe coefficients, equations, feature interactions and scans. A long analysis file may therefore contain more than one logical table.

The prediction table's row ordering follows the material data used in that run. It does not include an independent, trustworthy chemical identity simply because its predicted values look plausible. Reattach the `Name` and parsed composition through verified row correspondence or an explicitly constructed identifier table. Do not align a shuffled prediction table to an unshuffled material list.

### Actual ensemble results

| Final whole-data ensemble | Number of labeled rows | R2 | RMSE | MAE | Invalid count |
| --- | ---: | ---: | ---: | ---: | ---: |
| `bag0_Yield` | 269 | 0.862234 | 184.433034 | 144.752732 | 0 |
| `bag0_Fracture` | 153 | 0.692446 | 8.845106 | 7.130397 | 0 |

These are final **whole-data** metrics in the archived Designer tables. They include data used in individual-model construction and ensemble fitting. They must not be reported as outer-test accuracies. They are also not the oxygen-ion-conductor benchmark numbers from the preprint.

The summary headers include `R2`, `PearsonR` and `PearsonR2`. R2 based on residual sum of squares is not generally identical to squared correlation for arbitrary predictions. Use the metric actually named in the table, especially when comparing a raw single model, a recalibrated one-member stack and a multi-member ensemble.

### Why the consensus is not a simple average

For Yield, the final intercept is about `-25.3696283`. Several stacking weights are negative, including approximately `-0.2667594`, `-0.0986623` and `-0.0927490`; the remaining weights are positive. The archive is using a fitted linear combination with an intercept, not a probability-weighted vote.

Progressive accumulation indices start at zero, so accumulation level 0 means one member and level 9 means ten. A one-member stack can already differ from the corresponding raw single model because its intercept/slope are refitted. Increasing the number of members should be judged from the full accumulation table, not inferred from the last member alone. [E, single-target Designer stacking sections and analysis summaries.]

For new candidates, freeze this ensemble before evaluation. Re-fitting a stack on candidate placeholders is not deployment of the trained model.

## 2.8 Interaction summaries and feature-scan tables

### `.zc.txt`: occurrence and coefficient-based statistics

The file `bag0_Yield.zc.txt` contains columns such as `TermName`, `InteractionLevel`, `Min(|Z|)`, `Ave(|Z|)`, `Max(|Z|)`, `InHowManyModelsAppear`, `HowManyAppear` and `WhereToFind(z@address)`.

`InteractionLevel` is the number of base descriptors in the queried set. `InHowManyModelsAppear` counts models containing that set, whereas `HowManyAppear` counts term occurrences. A descriptor may occur more than once in a model. Addresses such as `005001.txt-1-2` identify the manifest and internal chain/term location; the numerical suffixes are code-style indices, not the one-based ordinal headings in a manuscript table. [P, SI S30-S36; GD, term-analysis routines.]

In the archived Yield ensemble, `Og`, `AveThermalExpansion` and the pair `AveThermalExpansion-Og` appear across all ten models. In Fracture, `Og` appears in all ten, while `AveShearModulus` and its interaction with `Og` have lower occurrence counts. This indicates recurring representation in these fitted equations; it does not establish an independent causal intervention on temperature or modulus.

**Do not silently interpret `Ave(|Z|)` as nonnegative.** The current aggregation uses signed stacking coefficients. The archived Yield row for `AveXensity` contains `Ave(|Z|) = -0.3241703455`, despite the absolute-value notation in the header. Source inspection also shows occurrence-dependent aggregation. Consequently the field is an implementation-specific signed statistic and requires reconciliation with the intended manuscript definition before being described as a conventional absolute importance measure. [GD 6006-6037; E, bag0_Yield.zc.txt.]

The file can contain additional chain-related sections after the initial descriptor rows. A script that simply concatenates every matching name can double-count the same logical interaction. Identify the intended section and its header first.

### Feature-scan / partial-dependence-style outputs

Files such as `_bag0_Yield.AveThermalExpansion-Og.txt` tabulate a model surface while other features are fixed at their reference means. They are the implementation's feature-scan/partial-dependence-style outputs. This is not necessarily an average over the empirical joint distribution of all other features, and it is not an experimentally realizable independent variation of alloy composition. [P, SI S36-S38; GD, feature-map routines.]

Typical columns include `vAveThermalExpansion`, `vOg`, the actual `AveThermalExpansion` and `Og` coordinates, the predicted target, `ModelStdDv`, sample counts and local summary statistics. The `v...` columns are normalized mesh coordinates; use the actual feature columns for a raw-unit axis. A `-` in an empty local-data summary is not a measured zero.

For interpretation, inspect both the predicted surface and the nearby-data counts. A visually smooth low-error surface can extend into feature combinations unsupported by the measured dataset. Independent descriptor variation may also break the correlations induced by a real chemical composition.

`bag0_Yield.onevtest.txt` and the corresponding Fracture file provide additional one-variable-analysis information. They are not an extra independent held-out test dataset, despite `test` in the filename. Read their internal headers and connect them to the feature-analysis operation rather than treating them as another benchmark.

### LaTeX and Mathematica equation exports

Designer emits `.latex.txt` and `.mathematica.txt` files for further symbolic analysis. Each file records the model components and coefficients, followed by the intercept. Use the raw coefficients with the corresponding component expressions to reconstruct the numerical model; standardized coefficients support interpretation rather than substitution into the raw-unit prediction equation.

Use `.latex.txt` for mathematical typesetting and `.mathematica.txt` for Mathematica expressions. The `L=` and `M=` descriptor labels determine their respective display names. Keep the variable definitions and original numerical transform constants with the exported model. [GD, equation-export writers.]

## 2.9 Three performance statements that must remain distinct

| Statement | Supported interpretation | Unsupported interpretation |
| --- | --- | --- |
| Regressor `TEST` metrics | Validation performance used inside model search | Fully independent external-test performance |
| Designer single-target bag metrics | Whole-data fit/evaluation of the archived ensemble | The preprint's outer-fold benchmark |
| Dual-table predictions | Screening of a prepared candidate list | Accuracy against 43,758 measured mechanical-property labels |

The preprint separately describes five-fold outer benchmarking: 20% is excluded for testing, while the other 80% is split 8:2 for training and validation, giving 64/16/20 overall. It also explicitly describes full-data ensemble models constructed without a separate test set. Preserve that distinction rather than replacing every occurrence of the word `test` in the manuscript. [P, p. 9 and p. 12.]

# 3. Example002: social-psychology satisfaction models

## 3.1 What this example demonstrates

For the related social-psychology study, see [doi:10.31234/osf.io/waqjn_v2](https://doi.org/10.31234/osf.io/waqjn_v2) [SP].

`Example002_SocialPsychology` demonstrates GoodRegressor on a numeric, non-materials dataset. The response is `satisfaction`. The folders `Jpn020` and `Kor020` analyze the rows coded `Nation = 1` and `Nation = 2`, respectively. This chapter retains these archive labels. It does not infer a questionnaire, sampling frame, scale validation, or population-wide conclusion from the filenames.

The archive contains **421 observations: 201 in Jpn and 220 in Kor**. Each observation has an identifier, 20 numeric attributes including `Nation`, and one target. The regression uses 19 substantive base attributes; `Nation` is a cohort filter, not a selected predictor. The supplied target values range from 0 to 1. The questionnaire and its scoring codebook are not included in the inspected example, so the values are described on their supplied coded scale, not translated into unprovided response categories. [E2, `001_Designer/ajsp_re_atomdat.txt`; `001_Designer/ajsp_re_dat.txt`; `002_Curator/input.txt`.]

There are two distinct analyses to keep separate. The **full-cohort analysis** builds ten models per cohort and fits a stacking ensemble to that cohort. It supports examination of fitted equations, recurring feature interactions, and model slices. The **outer-fold analysis** builds another ten models per cohort per excluded fold, fits the stack on the development rows, and evaluates the frozen stack on the excluded rows. The second analysis is the genuine held-out test workflow; the first is not.

### Directory-to-workflow map

| Folder | Main contents | What to read first |
| --- | --- | --- |
| `001_Designer` | Synthetic-record representation, attribute lookup table, preliminary descriptor report | Configuration, `ajsp_re_atomdat.txt`, `ajsp_re_dat.txt` |
| `002_Curator` | Prepared numeric table and interaction-expansion report | `input.txt`, configuration, `ajsp2gr.txt` |
| `003_Regressor/ajsp2gr` | 120 search directories, 2,400 stage reports, 120 selected model-chain files | A run's `output.txt`, then its selected chain |
| `004_Designer/Jpn020` and `Kor020` | Whole-cohort ensembles, interactions, feature scans, five training/test directory pairs | Main analysis report; each pair's `fold.txt` and test `satisfaction.se.txt` |
| `005_MLBM` | Separate machine-learning benchmark inputs, scripts, logs, predictions and SHAP outputs | `metrics_summary.tsv`, loader settings, then fold predictions |

There is **no `001_Parser` step** here. The observation records are already prepared in the representation expected by Designer. The path is therefore Designer -> Curator -> Regressor -> Designer, with a separate benchmark branch. Do not invent chemical formula parsing, atom substitution, or a materials-design objective for this example.

### Inventory of the regression work

| Analysis branch | Cohorts | Outer folds | Models per cohort/fold | Searches | Stage reports |
| --- | ---: | ---: | ---: | ---: | ---: |
| Full-cohort interpretation | 2 | Not applicable | 10 | 20 | 400 |
| Held-out outer-fold evaluation | 2 | 5 | 10 | 100 | 2,000 |
| Total | | | | **120** | **2,400** |

Every search has 20 stage summaries and a selected-model manifest. The full-cohort branch contains 20 searches and 400 stage reports; the outer-fold branch contains 100 searches and 2,000 stage reports. [E2, `003_Regressor/ajsp2gr`.]

## 3.2 Read the observation representation before the models

Designer normally looks up properties of constituents. This example reuses that interface by assigning each observation a synthetic constituent identifier and unit content. The `Atom` field in the lookup table is consequently an observation key, **not an element symbol or an atomic number**. Likewise, the Curator column `WifeName` contains numeric identifiers in this archive; it is not a predictor selected by GoodRegressor.

The first record in `001_Designer/ajsp_re_dat.txt`, with spacing simplified for display, is:

```text
1 || None || 1 1 || 0.5
```

The leading `1` identifies the record. `None` is the dummy structure label. The pair `1 1` means synthetic constituent `1` with content `1`. The final `0.5` is the supplied `satisfaction` response. Designer looks up the attribute vector for key `1` in `ajsp_re_atomdat.txt`. With exactly one unit-weight constituent, a correctly handled mean returns that observation's attribute value; it does not summarize a sample of respondents.

The small `ajsp_re_strdat.txt` declares the placeholder descriptor/site `WW`/`W` and a `None` row. `UseStrName = 0` and the attribute-based descriptor declarations mean this is not a substantive structural model. Do not interpret the repeated words `Atom`, `Materials`, `Valence`, `Cation` or `Structure` in the shared interface as psychological variables unless the example explicitly maps them to data.

### Literal fields and their roles

| Fields | Role in this example | Interpretation boundary |
| --- | --- | --- |
| `Atom` / `WifeName` | Observation lookup/row identifier | Preserve correspondence; do not intentionally add the identifier to `setX`. |
| `Nation` | Cohort selector: 1 for Jpn, 2 for Kor | Used by `SubFilterX`; not among the 200 GoodRegressor predictors. |
| `AgeWife`, `AgeHusb`, `DAgeHW` | Age-labelled attributes and the supplied difference field | Keep the existing coding and difference convention. |
| `AcaWife`, `AcaHusb` | Academic/education-labelled coded attributes | The archive does not supply their category labels. |
| `IncomeWife`, `IncomeHusb`, `DIncomeHW` | Income-labelled coded attributes and difference | Do not assign a currency or raw income unit to the fractional codes. |
| `MarrYr` | Marriage-duration-labelled attribute | Preserve the supplied numeric values and header. |
| `F_oneflesh`, `F_famprefer`, `F_desireforind`, `F_selfalloc`, `F_psychoindiv` | Five `F_*` score fields | Preserve these labels; their full instruments and scoring definitions are not supplied here. |
| `avoid`, `assert`, `yield`, `reconcil`, `integrate` | Five response-style-labelled score fields | These are input scores, not instructions to the optimization engine. |
| `satisfaction` | Target/response; also the model-group name | Use the coded scale in error metrics; it is not a probability model merely because values lie in [0,1]. |

The `L=` and `M=` strings in the Designer configuration provide LaTeX and Mathematica display labels, such as `AgeHusb` mapped to an `A_h` display. They do not replace the literal numeric field names. The unmodified reference configuration is included in `reference_inputs/example002` for matching these labels to the files.

### The numeric handoff is already prepared

`002_Curator/input.txt` has 422 lines including its header and 22 columns: identifier, 20 attributes and target. Its 421 identifiers and 20 attribute values correspond to the lookup table in the same row order. This prepared numeric table is the starting point for the Curator stage. The later Designer configuration uses the mean-only setting `2` for the singleton observation representation.

## 3.3 Curator: from 19 base predictors to 200 inputs

The Curator's `XSpace` contains all 20 attributes, including `Nation`. Its `InteractXSpace` contains the 19 non-`Nation` attributes. Thus the cohort indicator remains available for filtering without becoming part of the interaction pool. `NameDuplicateHandling = -1`, `ParseData = -`, `ParseGroup = 3`, `setY = satisfaction * *`, and `SubFilterX = -` describe the shared, unfiltered table at this stage. [E2, `002_Curator/basic_properties_GoodCurator.txt`, lines 1-14.]

The resulting selected predictor catalogue has the following composition:

| Component | Count | Example |
| --- | ---: | --- |
| Original non-cohort attributes | 19 | `AgeHusb`, `F_psychoindiv`, `assert` |
| Pair products | 171 | `(AgeHusb*assert)` |
| Admitted pair ratios | 10 | `(AgeWife/AgeHusb)` |
| Total supplied to `setX` | **200** | Literal headers are used as predictor names. |

There are 171 unordered pairs of 19 distinct base features. Ratio construction is more restricted than multiplication. In this input, the five eligible nonzero, fixed-sign fields are `AgeWife`, `AgeHusb`, `MarrYr`, `F_selfalloc`, and `F_psychoindiv`. The archive contains one oriented ratio for each of their ten pairs, in input order. It does not contain both directions for every pair, nor every possible division among all 19 features. Later symbolic construction can introduce additional transformations and composite interactions. [E2, Curator report and `setX` in `003_Regressor/ajsp2gr/Jpn020001/basic_properties_GoodRegressor.txt`; GC, interaction generation.]

`ajsp2gr.txt` is a **report containing a dataset block**, not just a rectangular matrix. The complete-cohort Regressor `input.txt` files have 421 data rows and 203 columns: `WifeName`, `satisfaction`, `Nation`, and the 200 predictors. The outer-fold versions add `Fold1` through `Fold5`, giving 208 columns. Those indicators are filters, not additional members of `setX`.

Use the extractor provided with the manual to isolate the Curator's tabular data block when creating a new handoff; do not feed the whole log to Regressor. Check the resulting header against `setX`. The saved Regressor inputs are the immediate ground truth for what the archived searches received. They should be retained with the model even when the Curator can regenerate equivalent features.

## 3.4 Regressor configuration: full-cohort and outer-fold runs

The names `Jpn020001` through `Jpn020010`, and their Kor counterparts, identify ten full-cohort searches. The additional suffix `_Fold3`, for example, identifies a search whose development data exclude outer fold 3. The `020` in the directory naming must not be substituted for the final number of model terms: the selected model may have considerably fewer than 20 terms.

### Settings that control this example

| Setting | Archived value | Meaning here |
| --- | --- | --- |
| `setY` | `satisfaction * *` | Model the supplied target with no filtering bounds. |
| `setX` | 200 literal names | The expanded pool just described, excluding identifier, cohort and fold indicators. |
| `NVariableLimit` | `-20 -19I ... -1I` | Twenty stages; decrement active additive terms while allowing hierarchical interaction expansion. |
| `PickUpPLevel` | `0.05` | Threshold for the implemented coefficient/model significance filters. |
| `TrainRatio` | `0.8` | Internal training/validation split after cohort and outer-fold filtering. |
| `NTestRobustness` | `100` | Evaluate the fitted predictions on 100 resampled subsets; not 100 refits. |
| `TargetCPUTimeBeforeSwap` | `-100 -1 !!` | A 100-second run-through time target with the faster swap branch and the encoded `!!` mode. No positive after-swap cap is specified. |
| Active transforms | 67 | Counted before the `BetaRegression` boundary; not the preprint's 109-transform study setting. |
| `BetaRegression` | `none` | Linear-coefficient symbolic regression, not a bounded beta-response fit. |

The time record does not promise a 100-second end-to-end calculation. The punctuation is interpreted by Regressor; it is unrelated to the standalone `!` that disables Designer's composition optimization. See Section 1.5 for the full record semantics. Domain restrictions of roots, logarithms and divisions still apply even though the response is bounded in the supplied data.

For the full Jpn cohort, the cohort filter is:

```text
SubFilterX(IfNone,"-";{Name/Min/Max:IfNoLimit,Min/Max=*}*NSubFilterX) Nation 1 1
```

The corresponding Kor filter uses `Nation 2 2`. The Jpn outer-fold 1 configuration additionally requires `Fold1 = 0`:

```text
SubFilterX(IfNone,"-";{Name/Min/Max:IfNoLimit,Min/Max=*}*NSubFilterX) Nation 1 1 Fold1 0 0
```

These excerpts illustrate changed records, not complete runnable configurations. Preserve the fixed record order, the full 200-name `setX` line, the complete transform block, and the terminating records from the supplied configuration.

### Actual sample counts

The full-cohort Jpn searches use 160 internal training rows and 41 validation rows; the Kor searches use 176 and 44. For outer-fold 1, Jpn has 160 development rows and splits them into 128 training plus 32 validation rows. Kor has 176 development rows and splits them into 140 plus 36. Other Jpn folds have 161 development rows and a 128/33 split. Integer counts, rather than rounded percentages alone, explain these differences.

The selected-report prediction IDs were checked against the cohort and fold input flags. In all 100 outer-fold searches, the selected model's combined training and validation rows excluded its designated outer-test rows. This verifies the stored membership, not an independent rerun of every search or a claim that every preprocessing choice was learned within each fold.

## 3.5 Read the stage summary before the final equation

`output.txt` contains the 20-stage summary and a `RECOMMEND` block. The recommended entry is the stage with the largest archived average resampled R2, rather than necessarily the final stage, the smallest equation, or the highest internal validation R2. In all 120 searches, the recommended entry agrees with the maximum of the 20 summary values and with the terminal file in its selected manifest.

For `Jpn020001`, selected stages along the trajectory are:

| Stage | Active terms | Average resampled R2 | Training R2 | Validation R2 |
| --- | --- | --- | --- | --- |
| 1 | 20 | 0.699345 | 0.677828 | 0.806252 |
| 3 | 18 | 0.715781 | 0.716075 | 0.764033 |
| 6 | 15 | 0.729790 | 0.710934 | 0.795420 |
| 7 | 14 | 0.719404 | 0.784609 | 0.500415 |
| 8 | 13 | 0.742815 | 0.779008 | 0.652128 |
| 9 | 12 | 0.735740 | 0.737718 | 0.748265 |
| 15 | 6 | 0.620513 | 0.640448 | 0.570932 |
| 20 | 1 | 0.525063 | 0.520547 | 0.545201 |

The chosen model is stage 8 with 13 active additive terms, saved as `output_satisfaction_8L13I.txt`. Its validation R2 is lower than the first-stage value, while its average resampled R2 is higher. That apparent contradiction disappears once the two criteria are distinguished: validation steers local search, while the reported robustness aggregate is used for this stage recommendation. The aggregate includes predictions for development rows used during fitting; it is not a replacement for the outer-fold results in Section 3.10.

### The 20 selected full-cohort models

| Search | Chosen stage | Terms | Average resampled R2 | Validation R2 |
| --- | --- | --- | --- | --- |
| `Jpn020001` | 8 | 13 | 0.742815 | 0.652128 |
| `Jpn020002` | 6 | 15 | 0.731188 | 0.669103 |
| `Jpn020003` | 4 | 17 | 0.730705 | 0.733155 |
| `Jpn020004` | 2 | 19 | 0.708874 | 0.746680 |
| `Jpn020005` | 7 | 14 | 0.712298 | 0.693070 |
| `Jpn020006` | 3 | 18 | 0.711572 | 0.833924 |
| `Jpn020007` | 3 | 18 | 0.698405 | 0.848335 |
| `Jpn020008` | 8 | 13 | 0.712315 | 0.733815 |
| `Jpn020009` | 3 | 18 | 0.724451 | 0.782643 |
| `Jpn020010` | 10 | 11 | 0.732896 | 0.714688 |
| `Kor020001` | 9 | 12 | 0.712223 | 0.763636 |
| `Kor020002` | 2 | 19 | 0.696470 | 0.731323 |
| `Kor020003` | 7 | 14 | 0.703180 | 0.563496 |
| `Kor020004` | 3 | 18 | 0.717993 | 0.582068 |
| `Kor020005` | 5 | 16 | 0.683172 | 0.647962 |
| `Kor020006` | 3 | 18 | 0.690716 | 0.624760 |
| `Kor020007` | 2 | 19 | 0.717689 | 0.727777 |
| `Kor020008` | 3 | 18 | 0.681834 | 0.716184 |
| `Kor020009` | 4 | 17 | 0.702164 | 0.789923 |
| `Kor020010` | 3 | 18 | 0.709091 | 0.750120 |

The selected Jpn models have 11-19 active terms, and the selected Kor models have 12-19. These are counts of the final additive components, not counts of distinct raw score fields, scalar operations, or interaction levels. A component can contain several transformed and nested base features.

![Jpn full-cohort stage trajectories](figures/example002_Jpn_depth.png)

*Figure 3.1. The ten Jpn trajectories and their mean, calculated from the archived stage-summary records. Increasing stage number decreases the allowed term count. The vertical quantity is an internal model-selection criterion, not held-out test accuracy.*

![Kor full-cohort stage trajectories](figures/example002_Kor_depth.png)

*Figure 3.2. The corresponding Kor trajectories. The curves summarize these saved searches only; they do not establish a general psychological complexity ranking between cohorts.*

## 3.6 Worked output: `Jpn020001`, stage 8

Open `003_Regressor/ajsp2gr/Jpn020001/output_satisfaction_8L13I.txt`. Long reports include search progress before the selected result, so locate `FINAL_ANSWER` rather than reading the first number that resembles a score. In this file, `FINAL_ANSWER` begins at line 711, the coefficient table at line 726, `TEST` at line 903, and `ROBUSTNESS` at line 952.

### Training, validation and resampled summaries

| Report block | Rows / repetitions | R2 | RMSE | MAE |
| --- | --- | ---: | ---: | ---: |
| `FINAL_ANSWER`, training | 160 observations | 0.779008 | 0.143884 | 0.110888 |
| `TEST`, meaning validation | 41 observations | 0.652128 | 0.178455 | 0.144479 |
| `ROBUSTNESS`, `Ave` | 100 subsets of 41 rows | 0.742815 | 0.150454 | 0.116277 |

All errors are on the supplied `satisfaction` scale. The model has 13 active terms, an intercept and 146 residual degrees of freedom: 160 -13 -1. Its report prints `full_Frequentist = 1`; this means that the implementation's significance conditions were satisfied. It does not establish causality, external validity, or absence of collinearity.

The numerical form is an intercept plus 13 coefficient-weighted components. The first printed component is `(erf-2(AgeHusb*assert))^1`, whose transform label should be looked up in the original configuration. Here `erf-2` means the `erf` row with `b = 0.01`, not an error function raised to the power minus 2. With `a = 0` and `c = 1`, this component is `erf(0.01 * AgeHusb * assert)`.

| First component field | Printed value, rounded | How to use it |
| --- | ---: | --- |
| Fitted coefficient | 6.238546 | Multiply the component value in the actual model. |
| Standardized coefficient | 2.407953 | A scale-adjusted contribution used in model reporting; not a probability or direct intervention effect. |
| Standard error | 1.026523 | The report's coefficient uncertainty estimate under its fit. |
| t value | 6.077359 | Statistic used by the implemented significance check. |
| p value | 1.10596e-8 | Printed numerical significance value, not a multiple-search causal guarantee. |
| Lower / upper coefficient bounds | 4.108264 / 8.368827 | The report's interval fields. |
| VIF | 103.715251 | The same output records substantial dependence among model components despite the significance flag. |

An individual coefficient cannot be read in isolation as the marginal relationship between one raw score and satisfaction: `AgeHusb` and `assert` can appear in other components as well. The full chain, all coefficients and the intercept are required for numerical evaluation.

### Do not confuse the two dispersion quantities

In this report, robustness R2 has `Min = 0.580058`, `Ave = 0.742815` and `Max = 0.880566`. The corresponding standard deviation is approximately `0.061903`. It measures variability of the metric across the resampled subsets.

By contrast, Designer's per-observation `ModelStdDv` or `bag0_satisfaction_stdev` describes dispersion among the individual-model predictions for one observation. It is a different quantity from the standard deviation of robustness R2. Neither quantity is automatically a calibrated prediction interval.

The report's subsequent `INTERACTION` section records the expanded representation used to reconstruct later hierarchical components. The selected report by itself is not necessarily a standalone predictor: retain its predecessor reports and its training-time transform configuration.

## 3.7 Selected manifests and portable reconstruction

A selected file such as `003_Regressor/ajsp2gr/Jpn020001.txt` is a small manifest, not the equation text. Its first records name the target, identify beta mode as `none`, and point to `basic_properties_GoodRegressor.txt`. The `chain_filenames` block then lists stage 1 through stage 8 for this particular selected model. A different search can stop at a different stage.

The supplied manifests contain absolute paths beginning with the original `/work/...` directory. Moving the archive to another computer does not make those paths valid automatically. In a working copy, replace the old root consistently in the configuration reference and every chain path; preserve the rest of the path and the stage order. Do not point a `_Fold1` manifest to the similarly named full-cohort directory merely because its files are easier to find.

`CodeOnly/Utility/find_maxR2.py` reads the recommendation in each run's summary and generates the selected-model file for Designer. Run it in a working copy with the intended search directories; the example manifests show the saved selections.

Designer also produces LaTeX and Mathematica exports of the selected models. Section 2.8 explains their component, coefficient and intercept fields.

## 3.8 Designer: full-cohort ensemble analysis

The main Jpn and Kor post-processing configurations load their respective ten selected manifests. Their `RegModelFilename` line ends with a standalone **`!`**. In this context, `!` selects collection/analysis rather than composition optimization. Preserve it: replacing synthetic observation identifiers with other identifiers would not constitute a valid intervention on a person.

The important input combination is:

```text
MaterialsFilename(...)/HeadExists(...)/PostTagNames(...) ajsp2gr_dat.txt 0 satisfaction
RegModelFilename(...) Jpn020001.txt ... Jpn020010.txt !
AnalyzeInteractionLevel 2 10
SkipMinMaxKurtoInAnalyzeInteraction(Yes1No0) 2
```

The ellipsis above abbreviates unchanged labels and the member list for discussion only. Use the full supplied reference input to run the program. `AnalyzeInteractionLevel` and the mean-only setting are distinct controls: the first requests interaction analysis, while the second ensures the observation-attribute representation is interpreted using the mean-only branch.

Other shared-interface records remain present, including `TargetVReg = 1`, `TargetStdDvReg = 0.5`, `ZDistWeight = 1`, and `UseOffset = 1`. Their presence does not mean this archive has optimized a respondent to a desired satisfaction value. The saved output used here is the direct model/ensemble analysis. In particular, `TargetStdDvReg` is not the optional outer-test `TestModelStDevLimit` record described below.

### Prediction table: one row per observation

The full-cohort `__R__.ajsp_re_anaoutput.txt` has 201 Jpn or 220 Kor data rows and 13 columns:

| Column group | Number | Meaning |
| --- | ---: | --- |
| `reg0_satisfaction` through `reg9_satisfaction` | 10 | The individual selected-model predictions in the manifest-list order. |
| `bag0_satisfaction` | 1 | The final fitted stacking-ensemble prediction. |
| `bag0_satisfaction_stdev` | 1 | Model dispersion for that observation. |
| `satisfaction` | 1 | Supplied observed target for comparison. |

The prediction table does not carry an independent identifier column. Reattach IDs from the matching `ajsp2gr_dat.txt` selection in the same verified order. For a fold, apply the fold's row-selection rule before attaching IDs; joining the first 41 original rows to fold 1's 41 predictions would be wrong.

The detailed analysis report also contains a `ModelName` metric table and a final `StackingEnsemble` coefficient block. These are different logical tables within the same text file. `R2`, `PearsonR` and `PearsonR2` are separate fields. Use residual-based `R2` when quoting the reported R2; do not substitute squared correlation.

### Fitting the ensemble is not taking a simple average

The stack has an intercept and ten coefficients. The Jpn final intercept is approximately -0.0489841; its member weights include positive and negative values, for example 0.3574244 for member 0 and -0.0582911 for member 3. Thus `bag0_satisfaction` is not the unweighted mean of `reg0` through `reg9`. Member order and the intercept matter. [E2, `004_Designer/Jpn020/ajsp_re_anaoutput.txt`, final 11-line block.]

`acclevel_0` refers to fitting the ensemble from its first member, while `acclevel_9` refers to all ten. This index is **not** Regressor's interaction-depth stage. Even a one-member stack may include a fitted intercept and calibration coefficient, so it should not automatically be treated as numerically identical to the raw member prediction.

| Full-cohort ensemble | Labeled rows | R2 | RMSE | MAE | Invalid |
| --- | --- | --- | --- | --- | --- |
| Jpn | 201 | 0.793644 | 0.138724 | 0.106879 | 0 |
| Kor | 220 | 0.766905 | 0.123113 | 0.096828 | 0 |

These are whole-cohort fitting/evaluation metrics. All observations in the respective cohort participate in this branch's construction or stack fitting, and invalid count is 0 in both final full-cohort summaries. They must not be captioned as independently tested predictive performance.

![Full-cohort ensemble accumulation](figures/example002_stacking.png)

*Figure 3.3. Whole-cohort R2 as members are accumulated. The figure describes ensemble fitting, not a learning curve on an untouched test set. The actual held-out fold results are reported separately with their retained coverage.*

## 3.9 Interaction tables and model-slice outputs

### Read `.zc.txt` as an equation-structure summary

`bag0_satisfaction.zc.txt` summarizes recurring base-feature sets across the selected models. `InteractionLevel = 1` denotes one base feature; level 2 denotes a pair co-occurring within a model component. `InHowManyModelsAppear` counts models containing that feature set, whereas `HowManyAppear` counts occurrences across components. These are not the same count.

The initial table contains the following illustrative entries, rounded from the files:

| Cohort | Feature set | Level | Models | Occurrences | Printed Ave(abs Z) |
| --- | --- | --- | --- | --- | --- |
| Jpn | `AgeHusb` | 1 | 10 | 40 | 0.720953 |
| Jpn | `assert` | 1 | 10 | 51 | 0.540328 |
| Jpn | `yield` | 1 | 10 | 47 | 0.444657 |
| Jpn | `MarrYr-assert` | 2 | 10 | 23 | 0.103017 |
| Jpn | `IncomeHusb-yield` | 2 | 9 | 31 | 0.508124 |
| Kor | `avoid` | 1 | 10 | 49 | 0.552505 |
| Kor | `AgeHusb` | 1 | 10 | 37 | 0.435802 |
| Kor | `F_psychoindiv-assert` | 2 | 10 | 16 | 0.075733 |

For example, Jpn `AgeHusb` appears in all ten models, with 40 occurrences across terms. That means repeated representation in this fitted symbolic ensemble; it does not by itself establish that changing age causes a particular satisfaction change. A frequent pair can have a smaller aggregated coefficient statistic than a less frequent pair. The source's chain-selection logic and its coefficient aggregation should not be replaced by an improvised sorting rule.

The headings `Min(|Z|)`, `Ave(|Z|)` and `Max(|Z|)` refer to the implementation's coefficient statistics. As discussed in Section 2.8, signed stacking weights enter the aggregation; do not assume the `Ave(|Z|)` label is a conventional nonnegative feature-importance definition. Additional chain sections later in the file repeat some entries. Use the initial table or a specifically identified chain section rather than adding every repeated row together. [GD, lines 6006-6037 and interaction-chain reporting.]

### One-variable and two-variable slices

A file such as `_bag0_satisfaction.F_oneflesh.txt` contains 101 mesh points, plus its header. `_bag0_satisfaction.F_oneflesh-reconcil.txt` in Jpn contains 10,201 mesh rows, corresponding to a 101-by 101 grid. Those counts are **not additional respondents or test observations**. The separate `bag0_satisfaction.onevtest.txt` also contains analysis sections; `test` in that filename does not mean an independent dataset.

| Field | What is plotted or summarized |
| --- | --- |
| `vF_oneflesh`, or another `v...` field | Normalized mesh coordinate. |
| `F_oneflesh`, or the corresponding bare field | The coordinate in the supplied numeric coding. Use this for an unnormalized axis. |
| `satisfaction` | Ensemble prediction at the displayed coordinates, with other features at reference values. |
| `ModelStdDv` | Dispersion across the member models at the mesh point. |
| `Nsample`, `Min`, `Ave`, `+-`, `Max` | Local observed-data counts and summary values associated with the mesh location. |
| `Nsample(R)`, `Min(R)`, `Ave(R)`, `+-(R)`, `Max(R)` | Corresponding summaries of regressed responses. |
| `MeanAE`, `MaxAE` | Local absolute-error summaries produced by the routine. |

For the first Jpn `F_oneflesh` mesh row, the coordinate is 0, the surface prediction is approximately 0.0768711, `ModelStdDv` is 0.0956302, and `Nsample = 2`. The local observed mean is 0.25; the local mean regressed value is approximately 0.152199. These differ because a model slice with other features fixed is not the same operation as averaging the actual observations or averaging their individually evaluated predictions. The next grid point has `-` in the local-data fields: that indicates no reported local summary, not a zero measurement. [E2, `004_Designer/Jpn020/_bag0_satisfaction.F_oneflesh.txt`, lines 1-3; GD, feature-map writer near 6321.]

These are **fixed-reference feature scans**, not necessarily population-averaged partial dependence. Varying `AgeHusb` while fixing `AgeWife` and `DAgeHW`, or varying a fractional category code continuously, can produce combinations outside the represented observations. Inspect local support, the supplied scale coding and numerical validity before interpreting a curve. The archive does not establish that these slices are feasible interventions or causal effects.

## 3.10 Genuine outer tests: select folds, freeze the stack, report coverage

This example provides a second layer of evaluation beyond the legacy `TEST` block. Its outer-test rows are excluded from the corresponding Regressor searches and from stack fitting. The `_Fold1` searches feed `Fold1_train` and `Fold1_test`; the full-cohort searches do not substitute for them.

### Decode `fold.txt`

The Jpn fold 1 development directory contains:

```text
NFolds TakeFold Train0Test1 Straight0Rotate1
5 1 0 1
```

The corresponding test directory contains:

```text
NFolds TakeFold Train0Test1 Straight0Rotate1
5 1 1 1
TestModelStDevLimit= 0.5
```

`NFolds = 5` and `TakeFold = 1` select the first outer fold. `Train0Test1 = 0` keeps its complement; `1` keeps the fold itself. `Straight0Rotate1 = 1` uses rotating, row-order assignment. With zero-based row position `i` within the cohort, the fold is `(i % 5) + 1`. This is not random shuffling and is not computed from the numeric identifier: IDs can have gaps. [GD, fold reader near 9378-9403 and row-selection code 6641-6656.]

Jpn test-fold sizes are 41,40,40,40,40; Kor has 44 in each fold. The development sets have 160 or 161 Jpn rows and 176 Kor rows. Within each development set, ten selected symbolic models are constructed with their own internal training/validation operation. Designer then fits the stack using those development rows.

### How the archived script keeps test coefficients fixed

`run_traintest.sh` prepares per-fold training directories using the fold-specific manifests. It then builds test directories from the training setup and copies the last 11 lines of the training analysis report to `satisfaction.se.txt`: one intercept and ten weights. The test evaluation loads these coefficients rather than fitting another stack to its held-out observations. [E2, `004_Designer/Jpn020/run_traintest.sh`, lines 456-467, and Kor counterpart.]

The test and training directories contain canonical names such as `Jpn020001.txt`. Their **contents**, not just the copied basename, identify `_Fold1` model dependencies. The audit checked that all ten test coefficient files match the final block in their respective training reports and that the selected Regressor memberships exclude those test IDs. Keep the member order unchanged when copying a coefficient file.

This supports the intended stored separation of model fitting and outer evaluation. It does not prove that every data-preparation choice was fitted inside the folds: a common curated feature representation exists before fold filtering. Record that representation and its provenance rather than claiming an uninspected fully nested preprocessing procedure.

### The dispersion screen changes the denominator

The test-side `TestModelStDevLimit = 0.5` is active. The evaluation branch inserts a nonfinite ensemble result when model dispersion exceeds this limit; genuinely nonfinite numerical predictions can also be invalid. The metric table excludes invalid predictions. The result is **conditional test performance on retained rows**, not necessarily performance on all rows assigned to that fold. [GD, lines 7126-7142.]

| Cohort / fold | Assigned | Retained | Rejected | R2 | RMSE | MAE |
| --- | --- | --- | --- | --- | --- | --- |
| Jpn / 1 | 41 | 40 | 1 | 0.562357 | 0.204776 | 0.173363 |
| Jpn / 2 | 40 | 40 | 0 | 0.386462 | 0.230755 | 0.187063 |
| Jpn / 3 | 40 | 38 | 2 | 0.638999 | 0.161933 | 0.123334 |
| Jpn / 4 | 40 | 40 | 0 | 0.311118 | 0.266086 | 0.211849 |
| Jpn / 5 | 40 | 36 | 4 | 0.307738 | 0.238199 | 0.194688 |
| Kor / 1 | 44 | 43 | 1 | 0.412812 | 0.193535 | 0.152152 |
| Kor / 2 | 44 | 42 | 2 | 0.183012 | 0.220406 | 0.173706 |
| Kor / 3 | 44 | 38 | 6 | 0.619740 | 0.144766 | 0.113808 |
| Kor / 4 | 44 | 44 | 0 | 0.262663 | 0.228516 | 0.167048 |
| Kor / 5 | 44 | 43 | 1 | 0.213300 | 0.208776 | 0.170210 |

The metrics in this table were independently recalculated from the stored prediction rows and agree with the archived summaries. The displayed R2 is the `R2` field, not `PearsonR2`. `Rejected` here means excluded from the final metric calculation, encompassing the archive's dispersion and numerical-invalid outcomes; it is not a manually removed respondent from the original dataset.

| Cohort | Total retained | Coverage | Mean fold R2 | Mean fold RMSE | Mean fold MAE |
| --- | --- | --- | --- | --- | --- |
| Jpn | 194/201 | 96.52% | 0.441335 | 0.220350 | 0.178059 |
| Kor | 210/220 | 95.45% | 0.338305 | 0.199200 | 0.155385 |

The mean row averages the five fold-specific metrics with equal fold weight; it is not R2 recomputed after pooling all predictions. Coverage is computed from total retained and assigned observations. Report both. A statement such as "Jpn outer-test R2 is 0.4413" is incomplete without the 0.5 dispersion rule and 194/201 retained coverage.

Do not recalibrate `satisfaction.se.txt` on test responses, tune the dispersion cutoff on the same test outcomes without declaring that selection, or replace invalid predictions with zero. The companion evidence retains the invalid counts and observation IDs. Diagnostic reconstruction from the saved member columns also shows that rejected predictions are not harmless blanks; finite rejected values can differ greatly from the observed scale, and some member outputs are nonfinite. A full-coverage comparison requires an explicitly specified failure-handling policy, not silently dropping difficult rows.

### Safe rerun boundary

The archived helper script contains machine-specific source paths and expects a usable `GoodDesigner.x`. It also has `RECREATE_TRAIN_DIRS = 1` and `RECREATE_TEST_DIRS = 1`, enabling replacement of existing fold directories. These are not read-only defaults. Do not run that script in the original archive to "check" the results. Prepare a separate working copy, inspect the path and directory-recreation controls, compile the appropriate executable, validate manifests and coefficient order, and only then launch the intended operation.

## 3.11 Read the benchmark branch as it was actually run

`005_MLBM` is a separate implementation, with one `satisfaction.xlsx` workbook for each cohort. The Jpn workbook has 201 observations and the Kor workbook 220; both have 22 columns. Their identifiers and target order match the corresponding cohort tables. Eight benchmark models have saved five-fold prediction files: Ridge, ElasticNet, RandomForest, XGBoost, LightGBM, MLP, PhySO and EQL. Do not add a PySR result merely because a script or environment contains optional Julia/PySR-related code.

### Inputs are not identical to GoodRegressor's inputs

The benchmark's `ID_LIKE_COLS` is `['CompositionName']`. That name is absent from these workbooks, whose identifier is **`WifeName`**. The loader therefore retains `WifeName`, `Nation`, and the 19 substantive attributes, giving **21 numeric predictors**. This is confirmed by both the source loader and the saved SHAP metadata's feature list. GoodRegressor, in contrast, selects 200 predictors formed from the 19 substantive attributes and excludes `WifeName` and `Nation`. [BM, `benchmark_v14.py`, line 112 and lines 1399-1465; E2, `005_MLBM/ajsp2gr_Jpn/out/shap/satisfaction_shap_RandomForest_metadata.json`, and Kor counterpart.]

Within each cohort `Nation` is constant, but `WifeName` is a varying row identifier. Its inclusion is an actual input difference that should be reviewed; this manual has not silently removed the column or recalculated the benchmark. Combined with GoodRegressor's conditional prediction coverage, it prevents treating the archived numbers as a fully matched, identical-input comparison without further qualification.

### Outer folds and inner selection

The benchmark script uses its `RoundRobinKFold` for the five outer folds. For estimators with an external hyperparameter space, it uses four-fold round-robin `RandomizedSearchCV` within the development data, with R2 scoring. Models whose specification has `param_space = None` use their own internal-search fit instead of that same external randomized-search loop. Numeric preprocessing is part of the model pipeline. [BM, `RoundRobinKFold` nearline 50; `nested_cv_oof` fromline 1218; `build_preprocessor`.]

The launcher schedules eight models times five folds and sets a configurable `n500` hyperparameter regime. The label `n500` is not the actual number of respondents: the workbooks contain 201 and 220 rows. Cluster partition names, module loads, Python setup and resource counts in the archived shell script describe its original execution environment, not universal hardware requirements for GoodRegressor.

### Archived five-fold summaries

| Model | Jpn R2 | Jpn RMSE | Jpn MAE | Kor R2 | Kor RMSE | Kor MAE |
| --- | --- | --- | --- | --- | --- | --- |
| Ridge | 0.34846 | 0.24156 | 0.17932 | 0.39706 | 0.19376 | 0.13888 |
| ElasticNet | 0.33686 | 0.24296 | 0.18712 | 0.40336 | 0.19288 | 0.13624 |
| RandomForest | 0.26590 | 0.25224 | 0.18136 | 0.38030 | 0.19602 | 0.13926 |
| XGBoost | 0.06002 | 0.29014 | 0.21826 | 0.20562 | 0.22392 | 0.16504 |
| LightGBM | -0.14168 | 0.32090 | 0.23628 | 0.15252 | 0.22984 | 0.16410 |
| MLP | -0.13638 | 0.31744 | 0.23834 | 0.30368 | 0.20940 | 0.15860 |
| PhySO | -0.01892 | 0.30424 | 0.25872 | 0.02296 | 0.25018 | 0.19844 |
| EQL | -0.96052 | 0.41092 | 0.30324 | -0.91168 | 0.34822 | 0.25166 |

These values reproduce `metrics_summary.tsv`, rounded here to five decimals. The companion postprocessor takes the final matching metric line from each fold log, then computes means and sample standard deviations; the log metrics are printed to limited precision. Recalculating from the saved CSV predictions can therefore differ slightly in the last decimal. Each model has five saved fold results, and those prediction files cover the complete assigned cohort rows in this archive.

### File dictionary for benchmark analysis

| File or directory | What it contains | Practical check |
| --- | --- | --- |
| `satisfaction.xlsx` | Numeric input data for this benchmark | Inspect the actual header against `ID_LIKE_COLS` before running. |
| `model_logs/*_fold*.log.txt` | Per-model/per-fold execution and metric logs | Keep one intended completed log per model/fold when aggregating. |
| `out/satisfaction_predictions_*_fold*.csv` | Two-column true/predicted responses for one outer fold | Restore IDs with that fold's row-order rule, not a sequential whole-cohort join. |
| `metrics_summary.tsv` | Means, sample standard deviations and fold counts | Check `n_folds_found = 5`; duplicated logs can be double-counted by the supplied merger. |
| `predictions_wide.tsv` | Per-model true/predicted columns after concatenating fold files | Its row order is fold concatenation, not necessarily the original workbook order. |
| `out/best_params` | Hyperparameter records from the fold fits | Retain them with the model/version and fold definition. |
| `out/shap` | Final-refit explanations and their metadata | These explain full-data refitted models, not a new independent test. |

### SHAP outputs: read the metadata, not only the picture

The archive contains SHAP CSVs, metadata and summary images for RandomForest, XGBoost and LightGBM in each cohort. The `--final-refit-only` mode reads saved fold hyperparameter records and fits a model to the full cohort; this step does not rerun the outer cross-validation. SHAP therefore belongs to a full-data explanation branch, separate from the held-out predictions.

The saved metadata explicitly says that SHAP values explain the **inner tree-regressor output after the target transformation**, not the inverse-transformed satisfaction scale. Its `sample_index` is a zero-based row position after loader filtering, not `WifeName`. Preserve that scale and identifier distinction when comparing a SHAP plot to a symbolic-model slice. The feature list also visibly contains `WifeName`; any interpretation of this benchmark should acknowledge that input before attributing the entire plotted importance to substantive score fields.

# Appendix A. Output-file dictionary

The table covers the recurring outputs present in Example001 and the principal temporary files encountered in the source. A prefix such as `BCCYield` or `Yield` changes with the input/task.

| File or section | How to use it |
| --- | --- |
| Parser `*_log.txt` | Parsing diagnostics; reconcile failures and successful material counts. |
| Parser `*_output.txt` | Intermediate parsed formula/block representation and aggregate diagnostics. |
| Parser `*_dat.txt` | Parsed-material handoff plus appended target/postscript values. |
| Designer `*_anaoutput.txt` | Descriptor/material table, and with loaded models, additional performance summaries. Identify each logical table. |
| Designer `__R__.*_anaoutput.txt` | Model/ensemble predictions and dispersion fields aligned to material row order. |
| Designer `*_designeroutput.txt` | Main model/stacking report and composition-search progress when enabled. |
| Curator report, e.g. `BCCYield.txt` | Preamble, rectangular `DATASET`, then summaries; not a pure TSV throughout. |
| Regressor `input.txt` | Prepared numeric table with the literal curated header. |
| Regressor `output_<target>_<stage>.txt` | Full search-stage record, final fit, training/validation rows, robustness and interactions. |
| Regressor `FINAL_ANSWER` | Locate the final chosen fit within one stage, not necessarily the selected depth for the ensemble. |
| Regressor `TEST` | Legacy name for validation data/results. |
| Regressor `ROBUSTNESS` | Fixed-model resampling diagnostic, with minimum, mean, maximum and standard deviation of the metrics. |
| Small `005001.txt`-style manifest | Target/configuration and preceding-stage dependencies for one selected model. |
| Target `.se.txt` | Saved progressive stacking coefficients in the loaded model order. |
| Target/index `.latex.txt` | LaTeX representation of the symbolic components, coefficients and intercept. |
| Target/index `.mathematica.txt` | Mathematica representation of the symbolic components, coefficients and intercept. |
| `bag*.zc.txt` | Descriptor/interacting-set occurrences and implementation-specific coefficient statistics. |
| `bag*.onevtest.txt` | One-variable/model-analysis information; not an independent test set. |
| `_bag*.<feature-or-pair>.txt` | Feature-scan mesh, predictions, dispersion and nearby-data summaries. |
| `<species>.nb.txt` | Optional explicit replacement-neighbor input; may alter apparent configuration behavior. |
| `fold.txt` | Optional row/fold analysis selection; not a complete nested-validation pipeline. |
| `*.temporary.txt`, `*.reg`, `*.pr.txt`, `smt`, worker files | Intermediate coordination/working artifacts; isolate job directories and retain final dependencies before cleanup. |

## A.1 Additional files in Example002

| File / pattern | Purpose | Distinction to retain |
| --- | --- | --- |
| `ajsp_re_atomdat.txt`, `ajsp2gr_atomdat.txt` | Observation attribute lookup | Synthetic observation keys, not chemical elements. |
| `ajsp_re_dat.txt`, `ajsp2gr_dat.txt` | Singleton record, unit content and target | Preserve the cohort-specific row order. |
| `Jpn020001.txt`, `Kor020001.txt` | Selected whole-cohort model-chain manifests | Require original configuration and all predecessor stages. |
| `*_Foldk.txt` | Model selected while excluding outer fold k | Verify file contents after renaming to canonical member names. |
| `fold.txt` | Designer development/test membership and optional dispersion limit | Rotating assignment uses cohort row position. |
| `satisfaction.se.txt` in `Foldk_test` | Frozen development-fitted stack | One intercept plus ten ordered weights; do not refit on test responses. |
| `__R__.ajsp_re_anaoutput.txt` | Whole-cohort model predictions | Full-cohort fit, not outer-test prediction. |
| `Foldk_test/__R__.ajsp2gr_anaoutput.txt` | Held-out predictions with validity/dispersion behavior | Report assigned and retained counts. |
| `metrics_summary.tsv` | Separate benchmark fold-mean metrics | Its predictors and coverage differ from the archived GoodRegressor setup. |
| `out/shap/*_metadata.json` | Explanation scale, feature list and row indexing | Full-refit, transformed-target explanations; not outer-test metrics. |

# Appendix B. Source map and documentation policy

## B.1 Where the implementation is defined

| Topic | Primary source location |
| --- | --- |
| Parser parameter order, atom blocks | GP 217-265 |
| Curator parameter order | GC 103-234 |
| Composition-based selection | GC `pick_data`, approximately 249-385 |
| Designer configuration and descriptor declarations | GD 3430-4050; descriptor continuation near 4107 |
| Atomic-statistic weighting and zero handling | GD `get_mamssk`, 2019-2093 |
| Ensemble prediction and model dispersion | GD `calc_bag_chain`, 3209-3233 |
| Model-chain reading | GD `read_reg_chain`, from 4944 |
| Interaction aggregation | GD 6006-6037 |
| Group-indexed target access | GD material-processing loop near 7783 |
| Equation-string formatting | GD `chk_eql`, `get_powerstring`, `get_powerstringM` and equation-export writers |
| Regressor parameter order / transform boundary | GR 1290-1563 |
| Scalar numerical functions | GR `get_v`, 900-940 |
| Fixed-model robustness and standard-deviation output | GR, robustness routine and final square-root conversion |

### Additional implementation references for Example002

| Topic | Primary source location |
| --- | --- |
| Mean-only mode and zero retention | GD 2019-2093 and 9455-9470 |
| Rotating fold membership | GD 6641-6656; fold reader near 9378-9403 |
| Test-side dispersion rejection | GD 7126-7142 |
| Feature-scan output headings | GD near 6321; nearby observed/regressed-column lookup |
| Transfer of frozen stack to test directory | E2, `004_Designer/Jpn020/run_traintest.sh`, 456-467 |
| Benchmark identifier exclusion list | BM line 112; numeric loader1399-1465 |
| Outer and inner benchmark split logic | BM `RoundRobinKFold`; `nested_cv_oof` from 1218 |
| SHAP scale and actual feature list | E2, both benchmark `out/shap/*_metadata.json` files |

## B.2 Where the preprint supplies the conceptual explanation

The main-text Methods on pp. 24-27 explain lexicographic sampling, swap, transit, pick and repeated model construction. SI S3-S9 explains the five-step workflow and ensemble post-processing. SI S10 lists the study's scalar transforms; it is not a substitute for counting the active transforms in an actual configuration. SI S30-S38 explains interaction sets/chains and feature scans. SI S41-S43 explains offset-anchored materials prediction and standardized target distance.

The outer train/validation/test protocol is described on main-text p. 9; the full-data final-ensemble construction without a separate test set is described on p. 12. The manuscript and its illustrative older figures contain legacy `test` usage in some individual-model contexts. This manual follows the author's clarification for program outputs while preserving the genuinely excluded outer-test meaning.

## B.3 Source-versus-documentation distinctions

Numerical examples refer to the supplied saved calculations. Input descriptions incorporate the revised settings and definitions documented here. When adapting an example, retain its input tables, selected-model files and stacking coefficients alongside the numerical results.

Use the filenames and function names in this source map to locate the relevant operations. Line numbers identify the reference source and can move when code is edited. Keep the code version, active transformation catalogue and configuration with any newly generated results.

# Appendix C. Calculation-speed features

GoodRegressor combines three acceleration mechanisms that address different sources of cost. **Cholesky updates** reduce repeated linear-algebra work when nearby candidate models are fitted. **On-the-fly learning** limits the feature/transform trials explored in a refinement pass according to the evolving fit. **Active MPI** reallocates completed ranks to ranks with unfinished work. These mechanisms complement the time-limited, lexicographically ordered initial search rather than replacing its model-selection objective.

The starting presets in Section 1.5, `-100 -1 !!` and `-1000 -1 !!`, enable the faster swap branch and on-the-fly mode. Active MPI is part of the parallel coordination workflow; `!!` is not a separate switch that turns MPI on.

## C.1 Cholesky decomposition and incremental updates

For a candidate model whose symbolic components form the columns of a design matrix X, linear-coefficient fitting solves a least-squares problem. The implementation standardizes the relevant quantities and factorizes the Gram matrix:

```text
G = transpose(X) * X = L * transpose(L)
G * beta = transpose(X) * y
```

The Cholesky factor L permits the coefficient solution through triangular systems. During a swap or scalar-transform trial, much of the candidate design matrix is unchanged. Rebuilding and refactorizing the entire Gram matrix for every such trial repeats work unnecessarily.

The implementation therefore caches the standardized design matrix and its Cholesky factor. If column j changes from x_j to x_j + d, the new Gram matrix can be written as:

```text
u = transpose(X) * d
s = transpose(d) * d
G_new = G + u * transpose(e_j) + e_j * transpose(u)
          + s * e_j * transpose(e_j)
```

Here e_j selects the changed column. The code represents this correction using a two-dimensional symmetric matrix, decomposes that small matrix, and applies rank updates or downdates to the existing factor. With several changed columns, it can apply successive column updates. When too many columns change, or cached information is unsuitable, it calculates a fresh factorization instead. [GR, `cholesky_update_one_column`, `cholesky_update_multiple_columns`, `fit_any`, approximately 2595-2869.]

This reuse accelerates the repeated fitting of closely related candidate equations. The faster-swap option is selected by a **negative first value** in `TargetCPUTimeBeforeSwap`; the code uses the magnitude as the run-through time target. The numerical fitting path also contains fresh-factorization checks. Cached updates do not remove the need for a numerically usable design matrix or make singular, redundant components harmless.

The practical benefit depends on the dataset, term count and how much of the candidate model changes. This mechanism reduces repeated work; it does not imply a fixed speed-up factor or an exhaustive search of all possible models.

## C.2 On-the-fly learning: adapt the explored feature and transform scope

The `!!` preset sets `useSomeShotMode = true`. In this mode, a refinement pass does not indiscriminately explore every feature position and every associated scalar-transformation trial. Instead, the scope of trials is adjusted while the model is being trained.

Descriptor swaps and scalar transformations use significance-ordered component lists. Descriptor replacement starts from the less significant side, while scalar-function switching starts from the more significant side. The current coefficient statistics establish a working boundary in that ordered list. As improved candidates are accepted and the coefficient statistics change, the boundary is updated for subsequent trials. The source records these bounds through variables such as `NNend` and log entries containing `NNendUpdated`. [GR, `swap_descriptors` and `switch_functype` branches using `useSomeShotMode`.]

This is **on-the-fly learning of the exploration scope**, not the addition of a separate neural network. The evolving model helps determine where to spend the next refinement effort. Features and function types can be considered in a restricted set of combinations in one pass rather than sweeping the complete feature-by-transformation space regardless of the current fit.

The original `Functypes` catalogue remains the definition of the available transformations. A shortened pass does not redefine that catalogue or permanently remove unvisited functions from the saved input. It changes which trials are attempted during refinement. This distinction matters when comparing runs with the same catalogue but different speed settings: their visited candidate sets, and hence selected equations, can differ.

Useful log phrases include:

```text
swap_descriptors: useSomeShotMode
switch_functype: useSomeShotMode
NNend
NNendUpdated
```

For ordinary use, the two presets in Section 1.5 are the simplest entry point. The 100-second version gives a shorter initial search budget; the 1000-second version allocates more time to the initial ordered search. Both retain the on-the-fly refinement mode.

## C.3 Active MPI: LUHMANN followed by HEGEL

The parallel workflow has two complementary phases. The names **LUHMANN-AUTOPOIESES** and **HEGEL** are literal labels in the implementation and logs.

**LUHMANN-autopoieses: independent regional work.** The lexicographically ordered model space is distributed among MPI ranks. Each rank explores its assigned region and develops its local candidate through search/refinement. During the independent refinement phase, the ranks work on their own models. Ranks need not finish at the same time because their candidate models and remaining trials can differ.

**HEGEL-MasterSlave: completed ranks help unfinished ranks.** At coordination points, the program gathers completion status. When some ranks have completed their own work while others have not, completed ranks do not simply remain idle. The unfinished ranks become masters of cooperative groups, and completed ranks are assigned to those groups as helpers, called slaves in the code and its logs. They take portions of the remaining descriptor-swap or function-switch trials, working with the appropriate master's candidate state.

The coordinating rank distributes available helpers among the unfinished groups. Within each group, candidate evaluations are shared and the search continues for the master's model. Group membership is reconsidered as ranks finish. The coordination cycle terminates when all ranks have completed the relevant work. [GR, completion-state gathering, `Hegel` grouping and LUHMANN/HEGEL dispatch, approximately 7780-8020.]

```text
Initial lexicographic distribution
    -> each rank explores and refines its regional candidate
    -> LUHMANN: independent local work
    -> gather completion status
    -> HEGEL: completed ranks assist unfinished masters
    -> regroup as work completes
    -> finish when all ranks are done
```

This active redistribution reduces the idle time that would otherwise arise when a fixed MPI partition leaves some ranks waiting for slower ranks. It concerns the remaining search/refinement tasks; it need not mean that the program moves an entire original lexicographic interval between ranks.

The log can show `LUHMANN-AUTOPOIESES`, `HEGEL`, and `Master` / `Slave(s)` group listings. These messages describe coordination, not a new target variable, model family or validation split. Run each independent job in its own working directory and provide the shared files required by the MPI workflow.

## C.4 Choosing the time record

| Record value | Intended starting use | Enabled speed options |
| --- | --- | --- |
| `-100 -1 !!` | Shorter initial exploration | Faster swap / Cholesky reuse branch and on-the-fly scope selection |
| `-1000 -1 !!` | Larger initial exploration budget | The same options with a longer run-through target |

Use one of these values for the ordered `TargetCPUTimeBeforeSwap` record. The second field does not impose a positive refinement-time cap, so the whole job can take longer than 100 or 1000 seconds. The dataset, depth schedule, rank count and subsequent Designer interaction analysis also affect the total time.

For Designer post-processing, an additional practical choice is **`AnalyzeInteractionLevel 1 1`** when detailed interaction chains are not needed. That setting reduces the requested interpretation work; it is separate from the Regressor's three acceleration mechanisms.

**End of manual.**
