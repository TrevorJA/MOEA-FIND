# MOEA-FIND Style and Syntax Guide

*Created 2026-04-14 after the user flagged non-academic patterns in
the manuscript and SI drafts. Based on the user's explicit rules
combined with an opus-driven evidence review of nine published
papers in the project literature folder (Borgomeo 2015, Zaniolo
2024, Wheeler 2025, Kirsch 2013, Reed 2013, Hadka and Reed 2015,
Hadjimichael 2020, Quinn 2018, Bonham 2024). Durable. All future
edits to `manuscript_main_draft.md` and `supporting_info_draft.md`
must satisfy every rule below before commit.*

## 1. Scope and priority

This document governs `drafts/manuscript_main_draft.md`,
`drafts/supporting_info_draft.md`, and any figure captions,
table captions, or abstract text destined for submission to *Water
Resources Research*. It does not govern planning notes, design
decision documents, or the shell-versus-interior diagnostic note,
which may use lists, inline code formatting, and informal
constructions where that serves clarity for the author's own
reference.

When a rule in this guide conflicts with a rule anywhere else in
the project, this guide wins. When the guide is silent on a
question, the default is to imitate Borgomeo (2015, WRR), Zaniolo
(2024, EMS), or Wheeler (2025, JHE), in that order of precedence.

## 2. Rules the user stated explicitly

The user directed, on 2026-04-14, that the following patterns be
eliminated from every subsequent edit.

Rule 2.1. Do not use bold to highlight topic phrases at the start
of body paragraphs, and do not start paragraphs with incomplete
bold leads such as "**Pareto front equals feasible set.** The ...".
Every body paragraph opens with a complete declarative sentence.
Bold is reserved for the paper title, figure and table labels as
required by journal house style, and for the word `Figure` or
`Table` when cited in the text.

Rule 2.2. Do not use triple-backtick code blocks, inline
backtick-quoted symbols, or code-style notation anywhere in body
text. Mathematical content is rendered in LaTeX. Symbols in inline
text are italicised. Display equations are numbered and referenced
as `Eq. (n)`.

Rule 2.3. Do not use informal language. This includes contractions,
rhetorical questions, colloquialisms, authorial asides, scare
quotes for emphasis, and hedging adverbs that do not carry
information ("arguably", "essentially", "effectively",
"ultimately", "quite", "rather", "fairly").

Rule 2.4. Define every technical term on first use, and avoid
vague or non-technical terminology. Terms of art that appear once
in the manuscript without a definition, a citation, or a
parenthetical gloss are disallowed. The earlier draft used
"image of a map" without defining what a map is, which is the
canonical example of the failure mode this rule exists to prevent.
Replace such terms with language drawn from the water resources,
systems engineering, and decision-making-under-deep-uncertainty
literatures whenever possible. When a mathematical term has no
established plain-language equivalent in those literatures, define
it explicitly with a one-sentence gloss at first use and cite a
standard reference.

Rule 2.5. Do not use bulleted or numbered lists in body text.
Continuous prose is required in the Introduction, Methods,
Results, Discussion, and Conclusions. Lists are permitted only for
the genuinely structured content where they are standard practice
in the predecessor papers, specifically (i) stepwise algorithm
descriptions in the Methods, where each step is a short numbered
item interleaved with a display equation, and (ii) notation tables
or figure-panel enumerations.

Rule 2.6. Section and subsection headers must be informative,
objective, and descriptive. They name the content of the section,
not a claim about it. No rhetorical titles, no questions, no
declarative statements, no section titles that preview a result.
Headers are short (typically two to seven words).

Rule 2.7. The manuscript must build narrative and methods logically
and sequentially. A term introduced in section N cannot be relied
upon in section M for M less than N. A result reported in section
N cannot be the basis for an argument in section M for M less than
N unless the earlier section explicitly forward-references the
later one. Where forward references are unavoidable, they are
parenthetical and brief.

Rule 2.8. Avoid superfluous or over-zealous language. The
manuscript does not claim novelty in the abstract when the result
has not been presented. It does not use words such as "seminal",
"critical" (unless the word carries technical meaning), "robust"
(except as a technical term in the decision-making sense),
"unprecedented", "powerful", "cutting-edge", "state-of-the-art",
"novel" used more than once, or "first in the literature" without
a specific comparison. The discussion reports what the method
delivers in neutral terms and acknowledges limitations without
euphemism.

## 3. Rules derived from the literature evidence

Every rule in this section is supported by at least one verbatim
example from the nine sampled papers. Page and paragraph locations
are given for every quotation.

### 3.1 Section headers

Standard top-level section sequence is Introduction, Methods, Case
Study (or Application), Results, Discussion, Conclusions. Borgomeo
(2015, p. 5382 contents) uses "1. Introduction, 2. Methods, 3.
Application, 4. Summary and Conclusions". Zaniolo (2024, contents)
uses "1. Introduction, 2. Methods, 3. Case study, 4. Results, 5.
Conclusions and usability". Wheeler (2025) uses "Introduction,
Methods, Results, Discussion, Conclusions" with no section numbers
per the JHE house style. Headers are two to seven words, Title
Case or Sentence case depending on journal house style, and
descriptive. Borgomeo subsections in the Methods read "2.1
Rationale", "2.2 Simulated Annealing", "2.3 Constructing Streamflow
Time Series With Simulated Annealing", "2.4 Model Verification
and Validation". None of the sampled papers uses a section header
that previews a result or makes a claim.

The authoritative section outline for the MOEA-FIND draft lives in
`drafts/manuscript_main_draft.md` and overrides any earlier outlines in
notes or planning documents. Subsection headers must be two to seven
words, descriptive, and Title Case. Section titles such as
"Contributions and the honest scope of the claim" violate Rule 2.6 and
must be replaced with the descriptive form ("Contributions").

### 3.2 Paragraph openings

Body paragraphs open with a complete declarative sentence. The
opus evidence base confirms that across ten paragraphs sampled in
Borgomeo section 2 there is not one paragraph opened with a bold
lead, a sentence fragment, or an interjection. Acceptable opening
patterns include a direct topic sentence (Borgomeo 2015 section
2.3 paragraph 1, "A streamflow time series matching a set of
desired properties is generated by iterative improvement",
p. 5385), a citation-anchored historical claim (Borgomeo 2015
methods paragraph 1, "The first stochastic model for synthetic
streamflow generation was developed by Sudler [1927]", p. 5383),
or an ordinal signpost tying the paragraph to a stepwise procedure
introduced in the previous paragraph (Borgomeo 2015 section 2.2
paragraph 6, "Five elements are required for the algorithm to be
implemented", p. 5386). Paragraph-initial transition words such
as "Thus", "Therefore", "Hence", "However", and "In contrast" are
rare; in the ten sampled Borgomeo paragraphs the count of each is
zero. When a transition is needed it is placed mid-sentence, as
in "Our method is constructed on this same premise; however, we
formulate the streamflow generation problem..." (Borgomeo 2015,
p. 5384).

### 3.3 Paragraph length and structure

Body paragraphs in the sampled predecessors are four to eight
sentences. Methods paragraphs tend toward the upper end because
they interleave display equations. Discussion and Conclusions
paragraphs tend toward the shorter end. A paragraph that is only
one or two sentences is a sign that the content belongs in the
preceding paragraph. A paragraph that exceeds ten sentences is a
sign that a topic shift has been missed and should be split.

### 3.4 Mathematical notation

The three direct predecessors share the same notation conventions.
Scalars are italic Latin or Greek letters. Indices are italic when
they index, upright when they abbreviate a word. Reed (2013,
p. 439) writes the multi-objective problem with an italic boldface
decision vector, "minimize F(x) = [f_1(x), f_2(x), ..., f_M(x)]; x
in Omega; subject to c_m(x) = 0 for all m in E; c_n(x) <= 0 for
all n in I". Borgomeo (2015, Eq. (2), p. 5385) writes the
single-objective simulated annealing objective as a numbered
display equation with italic subscripts and the explicit
summation operator sum_{k=1}^{K}. Wheeler (2025, Eq. (2),
p. 04024056-3) writes the multisite extension with a nested sum
over sites. All three papers number equations sequentially within
the paper, right-justified in parentheses, with no per-section
numbering.

The MOEA-FIND draft will therefore use LaTeX math in display
equations wrapped in `$$...$$` delimiters or the journal's
environment of choice, with sequential numbering. Inline math uses
single-dollar delimiters. Scalars are italicised. Vectors are
bolded italic as in Reed (2013). Sets are upright capitals (Omega,
X, D, S). Operators such as `sum`, `max`, `min`, and `exp` are
rendered as `\sum`, `\max`, `\min`, `\exp`. Absolute value is
`\lvert ... \rvert`. Set membership is `\in`. The notation `(k+1)`
is replaced everywhere by `K+1` with an explicit reminder that `K`
is the number of target drought characteristic dimensions and
`K+1` is the number of objectives Borg MOEA sees. Under no
circumstances will inline code formatting (backticks) appear in
the body text of the manuscript or SI.

### 3.5 Bold and italic usage

Bold in body text is reserved for `Fig.`, `Figure`, `Table`, and
the paper title. No body paragraph opens with a bold lead word
anywhere in the sampled corpus, and the MOEA-FIND draft inherits
this convention. Italic is used for variable names, for foreign
phrases such as `et al.` and `a priori`, for the first use of a
technical term being defined in place, and for titles of cited
works. Italic is not used for general emphasis.

### 3.6 Lists and enumerations

Lists appear in the sampled papers only in three settings.
Borgomeo (2015, section 2.3, p. 5386) uses a numbered list for the
six canonical steps of the simulated annealing loop. Zaniolo
(2024, section 2.3, pp. 10 to 12) uses a bulleted list with italic
labels such as *Drought frequency deviation:* to enumerate the
five components of the weighted-sum objective, each bullet
containing a full sentence and the associated display equation.
Wheeler (2025, Methods, p. 04024056-3) uses a numbered list of
eight algorithmic steps in the simulated annealing loop and a
bulleted list in the Results to enumerate the six Nile subbasins.
Not one of the sampled papers uses a bullet list in the
Introduction, Discussion, or Conclusions. The MOEA-FIND draft
inherits this rule. Where a list is genuinely necessary to
enumerate procedural steps, objective components, or geographic
entities, it is allowed in the Methods only and its items are
either short numbered steps or short bulleted lines with full
sentences.

### 3.7 Em-dashes and hyphenation

Em-dashes (U+2014) do not appear as a primary rhetorical device in
any of the sampled papers. Commas, semicolons, and parentheses
carry the rhetorical load. Hyphens in compound adjectives
("risk-based", "data-driven", "bottom-up") are universal and are
not em-dashes. The MOEA-FIND draft contains zero em-dashes in
body text as verified on 2026-04-14.

### 3.8 Abstract structure

All three direct predecessors use a single-paragraph abstract of
approximately 10 sentences, no bullets, no subsection headers. The
template supported by the evidence is: one motivation sentence
naming the water-system problem, one to two sentences on the
specific limitation in current practice, one sentence announcing
the contribution or proposed method, two to four sentences on the
mechanism, one to two sentences naming the demonstration case, one
to three sentences on the main findings, and an optional closing
sentence on significance or code availability. The opening
sentence of the MOEA-FIND abstract must identify a water resources
problem (for example "Robust evaluation of reservoir operating
policies requires synthetic streamflow ensembles that span the
range of physically plausible drought outcomes") rather than a
mathematical claim.

### 3.9 Figure and table captions

Captions are either short fragments of fewer than 20 words that
name the object of the figure (Borgomeo 2015 Figure 1 caption,
p. 5386, "Flow chart of the simulated annealing synthetic
streamflow generation method") or full sentences describing each
panel (Borgomeo 2015 Figure 10 caption, p. 5395, "Box-plots of the
(a) mean, (b) standard deviation, and (c) interannual lag-1
autocorrelation statistics for 100 simulated sequences using
simulated annealing (SA Hist) and an AR 1 process. The horizontal
black lines show the same statistics for the observed data").
Panel markers are `(a) ... (b) ...` or `Panel a. ... Panel b. ...`
in plain text. Cross-references to earlier figures and tables
("Roman numerals correspond to the map in Fig. 1", Wheeler 2025
Figure 5, p. 04024056-7) are permitted. No caption opens with a
bold lead or a claim.

### 3.10 Paragraph flow and transitions

Transitions rely on lexical repetition of domain terms and
ordinal signposting rather than on paragraph-initial connectives.
Borgomeo section 2 uses "First", "Second", "The third step",
"Finally" to walk through the simulated annealing algorithm and
uses citation-framed openings such as "As noted by Bardossy [1998]"
to introduce motivating references. Zaniolo section 2 uses
"First, ... Second, ... Third, ... Fourth, ... Finally" at
paragraph boundaries. Wheeler Methods uses "The second step is to
generate..." and "The third step is to invoke the simulated
annealing algorithm...". The MOEA-FIND draft will use the same
pattern in its Methods section.

## 4. Domain vocabulary

The following terms are standard in the sampled literature and can
be used in the MOEA-FIND draft without redefinition provided they
are cited correctly on first use. Terms marked with an asterisk are
used differently across the three direct predecessors and must be
explicitly defined on first use in the MOEA-FIND draft.

From the decision-making-under-deep-uncertainty tradition:
robustness (with a reference to the family of robustness metrics,
Hadjimichael 2020 p. 2), satisficing (Simon 1956), regret
(Savage 1951), expected value (Wald 1950), exploratory modelling
(Bankes 1993), scenario discovery, consequential scenarios
(Hadjimichael 2020), deep uncertainty and deeply uncertain factors
(Kwakkel et al. 2010; Lempert 2002; Lempert et al. 2003), Robust
Decision Making and RDM (Lempert et al. 2006, 2010). MORDM,
XLRM, and "many-objective robust decision making" are not
universal across the sampled papers and must be defined and cited
on first use in the MOEA-FIND draft.

From the synthetic streamflow tradition: synthetic hydrology
(Matalas 1967), synthetic streamflow generator, bootstrap,
moving blocks bootstrap, k-nearest-neighbour resampling (Nowak
et al. 2010), Kirsch-Nowak generator. These are treated as
standard named methods with a citation only.

From the drought characterisation tradition: Standardized
Precipitation Index or SPI (McKee et al. 1993), Standardized
Streamflow Index or SSI also known as Standardized Runoff Index
(Zaniolo 2024 section 2.1), Standardized Precipitation
Evapotranspiration Index or SPEI (Vicente-Serrano et al. 2010).
Drought frequency, intensity, duration, and severity are
standard but defined differently by Zaniolo (2024, p. 6, SSI
months) and Wheeler (2025, p. 04024056-8, Q75 years). The
MOEA-FIND draft must state its own definitions explicitly with a
citation to the convention it is following, which is the Zaniolo
SSI convention.

From the evolutionary multi-objective optimisation tradition:
Pareto front, Pareto optimal, Pareto dominated, decision vector,
decision variable, objective vector, feasible set, ε-dominance
(Laumanns et al. 2002), ε-box, ε-progress (Hadka and Reed 2013),
number of function evaluations or NFE (Reed 2013 notation index),
hypervolume indicator (Reed 2013), controllability (Hadka and
Reed 2015 section 2.3). These are defined on first use with a
citation to the primary reference. "Epsilon-box tiling" is not
standard terminology in the sampled papers and should be
introduced with a brief gloss the first time it appears in
section 3.

From the Delaware River Basin and water supply domain:
Flexible Flow Management Program or FFMP (which must be defined
on first use because it does not appear in any of the sampled
papers), Cannonsville, Pepacton, Neversink, Montague flow
target, Trenton flow target, New York City Department of
Environmental Protection.

Terms that do not appear in any sampled paper and therefore
require an explicit in-text definition on first use include
"image of a map", "codimension-1 affine subspace", "bijective
projection", "orthant occupancy", and "feasibility discovery".
Each of these must be introduced with a one-sentence gloss and,
where applicable, a citation to a standard mathematical reference.
Alternatively, they can be replaced with language from the water
resources tradition: for example, "image of a map" can become
"the set of drought characteristic vectors produced by the
generator over all feasible decision vectors", and "codimension-1
affine subspace" can become "a K-dimensional subset of the K+1
objective space defined by the constraint that the sum of all
objectives is a constant determined by the anti-ideal".

## 5. Common anti-patterns to check against every commit

The following specific anti-patterns were observed in the MOEA-FIND
draft and SI draft before 2026-04-14. Each is now disallowed.

5.1. Bold-led topic sentences. Any paragraph opening with
`**Topic.** The...` is disallowed. The content of the topic phrase
either becomes the subject of the opening sentence or is folded
into the preceding paragraph. Before committing an edit, search
the diff for the regex `^\*\*[A-Z][^*]{1,80}\.\*\*` and remove any
match.

5.2. Triple-backtick code blocks. Any text wrapped in triple
backticks is disallowed in the manuscript and SI drafts. Display
equations use LaTeX math delimiters. Pseudocode is rendered as a
numbered list of natural-language steps interleaved with display
equations, following the Borgomeo section 2.3 and Wheeler Methods
templates.

5.3. Inline backtick formatting. Any symbol or variable name
wrapped in single backticks is disallowed in the manuscript and
SI drafts. Variable names are italicised via LaTeX.

5.4. Em-dashes. The character U+2014 is disallowed in the
manuscript and SI drafts. The replacement is usually a comma, a
semicolon, or a parenthetical.

5.5. Undefined technical jargon. Any sentence that introduces a
technical term without a citation, a one-sentence gloss, or a
reference to an earlier definition is disallowed. The word
"image" is the canonical example: if it appears at all, it must
be defined as "the image of the generator map `g`, meaning the
set `\{g(x) : x \in X\}` of all drought characteristic vectors
reachable from the generator", or else replaced with plain
language.

5.6. Declarative section titles and section titles that preview
results. "Contributions and the honest scope of the claim",
"Automated feasible-region discovery", and "The L1 simplex theorem"
are disallowed. Replacements are "Contributions", "Method",
"L1 Manhattan-distance construction".

5.7. Hedging filler adverbs. "Arguably", "essentially",
"effectively", "ultimately", "quite", "rather", "fairly", "very",
"clearly", "obviously", and "of course" are disallowed in body
text unless the word carries specific meaning. "Approximately",
"roughly", and "on the order of" are acceptable when they qualify
a numeric estimate.

5.8. Self-aggrandising adjectives. "Seminal", "unprecedented",
"cutting-edge", "state-of-the-art", "powerful", "robust" (except
as a technical term), and "novel" are disallowed. "First in the
literature" is disallowed unless accompanied by a specific
comparison set and a clear statement of what the comparison
covers.

5.9. Bullet lists in Introduction, Discussion, or Conclusions.
Disallowed per Rule 2.5 and section 3.6 of this guide.

5.10. Unnumbered paragraph headers in body text. Any text of the
form `**Header.** text` or `**Header**: text` inside a paragraph
is disallowed. Sub-topics are introduced with a new paragraph and
a complete opening sentence.

5.11. Direct quotation of other publications. Verbatim quotations
from cited works are disallowed in the main text, the Supporting
Information, the Abstract, the Plain Language Summary, and the Key
Points. Every cited claim is paraphrased in the author's own words
and followed by a parenthetical citation. This rule is strict and
applies even when the original wording is more precise than the
paraphrase. Quotations are permitted only in the working notes under
`scratch/` or `literature/notes/` for author reference, never in a file
destined for submission. The reason is that academic papers in
*Water Resources Research* and its subfield do not use direct
quotation as a rhetorical device except in history-of-science or
controversy contexts; using it in a methods paper signals that the
author could not internalise the cited material enough to restate
it, and reviewers read it that way. Before every commit, search
the diff for the regex `"[A-Z][^"]{20,}"` and for smart-quote
variants, and remove or paraphrase any match that is not a
dataset name or a proper noun.

5.12. Verbose citation formatting. In-line citations use the
compact form `(Author, Year)` or, for a single author, `Author
(Year) showed...`. Journal names, volumes, page numbers, and
section numbers do not appear in the in-line citation. The full
reference lives in the bibliography only. Citing a specific page
or section is reserved for the rare case where the manuscript
disputes or refines a specific sentence in the cited work, and
even then the form is `(Author, Year, p. N)` with no journal name.
Forms like `Borgomeo, Farmer, and Hall (2015, *Water Resources
Research*)` or `Bonham et al. (2025, section 2.1.2)` are
disallowed in body text.

5.13. Quoted numeric results from cited publications. Transcribing
a specific metric value, coefficient, sample size, or experimental
count from a cited paper into the body text is disallowed. The
claim is rephrased as a qualitative takeaway supporting the
manuscript's own argument. For example, a statement such as
"Bonham et al. (2024) show that input-space space-filling metrics
predict robustness ranking preservation with coefficient of
determination between 0.77 and 0.91" is replaced by "Bonham et al.
(2024) report that input-space space-filling metrics are strong
predictors of robustness ranking preservation, but conduct the
entire analysis in an input coordinate system rather than in the
hazard characteristic coordinates of interest to the planner".
The reason is that transcribed numbers tie the manuscript's
argument to the specific experimental conditions under which the
original number was measured, limit the generality of the point,
and signal a compilation-style reading of the literature.
Numeric results quoted from the manuscript's own experiments are
obviously still required.

5.14. Forward references that break narrative sequencing. A term
introduced in section `N` must not be used in any earlier section
without either (a) a parenthetical one-sentence gloss at first use
or (b) a bracketed forward reference such as "(defined in section
M)" that is explicit about its purpose. Concrete failure modes to
avoid: naming a case study location in a methods subsection before
the case study has been introduced, naming a metric in the
problem statement before it has been defined, naming a figure
number in prose before that figure is introduced by an immediately
adjacent sentence, naming a mathematical object (Manhattan
distance, epsilon-box, projection pi) in a subsection that
precedes its definition. The reason is that reviewers read papers
linearly and every forward reference forces a mental context
switch that the author is expected to eliminate through ordering.
Before every commit, re-read each section in order and verify that
every proper noun, symbol, and method name has been introduced in
the current or a prior section.

5.15. Informal or quippy sentences. Short authorial sentences that
function as rhetorical flourishes rather than as content are
disallowed. Examples of the failure mode: "This paper is concerned
with the construction of that ensemble", "MOEA-FIND is a direct
response to that open question", "This paper presents X", "We now
turn to Y", "This is the central contribution". Replacements are
either deletion (when the sentence adds no information) or merger
into a preceding content-bearing sentence. Short sentences that
state a specific result or a specific definition are not informal
and are permitted. The distinction is whether the sentence carries
new information or only rhetorical weight. Transitional sentences
are a specific subset of this anti-pattern and should be removed
or folded into adjacent paragraphs.

5.16. Informal group designations for research communities.
Phrases that name research groups informally (e.g., "the Reed
group", "the Reed-group literature", "the Kasprzyk group") are
disallowed in body text. Replacement is either citation of the
specific papers the phrase is pointing at or a neutral descriptor
of the research tradition (e.g., "the many-objective robust
decision making literature", "exploratory modelling under deep
uncertainty", "the search-based synthetic streamflow generation
tradition"). The reason is that group-labelling is inside-baseball
language that signals a subfield-insider register, which a
*Water Resources Research* reviewer outside the group will read
as parochial.

5.17. Thin Introduction literature discussion. The Introduction
of a methods paper in this subfield establishes the gap by naming
specific prior works, summarising each in a sentence or two of the
author's own words, and tying each summary to the specific
limitation the manuscript addresses. A one-paragraph mention of a
predecessor that names the paper and its result but does not
explain the methodological choice the predecessor made is too
thin. A two or three sentence paraphrase per load-bearing
predecessor is the minimum; a richer treatment citing the specific
methodological lineage (simulated annealing, Kirsch-Nowak
bootstrap, conditioned Latin hypercube, exploratory modelling,
scenario discovery) is the norm. The author's own method is
justified by explicit contrast with the paraphrased predecessor,
not by novelty claims.

## 6. Pre-commit checklist

Before any edit to `drafts/manuscript_main_draft.md` or
`drafts/supporting_info_draft.md` is finalised, the author confirms
that the edit satisfies every item in this list.

The edit contains no bold-led topic sentences matching
`**[A-Z][a-z].*\.**`. The edit contains no triple-backtick fences
and no inline backtick-quoted variable names. The edit contains no
em-dashes (U+2014). Every new technical term is accompanied by a
citation or a one-sentence gloss on first use. Every new paragraph
opens with a complete declarative sentence. No new section or
subsection header previews a result, makes a claim, or asks a
question. The notation uses `K` and `K+1` consistently and
italicises variables. Lists do not appear outside the Methods, and
where they appear they enumerate procedural steps, objective
components, or named entities. The edit does not introduce any of
the disallowed adjectives listed in section 5.7 or 5.8. Cross
references to sections, figures, tables, and equations use the
forms `section N`, `Figure N`, `Table N`, `Eq. (N)` and are
accurate.

The checklist is run before every commit of the drafted files.
Failures are fixed in place, not noted as followups.

## 7. Maintenance

This guide is updated when the user issues a new style rule, when
a new class of violation is observed, or when a new literature
precedent supersedes a previous one. The current version is
2026-04-14. Updates are dated and the superseded text is preserved
as a dated entry below.
