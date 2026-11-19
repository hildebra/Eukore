"""The versions of the build's steps (since 2026-11-19): what a later run goes by to take a step's kept result or to
make it again (docs/building.md#when-a-step-runs-again-the-steps-versions).

A step that keeps its result for later runs (a genome's calls in WORK/done, a kept search part, a parsed RefSeq
release, the published reads' annotation, ...) has a version here, MAJOR.MINOR. Its kept result holds the MAJOR that
made it (in a field, its key or its file's name: each step's note below says where); a later run takes it where that
is the step's MAJOR and makes it again where not.

- MINOR is raised for a change to the step's code that leaves its results as they were: a crash or its speed fixed,
  a message, a refactor, a new option's default that does what the code did. Nothing is run again.
- MAJOR is raised (MINOR back to 0) for a change that changes its results: another call, rule, threshold, column or
  file. Every kept result of the step is made again by the next run (a genome's: unless --reversion keeps it), so
  that a package never mixes results of two MAJORs.
Each step's history is kept beside it: the version, the date, what changed.

The code's own text is never taken for this: until 2026-11-19 the process stage's WORK/process.nothing_left held a
digest of every module's source, so that any edit, a comment's too, sent a SLURM array's tasks through the plan again.
What a step's result is made from is still in its key: its inputs' files, the options, the models, and where a step's
key holds them the values of its rules (organelle.rules_text, markers.recipe, rrna's join limits, the thresholds in
organelle.search_key and in qc's parts' keys): those change only with an edit that changes the results.

A version's MAJOR 1 is that of the results a step made before it had a version: a step versioned since 2026-11-19
keeps the keys it had at MAJOR 1 (tag), so that its kept results stay good.
"""

import hashlib

STEPS = {
    # A genome's run: its BUSCOs' calls, copies and features, its rRNA, its organelle molecules (genomes.run and what it
    # calls: loci, odb, placement, align, rrna's extraction). WORK/done/<accession>'s done.json `genome_record` (and the
    # package's genomes.tsv); a done genome of an older MAJOR is run again unless --reversion stamps it with this one.
    # 1: the runs until 2026-11-09.
    # 2 (2026-11-09): the coding sequence checked against the protein confirmed and recovered (loci.cds_state,
    #    cds_check); a region's clean alignment first; a Duplicated locus' other complete copies kept whole
    #    (genomes.ALTERNATES); an annotated model with a stop inside, or a fusion, not taken.
    # 3 (2026-11-14): a copy split across contigs joined (loci.join_split, `model` joined), the first pass' miniprot at
    #    the join's --outs and its alignments of loci.OUTS taken from those.
    # 4 (2026-11-15): a call complete against its query protein's own HMM coverage (odb.relative): a lineage's whole
    #    genes no longer Fragmented where its proteins lack a stretch of the profile.
    # 4.0 (2026-11-19): versioned here (genomes.GENOME_RECORD).
    "genome": "4.0",

    # Which genomes and passes a process run takes (eukore_build.plan_of, genomes_to_run, the passes after the genomes):
    # WORK/process.nothing_left, a --genomes_only task's note that it found none left, holds every step's MAJOR
    # (digest; eukore_build.process_fingerprint). MAJOR raised where a change makes a process run take a genome or pass
    # it did not (another rule of what is run again), so that an array's later tasks plan again rather than end at once.
    # 1.0 (2026-11-19): versioned here (until then the note held a digest of every module's text).
    "process_plan": "1.0",

    # A locus' proteins' rows on its profile, from hmmalign (align.align_locus: WORK/align/<locus>.<the HMM's MD5><tag>,
    # a line a distinct protein). MAJOR raised where a row's reading from hmmalign's output changes.
    # 1.0 (2026-11-19): versioned here.
    "protein_rows": "1.0",

    # The checks across species' kept vsearch searches of a locus' copies (qc.searched_loci; WORK/qc/): in the parts'
    # key (qc.Params.text), with the options and files. MAJOR raised where the search's options or files change.
    # 1.0 (2026-11-19): versioned here (qc.SEARCH_RECORD).
    "qc_search": "1.0",
    # The 18S' kept cmsearch parts of the SSU models (qc.CmParams, regions.ssu_models): in their key. MAJOR raised where
    # util.CMSEARCH_OPTIONS or the part's table change.
    # 1.0 (2026-11-19): versioned here (qc.CM_RECORD).
    "qc_cmsearch": "1.0",
    # The 18S' kept cmalign parts of the variable regions (qc.CmalignParams): in their key. MAJOR raised where
    # regions.align_file's lines change.
    # 1.0 (2026-11-19): versioned here (qc.CMALIGN_RECORD).
    "qc_cmalign": "1.0",
    # The 18S' kept window-screen parts (qc.MapParams): in their key. MAJOR raised where chimeras.map_file's lines
    # change.
    # 1.0 (2026-11-19): versioned here (qc.MAP_RECORD).
    "qc_window_screen": "1.0",
    # The 18S' kept UCHIME parts (qc.UchimeParams): in their key. MAJOR raised where chimeras.uchime_file's lines
    # change.
    # 2 (2026-11-15): the whole database, a query's own record left out by --self.
    # 2.0 (2026-11-19): versioned here (qc.UCHIME_RECORD).
    "qc_uchime": "2.0",

    # The 18S' SSU model calls from the models' cmsearch hits (regions.ssu_models: WORK/cmalign/models.<key><tag>, the
    # key the models' file's, SSU_MODELS' and SSU_MARGIN's). MAJOR raised where the calls' rules change (with
    # qc_cmsearch's where cmsearch's own hits change).
    # 1.0 (2026-11-19): versioned here.
    "ssu_models": "1.0",
    # The 18S' variable regions placed by cmalign (regions.run: WORK/cmalign/ssu.<key><tag>, the key the model's, the
    # regions table's, the reference's, SLACK's and MAX_SSU's; the 18S cmalign failed on: unaligned.<key><tag>). MAJOR
    # raised where the regions' placing changes (with qc_cmalign's where cmalign's own parts change).
    # 1.0 (2026-11-19): versioned here.
    "ssu_regions": "1.0",

    # A read's kept audit of its ends and flanks (audit.audit_text: the field's first value, in the regions caches).
    # MAJOR raised where the ends' or flanks' rules change.
    # 1.0 (2026-11-19): versioned here (audit.AUDIT_RECORD).
    "rrna_audit": "1.0",
    # The saved extractions' 18S records of joined pieces extracted again (rrna.resplit: WORK/rrna/joined.<key>, the key
    # its limits' and the MAJOR's). MAJOR raised where joined's rule changes again.
    # 1.0 (2026-11-19): versioned here (rrna.JOIN_RECORD).
    "rrna_join": "1.0",
    # The published rRNA reads' regions, kept in WORK/published/regions.<MAJOR>.tsv (eukore_build.PUBLISHED_REGIONS)
    # and their batches (published_batches: named by it).
    # 2 (2026-11-06): ITSx's guessed 5.8S borders no longer taken.
    # 3 (2026-11-06): regions on the record, not the unit, and the read's strand.
    # 4 (2026-11-14): an 18S read's 18S alone, each read's audit of its ends and flanks.
    # 4.0 (2026-11-19): versioned here.
    "published_regions": "4.0",

    # Eukore's own marker set (SOURCES/eukore_markers/, markers.ready: its state's `recipe`, a digest of
    # markers.RECIPE, its MAJOR, and of the build's rules' constants: markers.recipe). MAJOR raised with a change of its
    # rules that none of those constants shows; a constant changed builds it again by itself.
    # 1.0 (2026-11-19): versioned here (markers.RECIPE).
    "markers": "1.0",

    # The organelle genes' kept results: RefSeq's parsed release (organelle.refseq_genes: WORK/organelle/refseq_<key>,
    # with organelle.rules_text) and each assembly's search (organelle.search_key, its #searched line's inputs). To be
    # raised where record_genes', the searches' or the models' logic changes what they hold.
    # 1: until 2026-11-10 (never written: a cache outlived such a change).
    # 2 (2026-11-10): /pseudo, the longest feature, /codon_start, /transl_table, the spacer across the origin.
    # 3 (2026-11-10), on two branches at once: /transl_except on one; on the other a kinetoplastid's 12S its LSU
    #    (rrna_gene) and miniprot's model in place of nhmmer's single hit that is not the gene's coding sequence
    #    (with_miniprot).
    # 4: the two together, so that a cache either branch wrote is made again.
    # 5 (2026-11-15, the audit's O1 and O2): a gene read along its molecule spans its first base to its last
    #    (turns_of), a candidate contig's genes kept only beside a protein-coding gene of its organelle.
    # 5.0 (2026-11-19): versioned here (organelle.CACHE_VERSION).
    "organelle": "5.0",
    # An assembly's search of its candidate contigs (organelle.candidate_molecules; `candidates<MAJOR>` in its
    # #searched line's models, organelle.search_models).
    # 1.0 (2026-11-19): versioned here (organelle.CANDIDATE_RECORD).
    "organelle_candidates": "1.0",
    # One organelle to an assembly's molecule in a species that can hold a plastid (organelle.ONE_ORGANELLE:
    # `one_organelle<MAJOR>` in its #searched line's models; a fungus' or an animal's search is not changed by it).
    # 1 (2026-11-19): a titled molecule no candidate of the other organelle, a mitochondrion's candidate needs a gene
    #    of organelle.MITO_OWN.
    # 1.0 (2026-11-19): versioned here.
    "organelle_one": "1.0",
    # The MIDORI2 records the reference screen screens (organelle.midori2_choice; in the screen's key,
    # eukore_build.reference_screen).
    # 1.0 (2026-11-19): versioned here (organelle.CHOICE_RECORD).
    "midori2_choice": "1.0",
    # The reference screen's kept result (WORK/reference_screen/<key>, eukore_build.reference_screen: the key with its
    # MAJOR; the package's manifest names it). MAJOR raised where the screen's rules change.
    # 1.0 (2026-11-19): versioned here (refscreen.VERSION).
    "reference_screen": "1.0",
}


def version(step):
    """A step's version, MAJOR.MINOR."""
    return STEPS[step]


def major(step):
    """A step's MAJOR: what its kept results hold, and are made again by where it is not this one."""
    return int(STEPS[step].partition(".")[0])


def tag(step, prefix="."):
    """A step's MAJOR as a key's part where it had none before 2026-11-19: none at MAJOR 1 (its kept results of before
    stay good), else prefix and the MAJOR (`.v2`)."""
    n = major(step)
    return "" if n == 1 else f"{prefix}v{n}"


def majors():
    """{step: its MAJOR} of every step."""
    return {step: major(step) for step in STEPS}


def digest():
    """Every step's MAJOR as a short digest: what a result that depends on all of them is kept by
    (eukore_build.process_fingerprint). A MINOR raised leaves it as it was."""
    return hashlib.sha1(repr(sorted(majors().items())).encode()).hexdigest()[:12]


def text():
    """Every step's version as one line (the build's banner, its logs)."""
    return ", ".join(f"{step} {v}" for step, v in STEPS.items())
