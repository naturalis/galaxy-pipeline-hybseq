# galaxy-pipeline-hybseq

## This is a tagged version of an attempt to implement Hybpiper/Galaxy in separate steps, abandoned by @Eremus007

HybSeq is a data generation technique that constitutes a hybrid between shotgun high-throughput sequencing
and amplicon sequencing. The technique results in markers suitable for phylogenomics (because they are 
single copy, universally present, and evolve at the optimal rate) by targeted sequencing using predefined
bait kits, such as [developed by Kew](https://pubmed.ncbi.nlm.nih.gov/31477409/) and others. As such, the
sequencing results in multiplexed FASTQ data that both needs to be demultiplexed back to the input species
as well as sorted by the various baits. The ideal end result is a set of as many multiple sequence alignments
(MSAs) as there were baits, each containing as many rows as there were species, as the input for a multi-marker
phylogenetic analysis yielding a robust tree. In this repository, the artefacts are developed to achieve
this ideal end state. The repository is inspired by the following prior art:

- [Nikolov et al., 2019](https://doi.org/10.1111/nph.15732) - a seminal paper describing how to resolve 
  the phylogenetic backbone of the Brassicaceae using HybSeq. The paper describes an informal reference
  pipeline.
- [brassicaceae-hybseq-pipeline](https://github.com/naturalis/brassicaceae-hybseq-pipeline) - an 
  implementation of the Nikolov et al. pipeline using snakemake. The implementation is not quite complete
  but very close, lacking the last few steps of MSA file conversion, concatenation, and phylogenetics.
- [galaxy-tool-hybseq](https://github.com/naturalis/galaxy-tool-hybseq) - a first attempt to port the
  snakemake pipeline to Galaxy. Many of the tools used by Nikolov et al. are already available in Galaxy's
  package management system (the 'tool shed') but some are missing. This repo attempts to wrap those missing
  tools, especially YASRA.
- [HybPiper](https://bsapubs.onlinelibrary.wiley.com/doi/10.3732/apps.1600016) - a more current HybSeq
  pipeline. It is very possible that this is preferable over the Nikolov et al. pipeline. This needs to be
  investigated.
- [Naturalis Galaxy Portal](https://galaxy.naturalis.nl/), a local install of the Galaxy workflow management
  system. This is where the wrappers for the different pipeline steps will be deployed eventually. The 
  portal is administered by Dick Groenenberg.
  
# Goals

- Ability to stage HybSeq data (from SRA?) into the Galaxy portal
- Implementation of all steps of a HybSeq pipeline as Galaxy tools, where needed as wrappers and scripts
  deposited in this repository
- A Galaxy [workflow](https://galaxyproject.org/learn/advanced-workflow/) that chains all tools with 
  sensible defaults
- A plausible phylogenetic tree result  
