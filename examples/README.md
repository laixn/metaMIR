# About
These are sample gene list files for use with metaMIR in various modes of operation


### example1.txt
A short analysis of genes from the [WNT](https://en.wikipedia.org/wiki/Wnt_signaling_pathway) signaling pathway, in this standard mode, the script will search for miRNA candidates that will target all or some subset of the listed genes.


### example2.txt
The 'core' analysis mode. Here, some of the supplied genes are prefixed with asterisks. The analysis is performed such that any combination list returned in the analysis must contain at least one of these genes. There is no imposed restriction on how many of the input genes can be marked in this way, however the more that are required to appear in the output, the more restricted the analysis will be and potentially may not return any results.

The analysis consists of genes from the [Hippo/YAP](https://en.wikipedia.org/wiki/Hippo_signaling_pathway) signaling pathway and transcription factors associated with neural plate border formation during the process of [neural crest development](https://en.wikipedia.org/wiki/Neural_crest). The analysis is set to ensure that Hippo/YAP genes must be included in the results returned.


### example3.txt
In this example, the script is run in differential analysis mode. Here, genes prefixed with a minus sign (-) in the analysis are not to be targeted by any returned miRNA. This type of analysis may be more suitable when searching for miRNA candidates that may target groups of genes of interest when specificity or off-target effects are an issue.

In this example, the bulk of genes in the default or 'forward' analysis are drawn from the transforming growth factor beta signaling pathway, or [TGFb signaling](), a pathway that can be overactive in some cancers. The 'negative' genes are example tumour suppressors. The hypothetical analysis here would search for potential therapeutic miRNAs which would decrease expression of TGFb components, while preferentially *not* targeting tumour suppressors.


### example4.txt
The final example contains both the core and differential analyses, where separate genes are marked with asterisks to ensure only results are returned that contain at least one, and other genes are prefixed with minus signs to avoid their targeting according to predictions.

The 'core' genes are from the Hippo signaling pathway, the general 'forward' analysis genes are from the TGFb signaling pathway, and the genes not to be targeted are again tumour suppressors.
