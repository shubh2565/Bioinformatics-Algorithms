# Gotoh Algorithm

The algorithm by Osamu Gotoh (1982) computes the optimal global alignment of two sequences when using an affine gap scoring. Here, the scoring of a long consecutive gap (insertion/deletion) is favored over a collection of small gaps with the same combined length. This incorporates the assumption that a single large insertion/deletion event is biologically more likely to happen compared to many small insertions/deletions. While sophisticated gap scoring models can be applied in the generic algorithm by Waterman-Smith-Beyer (1976), affine gap scoring used in Gotoh's algorithm enables a reasonable gap model with reduced runtime.

Source: [Freiburg RNA Tools](http://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Gotoh)
Original Paper: [ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/0022283682903989?via%3Dihub)
