# GENMOTIF

This is the source code for the paper:

J. Serra, A. Matic, J.L. Arcos, & A. Karatzoglou, "A genetic algorithm to discover flexible motifs with support", arXiv:1511.04986, 2015.

Abstract: Finding repeated patterns or motifs in a time series is an important unsupervised task that has still a number of open issues, starting by the definition of motif. In this paper, we revise the notion of motif support, characterizing it as the number of patterns or repetitions that define a motif. We then propose GENMOTIF, a genetic algorithm to discover motifs with support which, at the same time, is flexible enough to accommodate other motif specifications and task characteristics. GENMOTIF is an anytime algorithm that easily adapts to many situations: searching in a range of segment lengths, applying uniform scaling, dealing with multiple dimensions, using different similarity and grouping criteria, etc. GENMOTIF is also parameter-friendly: it has only two intuitive parameters which, if set within reasonable bounds, do not substantially affect its performance. We demonstrate the value of our approach in a number of synthetic and real-world settings, considering traffic volume measurements, accelerometer signals, and telephone call records.

http://arxiv.org/abs/1511.04986


To compile, you'll need to run "python build.py folder" and have gcc installed. There is also a python wrapper for running the compiled binary in the bin/ folder.

If using this code, parts of it, or developments from it, please cite the above reference. We do not provide any support or assistance for the supplied code nor we offer any other compilation/variant of it. In addition, we assume no responsibility regarding the provided code.

