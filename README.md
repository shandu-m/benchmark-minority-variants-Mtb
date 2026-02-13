# Benchmarking within-sample minority variant detection with short-read sequencing in *M. tuberculosis*

Shandukani Mulaudzi<sup>1</sup>, Sanjana Kulkarni<sup>1</sup>, Maximillian G. Marin<sup>1,2</sup>, Maha Farhat<sup>1,3</sup>

<sup>1</sup>Department of Biomedical Informatics, Harvard University, Boston, MA 02115, USA
<sup>2</sup>Department of Data Science, Dana-Farber Cancer Institute, Boston, MA 02215, USA
<sup>3</sup>Division of Pulmonary and Critical Care, Department of Medicine, Massachusetts General Hospital, Boston, MA 02114, USA

This is the project repo describing all sequence data simulation, sequence processing, variant calling and data analysis associated with our publication.

- `data_simulation`: Describes our data simulation process.
- `wgs_processing`: Describes how we processed the simulated sequencing data.
- `variant_calling`: Describes how the benchmarked variant callers were run on the simulated strains.
- `variant_filtering`: Describes how we used an error model and allele fraction adjustment to filter FreeBayes SNVs and INDEL variants respectively. The error model is provided here to be used on FreeBayes variant calls.
- `analysis`: All benchmarking analysis.
