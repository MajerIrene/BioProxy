# BioProxy #

Escherichia coli is a gram negative bacterium and it is the most studied bacterial model organism for many reasons: its ability to grow fast on cheap media; the fact that it can be genetically manipulated and modified much easier with respect to eukaryotic organisms or even other bacteria. When an organism is as studied as E. coli the amount of data available about it can be massive. For this reason we would like to perform a bioinformatic analysis to extract knowledge from one of the many RNA-seq dataset available for E. coli.

This project aims to establish a proxy for determining the activity of transcription factors based on the expression levels of the genes regulated by those factors. In prokaryotes genes with related functions are encoded together within the same genome block, which is called operon, and they are usually transcribed together in a so called polycistronic transcript, because they are under the control of the same promoter.

The null hypothesis posits that the abundance of transcripts serves as a reliable proxy for activity of the regulator. However, this assumption may not hold for all transcription factors: some of them are active only when phosphorylated, meaning that their transcripts may be constitutively expressed, while their activity depends on phosphorylation, correlating instead with the activity of the corresponding kinase.

We will also focus on developing a general model to relate transcription factor activity to transcript abundance. Significant deviations from this model could indicate transcription factors whose activity is not accurately reflected by transcript levels. One challenge of this approach is the fact that genes of the same operon may be highly co-expressed together, which is a problem for many models like linear regression, so we will try to address it during our analysis.

