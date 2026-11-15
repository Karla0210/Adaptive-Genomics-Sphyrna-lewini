
# Quality filtering

Scripts used for adapter removal, read trimming, and quality assessment of sequencing reads from *Sphyrna lewini*.

## Trimmomatic

Adapter sequences were removed using Trimmomatic v0.39 with Nextera adapter settings. Reads were truncated to 120 bp (`CROP:120`), and reads shorter than 120 bp were discarded (`MINLEN:120`).

## Quality assessment: FastQC and MultiQC

Sequence quality was assessed using FastQC, and the resulting reports were summarized using MultiQC.

The mean Phred quality scores remained high across read positions, generally ranging from Q36 to Q40, indicating high sequencing quality after trimming.

### Mean quality scores

<img width="1502" height="622" alt="image" src="https://github.com/user-attachments/assets/70654e27-1f51-4e32-b45a-0ebad473cb73" />


*Figure 1. Mean Phred quality scores across sequencing reads, visualized using MultiQC.*

