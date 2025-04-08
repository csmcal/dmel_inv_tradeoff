# Drosophila Inversion Tradeoff Analysis

This repository contains the analysis pipeline and tools used in McAllester & Pool 2025 to study inversion polymorphisms in *Drosophila melanogaster*. The code handles:

- Read alignment and quality control
- Inversion frequency calling
- Statistical analysis of frequency data

The SAIsim simulator used in this study is housed in a separate repository at [https://github.com/csmcal/SAIsim](https://github.com/csmcal/SAIsim).

## Repository Structure

- `analysis/`: Main analysis scripts and pipelines
  - `Alignment/`: Read alignment and processing tools
  - `Significance/`: Statistical analysis scripts
  - `Graphing/`: Data plotting scripts

## Citation
If you have any reason to reference this repository, please cite:

> McAllester, C. S., & Pool, J. E. (2025). The potential of inversions to accumulate balanced sexual antagonism is supported by simulations and Drosophila experiments. *eLife*.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.