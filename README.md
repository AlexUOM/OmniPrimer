# OmniPrimer: complete automation of primer design for the NHS

## Overview
This repository hosts the code of OmniPrimer, a bioinformatic automation tool aimed to optimise the long and time-consuming clinical procedure to design primers for the [Sheffield Diagnostic Genetics Service](https://www.sheffieldchildrens.nhs.uk/sdgs/) (part of the Sheffield Children’s NHS Foundation Trust). OmniPrimer effectively alleviates the burden on the NHS, saving hundreds of hours of clinicians' work. This translates into significant cost savings for the NHS, amounting to tens of thousands of pounds, and enhances patient outcomes by enabling clinicians to focus their efforts and attention directly on patient care.

## Table of Contents
- [Introduction](#introduction)
- [Repository Structure](#repository-structure)
- [Usage](#usage)
- [Dependencies](#dependencies)
- [About](#about)

## Introduction
Accurate and efficient primer design is essential in molecular diagnostics, particularly in diagnosing rare genetic diseases where precision is paramount. To ensure accurate diagnosis, primers must bind to a genomic region free from single nucleotide polimorphisms (SNP). This ensures that the amplification of the exon of interst takes place, improving the realibity of the genetic testing and subsequent diagnosis. However, generating clinical-grade diagnosis primers is currently a long and time-consuming process that clinicians must perform manually. 

OmniPrimer was born to offer a comprehensive and automated solution to clinicians, integrating the existing clinical methodologies to generate primers with cutting-edge sorting algorithms. The tool was single-handedly developed in Python as part of a work placement for the Sheffield Children’s NHS Foundation Trust.

## Repository Structure
The repository is organized as follows:
- **`src/`**: Contains the source code of OmniPrimer and related utils.
- **`data/`**: Includes sample input data and output from OmniPrimer.
- **`docs/`**: Documentation and additional resources related to OmniPrimer.
- **`README.md`**: This README file providing an overview of the repository and instructions for usage.

## Usage
To utilize the bioinformatic pipeline, follow these steps:
1. Clone the repository to your local machine.
2. Install the necessary dependencies (see Dependencies section).
3. Prior to running OmniPrimer, ensure the following input file is available:
   - `.txt` file containing the full intronic and exonic sequence of the gene of interest. This should be downloaded from [Ensembl](http://www.ensembl.org/index.html).
5. Navigate to the `src/` directory.
6. Run the main OmniPrimer script with appropriate input parameters as described below to initiate the primer design process:
   - Set your desired download directory path using `profile.set_preference("browser.download.dir")` to save primer pairs.
   - Set the Ensembl page of the gene of interest using `browser_primers.get()`. Simply replace the placeholder link with yours.
7. Access the generated primers and analysis reports in the designated output directory.

## Dependencies
Ensure you have the following dependencies installed:
- Python (version X.X)
- BioPython (version X.X)
- [Additional dependencies, if any]

## About
This bioinformatic pipeline was developed by [Your Name] and [Contributor Name] as part of [Project/Research Institution]. It is released under the [License Type] license. For any questions, feedback, or issues, please contact [Your Email] or open an issue on GitHub.

**Disclaimer**: This pipeline is intended for research purposes only and should not be used for clinical diagnosis without proper validation and regulatory approval.

---
**Note**: Replace placeholders (e.g., [Your Name], [License Type], etc.) with appropriate information specific to your project/
