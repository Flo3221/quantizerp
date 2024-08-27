# QuantizerP

## Description

QuantizerP is a project developed to illustrate the concepts presented in the paper "Simultaneous compensation of input delay and state/input quantization
	for linear systems via switched predictor feedback" ([SCL](https://www.sciencedirect.com/science/article/pii/S0167691124002007?ref=pdf_download&fr=RR-2&rr=8b9b8fc31926eec4) ). 
 
 This project simulates a linear time-delay system. For detailed equations and parameters, 
 refer to the ([Mathematical Details](https://tucgr-my.sharepoint.com/:b:/g/personal/fkoudohode_tuc_gr/EQMZt_JoHTxClA2cSibvOacBsN0FpqYPkZIoLCFf_xnM8w?e=PFkfER)).

## Requirements

To run this project, you will need:

- MATLAB R2023b or later.

## Installation

Follow these steps to set up the project:

1. Download the project files from [QuantizerP GitHub Repository](https://github.com/flo3221/quantizerp).
2. Extract the contents to a directory of your choice.
3. Open MATLAB and navigate to the project directory using the `cd` command:
    ```Matlab
    cd /path/to/quantizerp
    ```

## Usage

To use QuantizerP, follow these steps:

1. Open MATLAB and ensure you are in the project directory.
2. Run the main script or function:
    ```Matlab
    Quantized_predictor_feedback.m
    ```
3. Ensure that the `private` folder is in the same directory. This folder is used to solve initial-boundary value problems for first-order systems of hyperbolic partial differential equations (PDEs).

### Functions

QuantizerP includes the following key functions:

- `hpde.m` and `setup.m`: These functions solve the transport PDE described.
- `mu`: Implements the switching parameter $\mu(t)$.
- `quantizer`: Implements the quantizer function.

## Examples

Refer to the following scripts for examples of how to use QuantizerP:

- `Quantized_predictor_feedback.m`
- `fixed_mu_quantized_predictor_feedback.m` (for the fixed switching parameter case)

## Contributing

To contribute to QuantizerP, please follow these steps:

1. Fork the repository on GitHub.
2. Create a new branch for your feature or fix.
3. Make your changes and commit them.
4. Submit a pull request with a detailed description of your changes.

## License

This project is licensed under the CC BY-NC-ND  ([`LICENSE`](https://creativecommons.org/licenses/by-nc-nd/4.0/)).

## Contact

For questions or feedback, please contact [fkoudohode@tuc.gr](mailto:fkoudohode@tuc.gr).

## Cite this work
```
@article{fkoudohode2024,
title = {Simultaneous compensation of input delay and state/input quantization for linear systems via switched predictor feedback},
journal = {Systems & Control Letters},
volume = {192},
pages = {105912},
year = {2024},
issn = {0167-6911},
doi = {https://doi.org/10.1016/j.sysconle.2024.105912},
url = {https://www.sciencedirect.com/science/article/pii/S0167691124002007},
author = {F. Koudohode and N. Bekiaris-Liberis},
}
```
