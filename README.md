<br />
<div align="center">
  <a href="https://github.com/cfjimenezv07/HDFTS_Mortality">
    <img src="Kaustlogo.png" alt="KAUST Logo" height="150">
    <img src="MQ.png" alt="MQ Logo" height="150">
  </a>

<h3 align="center">Forecasting high-dimensional functional time series: Application to sub-national age-specific mortality</h3>
</div>

## Abstract
<p align="justify">
We study the modeling and forecasting of high-dimensional functional time series (HDFTS), which can be cross-sectionally correlated and temporally dependent. We introduce a decomposition of the HDFTS into two distinct components: a deterministic component and a residual component that varies over time. The decomposition is derived through the estimation of two-way functional analysis of variance. A functional time series forecasting method, based on functional principal component analysis, is implemented to produce forecasts for the residual component. By combining the forecasts of the residual component with the deterministic component, we obtain forecast curves for multiple populations. We apply the model to age- and sex-specific mortality rates in the United States, France, and Japan, in which there are 51 states, 95 departments, and 47 prefectures, respectively. The proposed method is capable of delivering more accurate point and interval forecasts in forecasting multi-population mortality than several benchmark methods considered.
</p>

### Main Results
The R script files in the `R Code` folder should be used in the following order:

#### 1. Core Methodology & Estimation (`Rcodes_paper/`)
* **FM_ANOVA.R**: Computes the Functional Mean Analysis of Variance (FM) decomposition on FTS to isolate deterministic and time-varying components.
* **FMP_ANOVA.R**: Computes the Functional Multi-way Panel Analysis of Variance based on the Functional Median Polish (FMP) framework to decompose FTS into deterministic and time-varying components.
* **MITS_class.R**: Defines custom object structures and classes to efficiently manipulate and organize multivariate panel functional time series streams.
* **aux_HNT.R**: Serves as the primary helper script housing backend functions and numerical routines for the core high-dimensional estimation framework.

#### 2. Forecasting Engine & Evaluation (`Rcodes_paper/`)
* **method.FPE.R**: Determines the optimal number of functional principal components to retain using the Functional Prediction Error criterion.
* **forecast_Arima.R**: Fits ARIMA models to the estimated functional principal component scores to project them across future horizons.
* **nonparametric_fof_regression.R**: Performs function-on-function non-parametric regressions to link historical panel patterns with future dynamics.
* **New_Point_forecast.R**: Coordinates the primary forecasting engine, computing the point forecasts (PFE) for the functional residuals obtained from FMP-ANOVA and FM-ANOVA.
* **Naive_point_forecasts.R**: Computes the naive point forecasts based on the assumption of cross-sectional independence.
* **sieve_code.R**: Implements the functional sieve bootstrap framework by Paparoditis and Han Lin Shang (JASA 2023) utilized to build simultaneous, uniform prediction bands.
* **Compute_coverages.R**: Computes the empirical coverage probabilities as well as the interval scores (IFE) for the interval forecasts.
* **Independence.R**: Tests the baseline structural assumption of cross-sectional independence among sub-national functional series.
* **HNT23_GSY19.R**: Executes specialized comparative modeling pipelines using recent methodological frameworks from relevant high-dimensional literature.
* **forecast_errors.R**: Computes structural error metrics like Point Forecast Error (PFE) and Interval Forecast Error (IFE) to quantify the precision of sub-national curve predictions.

#### 3. Plotting, Replication, & Graphics (`Rcodes_paper/`)
* **Plots_mortality.R**: Produces the descriptive demographic functional data visualizations (rainbow plots) for the United States datasets.
* **Figure_1.R to Figure_6.R**: Dedicated scripts containing the exact R code required to generate plots 1 through 6 featured in the paper.

#### 4. Datasets, Replicated Results, & Apps
* **datasets/ & dataset_entries/**: Specialized project directories containing all raw sub-national age-specific mortality datasets and regional names utilized in the case study.
* **Results_Figure_2 to Results_Figure_6 / Test_Results**: Target output directories housing the raw simulation outputs, generated metrics, and test results mapped directly to the paper's corresponding figures.
* **Shiny_app/**: A dedicated standalone folder containing all of the execution scripts and asset code needed to run the interactive web application presented in the paper.

## How to Cite
If you use this code or methodology in your research, please cite our paper:

**APA Style:**
> Jiménez-Varón, C. F., Sun, Y., & Shang, H. L. (2024). Forecasting High-Dimensional Functional Time Series: Application to Sub-National Age-Specific Mortality. *Journal of Computational and Graphical Statistics*, 33(4), 1160–1174. https://doi.org/10.1080/10618600.2024.2319166

**BibTeX:**
```bibtex
@article{jimenezvaron2024forecasting,
  author    = {Jim{\'e}nez-Var{\'o}n, Cristian F. and Sun, Ying and Shang, Han Lin},
  title     = {Forecasting High-Dimensional Functional Time Series: Application to Sub-National Age-Specific Mortality},
  journal   = {Journal of Computational and Graphical Statistics},
  volume    = {33},
  number    = {4},
  pages     = {1160--1174},
  year      = {2024},
  doi       = {10.1080/10618600.2024.2319166},
  url       = {[https://doi.org/10.1080/10618600.2024.2319166](https://doi.org/10.1080/10618600.2024.2319166)}
}
