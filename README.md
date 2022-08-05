# Determining Ice Particles Characteristics in Sky Clouds (project finalized in 2016)

### Author:
This package is developed by Ehsan Erfani as part of a chapter of my PhD dissertation and a peer-reviewed paper. 

### Project Description:
Each ice particle in sky clouds has a unique shape, and this makes it very difficult to determine ice particle characteristics 
(e.g. area and mass) from particle size.

Here, we derived ice particle mass and area from their size using a nonlinear regression model (2nd-order polynomial fit
in log-log space). We, then, classified ice particles based on features such as temperature and cloud type, and derived mass and area for each group.

The steps include data gathering (data from instruments on an airplane in multiple field campaigns), data cleaning (imputing missing values, correcting wrong values, processing initial data, validating initial data, ...), extensive data analysis, classification of ice particles based on their features, implementing
nonlinear regression model, calculating accuracy score (such as R-squared), quantifying uncertainties, and validating
the regression model.

The results were presented in multiple national/international conferences and a paper was published in a peer-reviewed
journal (see the paper: Erfani and Mitchell, 2016).

I also collaborated with scientists from the National Center for Atmospheric Research (NCAR) to use my model in order to improve global climate models. We implemented our results into a global climate model and therefore increased the accuracy in the simulation of ice particles by 20% (for more details, see the paper: Eidhammer et al., 2017).

Selected MATLAB codes, plots, and papers are uploaded. More details will be provided upon request.

### The analyses are featured in:
Erfani, E. and Mitchell, D. L.: Developing and bounding ice particle mass- and area-dimension expressions for use in atmospheric models and remote sensing, Atmos. Chem. Phys., 16, 2016.

Eidhammer, T., Morrison, H., Mitchell, D., Gettelman, A., & Erfani, E. (2017). Improvements in global climate model microphysics using a consistent representation of ice particle properties. Journal of Climate, 30(2), 609-629.
