# Determining Ice Particles Characteristics in Clouds (project finalized in 2016)

Each ice particle has a unique shape, and this makes it very difficult to determine ice particle characteristics 
(e.g. area and mass) from particle size.

Here, we derived ice particle mass and area from their size using a nonlinear regression model (2nd-order polynomial fit
in log-log space). We, then, classified ice particles based on features such as temperature and cloud type, and derived mass and area for each group.

The steps include data gathering, data cleaning (imputing missing values, correcting wrong values, processing initial data, 
validating initial data, ...), extensive data analysis, classification of ice particles based on their features, implementing
nonlinear regression model, calculating accuracy score (such as R-squared), quantifying uncertainties, and validating
the regression model.

The results were presented in multiple national/international conferences and a paper was published in a peer-reviewed
journal (see the paper: Erfani and Mitchell, 2016).

This model has been used by scientists from the National Center for Atmospheric Research (NCAR) in order to improve global
climate models. They implemented our results into a global climate model and therefore increased the accuracy in the
simulation of ice particles by 20% (for more details, see the paper: Eidhammer et al., 2017).

Selected MATLAB codes, plots, and papers are uploaded. More details will be provided upon request.
