% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EC_modelling_fda.R
\name{modelling_fda}
\alias{modelling_fda}
\title{Function to run a Flexible Discriminant Analysis (FDA) within EcoCommons.
It is an statistic regression, a classification model based on a mixture
of linear regression models, which uses optimal scoring to transform the
response variable so that the data are in a better form for linear separation,
and multiple adaptive regression splines to generate the discriminant surface}
\usage{
modelling_fda(a, response_info, predictor_info, dataset_info)
}
\arguments{
\item{a}{List created from a json file, containing source_file$params}

\item{response_info}{Response object; a nested named list created on EC_build_response}

\item{predictor_info}{Predictor object; a nested named list created on EC_build_predictor}

\item{dataset_info}{Dataset object; a nested named list created on EC_build_dataset}
}
\description{
Function to run a Flexible Discriminant Analysis (FDA) within EcoCommons.
It is an statistic regression, a classification model based on a mixture
of linear regression models, which uses optimal scoring to transform the
response variable so that the data are in a better form for linear separation,
and multiple adaptive regression splines to generate the discriminant surface
}
