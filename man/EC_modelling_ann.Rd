% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EC_modelling_ann.R
\name{EC_modelling_ann}
\alias{EC_modelling_ann}
\title{Function to run an Artificial Neural Network (ANN) within EcoCommons. It is a
machine learning model that uses a hidden layer from linear combinations of
predictor variables to predict species occurrence probabilities}
\usage{
EC_modelling_ann(a, response_info, predictor_info, dataset_info)
}
\arguments{
\item{a}{List created from a json file, containing source_file$params}

\item{response_info}{Response object; a nested named list created on EC_build_response}

\item{predictor_info}{Predictor object; a nested named list created on EC_build_predictor}

\item{dataset_info}{Dataset object; a nested named list created on EC_build_dataset}
}
\description{
Function to run an Artificial Neural Network (ANN) within EcoCommons. It is a
machine learning model that uses a hidden layer from linear combinations of
predictor variables to predict species occurrence probabilities
}
