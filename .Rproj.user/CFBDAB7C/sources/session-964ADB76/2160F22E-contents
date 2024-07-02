# bancr

<!-- badges: start -->
[![natverse](https://img.shields.io/badge/natverse-Part%20of%20the%20natverse-a241b6)](https://natverse.github.io)
[![Docs](https://img.shields.io/badge/docs-100%25-brightgreen.svg)](https://flyconnectome.github.io/bancr/reference/)
[![R-CMD-check](https://github.com/flyconnectome/banc/workflows/R-CMD-check/badge.svg)](https://github.com/flyconnectome/banc/actions)
[![R-CMD-check](https://github.com/flyconnectome/bancr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/flyconnectome/bancr/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of **bancr** is to support analysis of the Brain And
Nerve Cord dataset aka (BANC), especially autosegmentation data. Those 
data are made available by the banc project led by Wei-Chung Allen Lee (Harvard) and  collaborators including Zetta.ai. 

To access banc resources, you must have permissions to access the [banc
autosegmentation
dataset](https://banc-reconstruction.slack.com/archives/C01RZP5JH9C/p1616522511001900)
and have [confirmed your
acceptance](https://banc-reconstruction.slack.com/archives/C01RZP5JH9C/p1617404290005300)
of the banc proofreading and data ownership guidelines. At this point you should
have a linked Google account that will be authorised (see below) for access to
banc online resources.

Broadly speaking the **bancr** package is a thin wrapper over the 
[fafbseg](https://github.com/natverse/fafbseg) package setting up necessary 
default paths etc. It is based on anothe wrapper for a separate project, 
[fancr](https://github.com/flyconnectome/fancr).

## Installation

You can install the development version of bancr from github:

```r
remotes::install_github('flyconnectome/bancr')
```

To do anything useful with the bancr package, you need authorisation to access
banc resources. To prove your authorisation for programmatic access you must
generate and store a token in your web browser after logging in to an approved
Google account. This should be streamlined by running the following command in R
(which will also set you up for Pythonic access via cloudvolume.)

```r
# set up token - will open your browser to generate a new token
banc_set_token()
# if you already have one do 
# banc_set_token("<my token>")
```

To check that everything is set up properly, try:

```r
dr_banc()

banc_xyz2id(cbind(34495, 82783, 1954), rawcoords=TRUE)
svids=banc_leaves("720575941650432785")
head(svids)
```

### Updating

You can just repeat the install instructions, but this ensures
that all dependencies are updated:

```r
remotes::install_github('flyconnectome/bancr')
```
