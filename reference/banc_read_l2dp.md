# Read L2 skeleton or dotprops for BANC neurons using fafbseg-py

`banc_read_l2skel` reads one or more neurons as simplified L2 skeletons.

`banc_read_l2dp` reads one or more neurons as simplified dotprops
format. See
[`read_l2skel`](https://rdrr.io/pkg/fafbseg/man/read_l2skel.html).

## Usage

``` r
banc_read_l2dp(id, OmitFailures = TRUE, dataset = NULL, ...)

banc_read_l2skel(id, OmitFailures = TRUE, dataset = NULL, ...)
```

## Arguments

- id:

  One or more flywire ids

- OmitFailures:

  Whether or not to drop neurons that cannot be read from the results
  (rather than erroring out). Default `TRUE`.

- dataset:

  An optional CAVE dataset name (expert use only, by default will choose
  the standard banc dataset). See details.

- ...:

  Additional arguments passed to the `fafbseg.flywire.l2_skeleton` or
  `fafbseg.flywire.l2_dotprops`functions.

## Value

a [`neuronlist`](https://rdrr.io/pkg/nat/man/neuronlist.html) containing
one or more [`neuron`](https://rdrr.io/pkg/nat/man/neuron.html) or
[`dotprops`](https://rdrr.io/pkg/nat/man/dotprops.html) objects. Note
that neurons will be calibrated in `nm` while dotprops will be
calibrated in microns.

## Details

`banc_read_l2dp` uses a special data structure for rapid download of the
dotprops version of neurons required for NBLASTing. It leverages the
python navis / fafbseg-py packages and you will need to install these,
typically using the
[`simple_python`](https://rdrr.io/pkg/fafbseg/man/simple_python.html)
function.

`banc_read_l2skel` treats the dataset argument a little differently than
`banc_read_l2dp` because it actually needs to identify two data sources
a CAVE data

See [`read_l2skel`](https://rdrr.io/pkg/fafbseg/man/read_l2skel.html)
for additional details of

## Examples

``` r
if (FALSE) { # \dontrun{
# one time install of necessary python packages
fafbseg::simple_python(pkgs="fafbseg")

dna02=c("720575941478275714", "720575941512946243")
dna02.latest=banc_latestid(dna02)
dna02.dps <- banc_read_l2dp(dna02.latest)

# plot those
nclear3d()
plot3d(dna02.dps, lwd=3)
# nb dotprops are always in microns
wire3d(banc.surf/1e3, col='grey')

nclear3d()
dna02.skel <- banc_read_l2skel(dna02.latest)
plot3d(dna02.skel, lwd=2)
# nb neuron skeletons are in nm
wire3d(banc.surf, col='grey')
} # }
```
