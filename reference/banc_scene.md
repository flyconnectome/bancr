# Return a sample Neuroglancer scene URL for BANC dataset

Return a sample Neuroglancer scene URL for BANC dataset

## Usage

``` r
banc_scene(
  ids = NULL,
  open = FALSE,
  layer = NULL,
  url =
    "https://spelunker.cave-explorer.org/#!middleauth+https://global.daf-apis.com/nglstate/api/v1/6431332029693952"
)
```

## Arguments

- ids:

  A set of root ids to include in the scene. Can also be a data.frame.

- open:

  Whether to open the URL in your browser (see
  [`browseURL`](https://rdrr.io/r/utils/browseURL.html))

- layer:

  the segmentation layer for which `ids` intended. Defaults to
  'segmentation proofreading', but could point to another dataset layer.

- url:

  a spelunker neuroglancer URL.

## Value

A character vector containing a single Neuroglancer URL (invisibly when
`open=TRUE`).

## Details

See [banc
slack](https://banc-reconstruction.slack.com/archives/C01RZP5JH9C/p1616522511001900)
for details.

## See also

[`bancsee`](https://flyconnectome.github.io/bancr/reference/bancsee.md)

## Examples

``` r
if (FALSE) { # \dontrun{
browseURL(banc_scene())
banc_scene(open=T)
banc_scene("720575941545083784", open=T)
} # }
```
