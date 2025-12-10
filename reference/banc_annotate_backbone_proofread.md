# Annotate positions as backbone proofread

Mark specific positions as backbone proofread in the CAVE annotation
system.

## Usage

``` r
banc_annotate_backbone_proofread(
  positions,
  user_id,
  units = c("raw", "nm"),
  proofread = TRUE,
  datastack_name = NULL
)
```

## Arguments

- positions:

  3D coordinates in BANC space

- user_id:

  Integer user ID for the annotation

- units:

  Character, coordinate units - either "raw" or "nm"

- proofread:

  Logical, whether to mark as proofread (default TRUE)

- datastack_name:

  Optional datastack name

## Value

Data frame of added annotations

## Examples

``` r
if (FALSE) { # \dontrun{
# Add an annotation to a point in raw voxel space
banc_annotate_backbone_proofread(c(117105, 240526, 5122), user_id = 355, units = "raw")

# Add an annotation to a point in nm
banc_annotate_backbone_proofread(c(468420, 962104 ,230490), user_id = 355, units = "nm")

# deannotate a point, only from points added with given user_id. Use user_id = NULL to remove from full pool
banc_deannotate(c(468420, 962104, 230490), user_id = 355, units = "nm", table = "backbone_proofread")
} # }
```
