# Check if Elastix Transform is possible

Verifies if the specified transform file exists.

## Usage

``` r
check_if_possible(file, on_error = "raise")
```

## Arguments

- file:

  Path to the Elastix transform file.

- on_error:

  Action to take on error: "raise" to stop execution, or any other value
  to return an error message.

## Value

NULL if file exists, or an error message if the file doesn't exist and
on_error is not "raise".
