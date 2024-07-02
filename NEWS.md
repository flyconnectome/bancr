# bancr 0.4.0

* Update spine URL for `banc4to3()` to cope with server move (#5, @jefferis)
* banc surf mesh by @dokato in https://github.com/flyconnectome/bancr/pull/4
* banc to manc registration added by @dokato in https://github.com/flyconnectome/bancr/pull/3

## New Contributors
* @dokato made their first contribution in https://github.com/flyconnectome/bancr/pull/4

**Full Changelog**: https://github.com/flyconnectome/bancr/compare/v0.3.0...v0.4.0

# bancr 0.3.0

* switch to CAVE infrastructure rather than Zetta. See https://global.daf-apis.com/info/ (#1)
* includes update return format for `banc_change_log()`

You may find that you need to generate a new token with `banc_set_token()`

# bancr 0.2.0

* Fast `banc_xyz2id()` mapping using Itanna (spine) services (now thousands per second rather than 10s)
* Added `banc_change_log()`
* Added `dr_banc()` function
* Added a `NEWS.md` file to track changes to the package.
