.onLoad <- function(libname, pkgname) {

  # make banc4<->banc3 bridging registrations available
  #register_banc3to4()
  invisible()
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    'Use `with_banc()` to wrap many additional fafbseg::flywire_* functions for use with the BANC\n',
    'Alternatively `choose_banc()` to set all flywire_* functions to target the BANC!\n',
    'Use dr_banc() to get a report on your installation.\n',
    'Trouble? Visit https://flyconnectome.github.io/bancr/SUPPORT.html or #code on banc Slack')
}
