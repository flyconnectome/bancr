.onLoad <- function(libname, pkgname) {

  # make banc4<->banc3 bridging registrations available
  #register_banc3to4()
  invisible()
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    'Use `with_banc()` to wrap many additional fafbseg::flywire_* functions for use with the BANC
',
    'Alternatively `choose_banc()` to set all flywire_* functions to target the BANC!
',
    'Use dr_banc() to get a report on your installation.
',
    'Use register_banc_coconat(showerror=FALSE) to register bancr with coconatfly.
',
    'Trouble? Visit https://flyconnectome.github.io/bancr/SUPPORT.html or #code on banc Slack')
}
