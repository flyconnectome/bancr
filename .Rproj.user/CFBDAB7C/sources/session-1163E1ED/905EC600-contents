.onLoad <- function(libname, pkgname) {

  # make banc4<->banc3 bridging registrations available
  #register_banc3to4()
  invisible()
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    'Do `choose_banc()` to use many fafbseg::flywire_* functions!\n',
    'Use dr_banc() to get a report on your installation.\n',
    'Trouble? Visit https://flyconnectome.github.io/bancr/SUPPORT.html or #code on banc Slack')
}
