.onAttach = function(libname, pkgname){
  packageStartupMessage("This is serac, a package for ShortlivEd RAdionuclide Chronology of recent sediment cores.")
  packageStartupMessage("   If you need a little help to get started, type in help_serac() or vignette('serac')")
  packageStartupMessage("   If you're using serac for the first time, you may want to use the function user_infos() to enter metadata.")
  packageStartupMessage("   We recommend you re-install the package before each use to get all recent updates, using:")
  packageStartupMessage("      devtools::install_github('rosalieb/serac', build_vignettes = TRUE)")
  packageStartupMessage("   Citation: Bruel and Sabatier (2020) - 10.1016/j.jenvrad.2020.106449")
}
