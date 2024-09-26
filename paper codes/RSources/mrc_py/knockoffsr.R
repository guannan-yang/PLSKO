#' Setup knockoffsr
#'
#' This function initializes Julia and the Knockoffs.jl package.
#' This will install Julia and the required packages
#' if they are missing.
#'
#' @param julia_dir string, path to folder containing Julia executable
#' @param pkg_check logical, check if Knockoffs.jl package exist and install if necessary
#' @param ... Parameters are passed down to JuliaCall::julia_setup
#'
#' @examples
#'
#' knockoffsr::knockoff_setup()
#'
#' @export
knockoff_setup <- function(julia_dir="", pkg_check=TRUE, ...){
  if(julia_dir==""){
    julia <- JuliaCall::julia_setup(installJulia=TRUE,...)
  } else {
    julia <- JuliaCall::julia_setup(JULIA_HOME = julia_dir, ...)
  }
  if(pkg_check) JuliaCall::julia_install_package_if_needed("Knockoffs")
  JuliaCall::julia_library("Knockoffs")

  functions <- JuliaCall::julia_eval("filter(isascii, replace.(string.(propertynames(Knockoffs)),\"!\"=>\"_bang\"))")
  ko <- julia_pkg_import("Knockoffs",functions)
  ko
}

# from: https://github.com/SciML/diffeqr/blob/master/R/diffeqr.R
julia_function <- function(func_name, pkg_name = "Main",
                           env = emptyenv()){
  fname <- paste0(pkg_name, ".", func_name)
  force(fname)
  f <- function(...,
                need_return = c("R", "Julia", "None"),
                show_value = FALSE){
    if (!isTRUE(env$initialized)) {
      env$setup()
    }
    JuliaCall::julia_do.call(func_name = fname, list(...),
                             need_return = match.arg(need_return),
                             show_value = show_value)
  }
  force(f)
  env[[func_name]] <- f
}

julia_pkg_import <- function(pkg_name, func_list){
  env <- new.env(parent = emptyenv())
  env$setup <- function(...){
    JuliaCall::julia_setup(...)
    JuliaCall::julia_library(pkg_name)
    env$initialized <- TRUE
  }
  for (fname in func_list) {
    julia_function(func_name = fname,
                   pkg_name = pkg_name,
                   env = env)
  }
  env
}

