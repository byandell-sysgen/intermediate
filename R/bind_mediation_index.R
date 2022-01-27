#' Bind mediation objects
#' 
#' @export
bind_mediation_index <- function(object, id_name = "mediator_id", ...) {
  # Remake one element as mediation_index object.
  out <- object[[1]]
  
  object <- purrr::transpose(object)
  # Bind elements that are data frames.
  for(i in c("test","fit","fitsLR","best","joint")) {
    out[[i]] <-
      dplyr::bind_rows(
        object[[i]],
        .id = id_name)
  }
  # Bind elements that are vectors; set to NULL if empty.
  for(i in c("normF","driver_names","driver_levels")) {
    isnt <- !sapply(object[[i]], is.null)
    out[[i]] <- {
      if(any(isnt)) {
        if(i == "normF")
          as.data.frame(t(as.data.frame(object[[i]][isnt])))
        else
          object[[i]][isnt]
      }
      else
        NULL
    }
  }

  out
}