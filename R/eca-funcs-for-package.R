



#' find path for the "pseudo2" executable, or compile it with admb if it does not exist
#'
#' it seems that admb wants to do everything in the current working directory
Pseudo2BinaryPath <- function() {

	Should.Be <- file.path(".", "pseudo2")  # this should be the name of the executable
	if(file.exists(Should.Be)) {  # if pseudo2 exists in the current working directory, return it.
		return(Should.Be)  # note, I could do some more checking to make sure it is executable, etc, but I won't worry about it
	}
	else {
		warning("Using admb to compile pseudo2.  This should happen only once per session per working directory.", call. = FALSE, immediate. = TRUE, domain = NULL)
		if(file.exists("./pseudo2.tpl")) warning("*** pseudo2.tpl already exists in current working directory.  will not overwrite, but maybe it should be overwritten?\n\n", call. = FALSE, immediate. = TRUE)
		TPL.file <- file.path(system.file(package="psiblings"), "admbfiles", "pseudo2.tpl") # get the pseudo2.tpl file path
		if(!file.exists(TPL.file)) stop(paste("Something is wrong with your psiblings installation.  Can't find file", TPL.file))
		file.copy(TPL.file, "pseudo2.tpl", overwrite=F)
		system2("admb", args = "pseudo2")
		return(Should.Be)  # could do more checking to make sure it was properly compiled.
	}	
}