function(generate_version_h)
    # in case Git is not available, we default to "unknown"
    set(GIT_HASH "unknown")

    # find Git and if available set GIT_HASH variable
    find_package(Git QUIET)
    if(GIT_FOUND)
        execute_process(
                COMMAND ${GIT_EXECUTABLE} --no-pager show -s --pretty=format:%h -n 1
                OUTPUT_VARIABLE GIT_HASH
                OUTPUT_STRIP_TRAILING_WHITESPACE
                WORKING_DIRECTORY
                ${PROJECT_SOURCE_DIR}
                ERROR_QUIET
        )
    endif()

    # generate file version.h based on version.h.in
    configure_file(
            cmake/version.h.in
            generated/version.h
            @ONLY
    )
endfunction()
