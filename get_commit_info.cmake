if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
    # The most common situation is that code is being compiled from a git repository.
    # This assumes access to the git command and .git folder.
    find_program(GIT_EXECUTABLE git)
    if(GIT_EXECUTABLE)
        execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_COMMIT
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    endif()

    message(STATUS "Git commit info found using git command: ${GIT_COMMIT}")
elseif(EXISTS "${CMAKE_SOURCE_DIR}/git_info.txt")
    # git_info.txt is generated with the script write_git_info.sh
    # The intended use case for this is when the code is planned to be deployed
    # in a system which has no access to the command git and/or the .git folder.
    file(READ "${CMAKE_SOURCE_DIR}/git_info.txt" GIT_INFO_CONTENT)
	string(REGEX MATCH "GIT_COMMIT=\"([^\"]+)\"" _ "${GIT_INFO_CONTENT}")
	set(GIT_COMMIT "${CMAKE_MATCH_1}")

	message(STATUS "Git commit info found in git_info.txt: ${GIT_COMMIT}")
else()
    set(GIT_COMMIT "UNKNOWN")

	message(STATUS "Git commit info not found, setting as ${GIT_COMMIT}")
endif()