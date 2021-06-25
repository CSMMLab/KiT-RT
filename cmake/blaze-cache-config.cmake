set(flag 1)
if (WIN32)
    execute_process(COMMAND wmic cpu get L3CacheSize
        OUTPUT_VARIABLE tmp
        RESULT_VARIABLE flag
        ERROR_QUIET)
    if (flag)
        execute_process(COMMAND wmic cpu get L2CacheSize
            OUTPUT_VARIABLE tmp
            RESULT_VARIABLE flag
            ERROR_QUIET)
    endif (flag)
    if (flag)
        execute_process(COMMAND wmic cpu get L1CacheSize
            OUTPUT_VARIABLE tmp
            RESULT_VARIABLE flag
            ERROR_QUIET)
    endif (flag)
endif (WIN32)

if (UNIX)
    execute_process(COMMAND cat /sys/devices/system/cpu/cpu0/cache/index3/size
        OUTPUT_VARIABLE tmp
        RESULT_VARIABLE flag
        ERROR_QUIET)
    if (flag)
        execute_process(COMMAND cat /sys/devices/system/cpu/cpu0/cache/index2/size
            OUTPUT_VARIABLE tmp
            RESULT_VARIABLE flag
            ERROR_QUIET)
    endif (flag)
    if (flag)
        execute_process(COMMAND cat /sys/devices/system/cpu/cpu0/cache/index1/size
            OUTPUT_VARIABLE tmp
            RESULT_VARIABLE flag
            ERROR_QUIET)
    endif (flag)
endif (UNIX)

if (APPLE)
    execute_process(COMMAND sysctl -n hw.l3cachesize
        OUTPUT_VARIABLE tmp
        RESULT_VARIABLE flag
        ERROR_QUIET)
    if (flag)
        execute_process(COMMAND sysctl -n hw.l2cachesize
            OUTPUT_VARIABLE tmp
            RESULT_VARIABLE flag
            ERROR_QUIET)
    endif (flag)
    if (flag)
        execute_process(COMMAND sysctl -n hw.l1icachesize
            OUTPUT_VARIABLE tmp
            RESULT_VARIABLE flag
            ERROR_QUIET)
    endif (flag)

    if (flag EQUAL 0)
        math(EXPR tmp ${tmp}/1024)  # If successful convert to kibibytes to comply with rest
    endif (flag EQUAL 0)
endif (APPLE)

if (flag)
    message(WARNING "Cache size not found automatically. Using default value as cache size.")
    set(tmp "3072")
endif (flag)

string( REGEX MATCH "([0-9][0-9]+)" tmp ${tmp} )
math( EXPR BLAZE_CACHE_SIZE ${tmp}*1024 )
add_compile_definitions( BLAZE_CACHE_SIZE=${BLAZE_CACHE_SIZE}UL )
message( STATUS "Blaze: Automatic Cache Size Configuration = ${BLAZE_CACHE_SIZE} KiB" )
